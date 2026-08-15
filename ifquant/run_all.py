"""Run L1-L5 across every field in a folder and try to infer the group structure.

    python -m ifquant.run_all

The control/treatment assignment for this dataset has never been supplied, so this
script does not assume one. It measures every field independently and then reports
the evidence for a grouping (acquisition timing, auto-exposure level, cell density,
measured intensity) so a human can confirm it.
"""
from __future__ import annotations

import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ifquant.core import load_linear, exposure_time, estimate_flatfield, apply_flatfield
from ifquant.segment import (segment_nuclei_classical, segment_cells_seeded,
                              drop_edge_cells, check_alignment)
from ifquant.measure import measure_cells

SRC = Path(r"C:\workspace\bio img process\HMC3 activated marker IF GMy88 R INOX")
OUT = Path(r"C:\workspace\bio img process\results")

# Each field is a consecutive triplet: blue(DAPI), green(MyD88), red(iNOS).
# Verified by dominant-channel check -- the 2nd of every triplet is green-dominant.
GREEN = [766, 769, 796, 799, 802, 805, 808, 811, 814, 817, 820, 823]
NUC_DIAM = 215.0


def field_paths(green_num: int):
    return SRC / f"PB154{green_num - 1}.JPG", SRC / f"PB154{green_num}.JPG"


def main():
    OUT.mkdir(exist_ok=True)
    log = []
    def p(s):
        print(s, flush=True); log.append(s)

    p("L2: estimating retrospective flat-field from all %d green images" % len(GREEN))
    flat = estimate_flatfield([field_paths(n)[1] for n in GREEN], "G")
    p(f"L2: vignette falloff centre->corner = {100*(1-flat.min()):.0f}%\n")

    rows, frames = [], []
    for n in GREEN:
        t0 = time.time()
        bpath, gpath = field_paths(n)
        blue, _ = load_linear(bpath, "B")
        green, sat = load_linear(gpath, "G")

        shift, _ = check_alignment(blue, green)
        if np.abs(shift).max() > 1:
            blue = ndi.shift(blue, -np.asarray(shift), order=1, mode="nearest")

        green_c = apply_flatfield(green, flat)
        blue_c = apply_flatfield(blue, flat)

        nuclei = segment_nuclei_classical(blue_c, min_diameter=NUC_DIAM * 0.55)
        cells = segment_cells_seeded(green_c, nuclei)
        cells, nuclei = drop_edge_cells(cells, nuclei)

        # Offset is left at 0 here: it is unidentifiable without a secondary-only
        # control, and applying a guessed constant would distort between-field
        # comparison. Field ranking is unaffected by a common offset.
        df = measure_cells(green_c, cells, nuclei, sat, offset=0.0,
                           label_name=gpath.stem, group="unknown")
        frames.append(df)
        k = df[df["keep"]]

        rows.append({
            "field": gpath.stem,
            "exposure_s": exposure_time(gpath),
            "time": str(pd.Timestamp(  # EXIF capture time
                __import__("PIL").Image.open(gpath)._getexif()[36867].replace(":", "-", 2)))[11:],
            "n_nuclei": int(len(np.unique(nuclei)) - 1),
            "n_kept": len(k),
            "pct_saturated_cells": 100 * df["flag_saturated"].mean(),
            "cyto_mean": k["cyto_mean"].mean(),
            "cyto_median": k["cyto_mean"].median(),
            "cell_mean": k["cell_mean"].mean(),
            "raw8_mean": float(np.asarray(__import__("PIL").Image.open(gpath))[..., 1].mean()),
        })
        p(f"  {gpath.stem}  exp={rows[-1]['exposure_s']:.4f}s  "
          f"nuclei={rows[-1]['n_nuclei']:3d}  kept={rows[-1]['n_kept']:3d}  "
          f"cyto_mean={rows[-1]['cyto_mean']:.4f}  ({time.time()-t0:.0f}s)")

    per_field = pd.DataFrame(rows)
    per_cell = pd.concat(frames, ignore_index=True)
    per_field.to_csv(OUT / "per_field_all.csv", index=False)
    per_cell.to_csv(OUT / "per_cell_all.csv", index=False)

    p("\n=== per-field summary (offset = 0, linear radiance / s) ===")
    p(per_field.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    p("\n=== grouping evidence ===")
    p("  no control/treatment assignment was supplied; these are the separable signals\n")

    p("  (a) auto-exposure level -- the camera shortens exposure on brighter samples,")
    p("      so exposure is itself a coarse readout of sample brightness:")
    for e, g in per_field.groupby("exposure_s"):
        p(f"      {e:.4f}s  n={len(g):2d}  fields: {', '.join(g.field.str[-3:])}")
        p(f"                 cyto_mean {g.cyto_mean.mean():.4f} +- {g.cyto_mean.std():.4f}")

    p("\n  (b) acquisition gaps > 3 min (candidate slide/coverslip changes):")
    t = pd.to_datetime(per_field["time"], format="%H:%M:%S")
    gaps = t.diff().dt.total_seconds() / 60
    for i, gp in enumerate(gaps):
        if i > 0 and gp > 3:
            p(f"      {gp:5.1f} min gap before {per_field.field.iloc[i]}")

    p("\n  (c) is cyto_mean bimodal across fields?")
    v = np.sort(per_field["cyto_mean"].to_numpy())
    p(f"      sorted: {np.array2string(v, precision=3)}")
    p(f"      largest consecutive jump: {np.diff(v).max():.4f} "
      f"between {v[np.argmax(np.diff(v))]:.3f} and {v[np.argmax(np.diff(v))+1]:.3f}")
    p(f"      spread/median = {(v.max()-v.min())/np.median(v):.2f}")

    figure(per_field, per_cell, OUT / "qc_all_fields.png")
    (OUT / "log_all_fields.txt").write_text("\n".join(log), encoding="utf-8")
    p(f"\nwrote {OUT}")


def figure(per_field, per_cell, path):
    fig, ax = plt.subplots(2, 2, figsize=(15, 10))
    order = per_field.sort_values("cyto_mean")
    colors = ["tab:red" if e < 0.7 else "tab:blue" for e in order.exposure_s]
    ax[0, 0].barh(order.field.str[-3:], order.cyto_mean, color=colors)
    ax[0, 0].set_xlabel("cyto_mean (linear/s)")
    ax[0, 0].set_title("per-field MyD88, sorted\nred = 0.625s exposure, blue = 0.769s")

    data = [per_cell[(per_cell.image == f) & per_cell.keep].cyto_mean.to_numpy()
            for f in per_field.field]
    ax[0, 1].violinplot(data, showmeans=True)
    ax[0, 1].set_xticks(range(1, len(per_field) + 1))
    ax[0, 1].set_xticklabels(per_field.field.str[-3:], rotation=45)
    ax[0, 1].set_ylabel("cyto_mean"); ax[0, 1].set_title("per-cell distribution by field")
    ax[0, 1].grid(alpha=.3)

    ax[1, 0].scatter(per_field.exposure_s, per_field.cyto_mean, c=colors, s=90)
    for _, r in per_field.iterrows():
        ax[1, 0].annotate(r.field[-3:], (r.exposure_s, r.cyto_mean),
                          fontsize=8, xytext=(5, 3), textcoords="offset points")
    ax[1, 0].set_xlabel("EXIF exposure (s)"); ax[1, 0].set_ylabel("cyto_mean")
    ax[1, 0].set_title("exposure vs measured intensity"); ax[1, 0].grid(alpha=.3)

    ax[1, 1].scatter(per_field.n_nuclei, per_field.cyto_mean, c=colors, s=90)
    for _, r in per_field.iterrows():
        ax[1, 1].annotate(r.field[-3:], (r.n_nuclei, r.cyto_mean),
                          fontsize=8, xytext=(5, 3), textcoords="offset points")
    ax[1, 1].set_xlabel("nuclei per field"); ax[1, 1].set_ylabel("cyto_mean")
    ax[1, 1].set_title("cell density vs intensity"); ax[1, 1].grid(alpha=.3)

    plt.tight_layout(); plt.savefig(path, dpi=100); plt.close()


if __name__ == "__main__":
    main()
