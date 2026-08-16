"""Run L1-L5 on one control/treatment pair and write per-cell CSV + QC figure.

    python -m ifquant.run_pair --backend classical
    python -m ifquant.run_pair --backend cellpose
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ifquant.core import (load_linear, exposure_time, iso_speed, exif_audit,
                          estimate_flatfield, apply_flatfield, offset_from_darkest)
from ifquant.segment import (segment_nuclei_classical, segment_nuclei_cellpose,
                              segment_cells_seeded, segment_cells_cellpose,
                              drop_edge_cells, check_alignment)
from ifquant.measure import measure_cells, summarise, offset_sensitivity

SRC = Path(r"C:\workspace\bio img process\HMC3 activated marker IF GMy88 R INOX")
OUT = Path(r"C:\workspace\bio img process\results")
GREEN_SET = ["PB154766", "PB154769", "PB154796", "PB154799", "PB154802", "PB154805",
             "PB154808", "PB154811", "PB154814", "PB154817", "PB154820", "PB154823"]
PAIR = {"control": ("PB154765", "PB154766"), "treatment": ("PB154768", "PB154769")}
NUC_DIAM = 215.0


def process(group, blue_stem, green_stem, flat, backend, log):
    t0 = time.time()
    blue, _ = load_linear(SRC / f"{blue_stem}.JPG", "B")
    green, sat = load_linear(SRC / f"{green_stem}.JPG", "G")
    bp, gp = SRC / f"{blue_stem}.JPG", SRC / f"{green_stem}.JPG"
    log(f"[{group}] {blue_stem}(B) exp={exposure_time(bp):.4f}s iso={iso_speed(bp):.0f}  "
        f"{green_stem}(G) exp={exposure_time(gp):.4f}s iso={iso_speed(gp):.0f}")

    # The two channels are separate exposures with a filter-cube change in between,
    # so the field can shift. Register blue onto green before using blue's masks
    # to sample green.
    shift, err = check_alignment(blue, green)
    log(f"[{group}] channel misalignment: dy={shift[0]:+.1f} dx={shift[1]:+.1f} px -> corrected")
    if np.abs(shift).max() > 1:
        blue = ndi.shift(blue, -np.asarray(shift), order=1, mode="nearest")

    green_c = apply_flatfield(green, flat)          # L2
    blue_c = apply_flatfield(blue, flat)

    if backend == "cellpose":
        nuclei = segment_nuclei_cellpose(blue_c, diameter=NUC_DIAM)
        cells = segment_cells_cellpose(green_c, blue_c, diameter=NUC_DIAM * 1.6)
    else:
        nuclei = segment_nuclei_classical(blue_c, min_diameter=NUC_DIAM * 0.55)
        cells = segment_cells_seeded(green_c, nuclei)
    cells, nuclei = drop_edge_cells(cells, nuclei)
    log(f"[{group}] segmented: {len(np.unique(nuclei))-1} nuclei "
        f"(after edge removal), {time.time()-t0:.0f}s")

    return dict(group=group, blue=blue_c, green=green_c, sat=sat,
                nuclei=nuclei, cells=cells, stem=green_stem)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="classical", choices=["classical", "cellpose"])
    ap.add_argument("--column", default="cyto_mean")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    lines = []
    def log(s):
        print(s, flush=True)
        lines.append(s)

    log(f"=== backend={args.backend}  metric={args.column} ===")

    # L0: the camera was on P + ISO AUTO, so verify what it actually chose per frame
    # before comparing anything across frames.
    audit = pd.DataFrame(exif_audit([SRC / f"{s}.JPG" for s in GREEN_SET]))
    log(f"L0: exposure {audit.exposure_s.min():.4f}-{audit.exposure_s.max():.4f}s, "
        f"ISO {audit.iso.min():.0f}-{audit.iso.max():.0f}, "
        f"WhiteBalance={sorted(set(audit.white_balance.dropna()))} (0=auto, 1=manual)")
    if audit.white_balance.eq(0).any():
        log("L0: WARNING auto white balance -- per-channel gains vary between frames "
            "and cannot be corrected downstream; treat cross-image G comparisons as "
            "semi-quantitative only.")

    log("L2: estimating retrospective flat-field from %d green images..." % len(GREEN_SET))
    flat = estimate_flatfield([SRC / f"{s}.JPG" for s in GREEN_SET], "G")
    log(f"L2: vignette falloff centre->corner = {100*(1-flat.min()):.0f}%")

    res = {g: process(g, *PAIR[g], flat=flat, backend=args.backend, log=log)
           for g in ("control", "treatment")}

    # L3 offset: no secondary-only control available -> use darkest-percentile
    # estimate and report the full sensitivity sweep.
    off = min(offset_from_darkest(res[g]["green"], 1.0) for g in res)
    log(f"L3: offset estimate (1st pct of flat-fielded field) = {off:.4f} linear/s")

    frames = []
    for g, r in res.items():
        df = measure_cells(r["green"], r["cells"], r["nuclei"], r["sat"],
                           offset=off, label_name=r["stem"], group=g)
        frames.append(df)
    per_cell = pd.concat(frames, ignore_index=True)
    per_cell.to_csv(OUT / f"per_cell_{args.backend}.csv", index=False)

    log("\n=== L5 per-cell summary (offset-subtracted, linear radiance / second) ===")
    tab = []
    for col in ["nuc_mean", "cell_mean", "cyto_mean", "cyto_integrated"]:
        row = {"metric": col}
        for g in ("control", "treatment"):
            s = summarise(per_cell[per_cell.group == g], col)
            row[g] = s["mean"]
            row[f"{g}_n"] = int(s["n_cells_kept"])
        row["ratio_T/C"] = row["treatment"] / row["control"]
        tab.append(row)
    tab = pd.DataFrame(tab)
    log(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    tab.to_csv(OUT / f"summary_{args.backend}.csv", index=False)

    log("\n=== L3 offset sensitivity (the parameter we cannot measure) ===")
    # Physical upper bound: an offset that drives the dimmest control cells negative
    # cannot be right, so sweep only up to that. This is the tightest constraint the
    # data alone can give -- narrowing it further needs a secondary-only control.
    raw_ctrl = measure_cells(res["control"]["green"], res["control"]["cells"],
                             res["control"]["nuclei"], res["control"]["sat"], offset=0.0)
    off_max = float(np.percentile(raw_ctrl[raw_ctrl["keep"]][args.column], 1))
    log(f"L3: offset must be < {off_max:.4f} or the dimmest control cells go negative")
    sweep = offset_sensitivity(
        res["control"]["green"], res["control"]["cells"], res["control"]["nuclei"], res["control"]["sat"],
        res["treatment"]["green"], res["treatment"]["cells"], res["treatment"]["nuclei"], res["treatment"]["sat"],
        offsets=np.linspace(0.0, off_max, 6), column=args.column)
    log(sweep.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    sweep.to_csv(OUT / f"offset_sweep_{args.backend}.csv", index=False)

    log("\n=== L6 statistics ===")
    log("  n here = 1 field per group. Cell-level spread is shown for QC only;")
    log("  a cell-level test would be pseudoreplication. The unit of replication is")
    log("  the coverslip/experiment, so run every field and aggregate before testing.")
    for g in ("control", "treatment"):
        s = summarise(per_cell[per_cell.group == g], args.column)
        log(f"  {g:<10} n_cells={int(s['n_cells_kept']):3d}  "
            f"mean={s['mean']:.4f}  sd={s['sd']:.4f}  sem={s['sem']:.4f}")

    qc_figure(res, per_cell, flat, args, OUT / f"qc_{args.backend}.png", args.column)
    (OUT / f"log_{args.backend}.txt").write_text("\n".join(lines), encoding="utf-8")
    log(f"\nwrote {OUT}")


def qc_figure(res, per_cell, flat, args, path, column):
    from skimage.segmentation import find_boundaries
    fig, ax = plt.subplots(3, 3, figsize=(17, 14))
    for i, g in enumerate(("control", "treatment")):
        r = res[g]
        sl = (slice(800, 2000), slice(800, 2400))
        gr = r["green"][sl]
        ax[i, 0].imshow(r["blue"][sl], cmap="Blues_r"); ax[i, 0].set_title(f"{g}: DAPI (flat-fielded)")
        ax[i, 1].imshow(gr, cmap="Greens_r",
                        vmin=np.percentile(gr, 2), vmax=np.percentile(gr, 98))
        b = find_boundaries(r["cells"][sl], mode="thick")
        n = find_boundaries(r["nuclei"][sl], mode="thick")
        ov = np.zeros(b.shape + (4,)); ov[b] = [1, 0.3, 0, 1]; ov[n] = [0, 0.6, 1, 1]
        ax[i, 1].imshow(ov)
        ax[i, 1].set_title(f"{g}: cell (orange) / nucleus (blue) boundaries")
        d = per_cell[(per_cell.group == g) & per_cell.keep]
        ax[i, 2].hist(d[column], bins=45, color="tab:green" if i else "tab:blue", alpha=.75)
        ax[i, 2].axvline(d[column].mean(), color="k", ls="--")
        ax[i, 2].set_title(f"{g}: {column}  n={len(d)}  mean={d[column].mean():.3f}")
        ax[i, 2].set_xlabel("linear radiance / s")

    ax[2, 0].imshow(flat, cmap="magma"); ax[2, 0].set_title("L2 retrospective flat-field")
    ax[2, 0].set_xticks([]); ax[2, 0].set_yticks([])

    k = per_cell[per_cell.keep]
    data = [k[k.group == g][column].to_numpy() for g in ("control", "treatment")]
    ax[2, 1].violinplot(data, showmeans=True)
    ax[2, 1].set_xticks([1, 2]); ax[2, 1].set_xticklabels(["control", "treatment"])
    ax[2, 1].set_ylabel(column); ax[2, 1].set_title("per-cell distribution"); ax[2, 1].grid(alpha=.3)

    sweep = pd.read_csv(OUT / f"offset_sweep_{args.backend}.csv")
    ax[2, 2].plot(sweep["offset"], sweep["ratio_T_over_C"], "o-")
    ax[2, 2].axhline(1.0, color="r", ls=":")
    ax[2, 2].set_xlabel("assumed additive offset"); ax[2, 2].set_ylabel("ratio T/C")
    ax[2, 2].set_title("offset sensitivity of the effect size"); ax[2, 2].grid(alpha=.3)

    for a in ax[:2, :2].ravel(): a.set_xticks([]); a.set_yticks([])
    plt.tight_layout(); plt.savefig(path, dpi=100); plt.close()


if __name__ == "__main__":
    main()
