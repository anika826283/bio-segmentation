"""The full study: 6 conditions x 2 technical replicates, morphology-led.

    python -m ifquant.run_study

Design constraints this script is built around (see ifquant/groups.py):
  - biological n = 1 per condition, so NO statistical test is run
  - rep1/rep2 are technical replicates, used as a REPRODUCIBILITY CHECK:
    a metric is only reportable if the two replicates agree on the ranking
  - the green channel's absolute intensity is exposure-confounded and is
    included only to demonstrate that it fails the reproducibility check
  - the red channel was shot at 1.000 s for every field, so it is clean
"""
from __future__ import annotations

import argparse, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ifquant.core import load_linear, exposure_time, estimate_flatfield, apply_flatfield
from ifquant.segment import (segment_nuclei_classical, segment_cells_seeded,
                              drop_edge_cells, check_alignment)
from ifquant.morphology import measure_field, dapi_bleed_check
from ifquant.groups import SRC, CONDITIONS, all_fields

SRC = Path(SRC)
OUT = Path(r"C:\workspace\bio img process\results")
NUC_DIAM = 215.0

# Metrics carried into the summary, with how far each can be trusted.
METRICS = [
    ("blue_nuc_mean",       "control",   "DAPI per nucleus -- should be CONSTANT"),
    ("n_nuclei",            "count",     "cells per field"),
    ("nuc_area",            "shape",     "nucleus size"),
    ("nuc_circularity",     "shape",     "nucleus roundness (1.0 = circle)"),
    ("nuc_solidity",        "shape",     "nucleus convexity"),
    ("nuc_aspect",          "shape",     "nucleus elongation"),
    ("cell_area",           "shape*",    "cell territory size (segmentation-dependent)"),
    ("cell_nuc_area_ratio", "shape*",    "cytoplasm-to-nucleus size"),
    ("red_cyto_mean",       "intensity", "iNOS, cytoplasm -- CLEAN channel"),
    ("red_cell_cv",         "intensity", "iNOS punctateness (scale-invariant)"),
    ("red_nc_ratio",        "intensity", "iNOS nuclear:cytoplasmic"),
    ("green_nc_ratio",      "intensity", "MyD88 nuclear:cytoplasmic (exposure-immune)"),
    ("green_cell_cv",       "intensity", "MyD88 punctateness (scale-invariant)"),
    ("green_cyto_mean",     "CONFOUNDED", "MyD88 absolute -- exposure-confounded"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-exposure", action="store_true",
                    help="treat every image as if it had the same exposure, i.e. skip "
                         "dividing by EXIF exposure time")
    args = ap.parse_args()
    corr = not args.no_exposure
    tag = "" if corr else "_noexp"

    OUT.mkdir(exist_ok=True)
    log = []
    def p(s):
        print(s, flush=True); log.append(s)

    p(f"=== exposure correction: {'ON (divide by EXIF exposure)' if corr else 'OFF (assume fixed exposure)'} ===")
    fields = all_fields()
    p("L2: per-channel retrospective flat-fields")
    flat = {}
    for ch, idx in [("B", 2), ("G", 3), ("R", 4)]:
        paths = [SRC / f"PB154{f[idx]}.JPG" for f in fields]
        flat[ch] = estimate_flatfield(paths, ch)
        p(f"    {ch}: vignette falloff = {100*(1-flat[ch].min()):.0f}%")

    frames, rows = [], []
    for cond, rep, bn, gn, rn in fields:
        t0 = time.time()
        bp, gp, rp = (SRC / f"PB154{n}.JPG" for n in (bn, gn, rn))
        blue, sat_b = load_linear(bp, "B", normalise_exposure=corr)
        green, sat_g = load_linear(gp, "G", normalise_exposure=corr)
        red, sat_r = load_linear(rp, "R", normalise_exposure=corr)

        shift, _ = check_alignment(blue, green)
        if np.abs(shift).max() > 1:
            blue = ndi.shift(blue, -np.asarray(shift), order=1, mode="nearest")
            sat_b = ndi.shift(sat_b.astype(np.float32), -np.asarray(shift), order=0) > 0.5

        blue = apply_flatfield(blue, flat["B"])
        green = apply_flatfield(green, flat["G"])
        red = apply_flatfield(red, flat["R"])

        nuclei = segment_nuclei_classical(blue, min_diameter=NUC_DIAM * 0.55)
        cells = segment_cells_seeded(green, nuclei)
        cells, nuclei = drop_edge_cells(cells, nuclei)

        df = measure_field(nuclei, cells,
                           {"blue": blue, "green": green, "red": red},
                           {"green": sat_g, "red": sat_r},
                           cond, rep, gp.stem)
        frames.append(df)

        ks, ki = df[df.keep_shape], df[df.keep_intensity]
        row = {"condition": cond, "rep": rep, "field": gp.stem,
               "exp_B": exposure_time(bp), "exp_G": exposure_time(gp),
               "exp_R": exposure_time(rp),
               "n_nuclei": int(len(np.unique(nuclei)) - 1),
               "n_shape": len(ks), "n_intensity": len(ki),
               "pct_saturated": 100 * df.flag_saturated.mean(),
               "dapi_green_corr": dapi_bleed_check(df)}
        for m, tier, _ in METRICS:
            if m == "n_nuclei":
                continue
            row[m] = (ks if tier.startswith("shape") else ki)[m].median()
        rows.append(row)
        p(f"  {cond:<17} rep{rep}  {gp.stem}  n={row['n_nuclei']:3d}  ({time.time()-t0:.0f}s)")

    per_cell = pd.concat(frames, ignore_index=True)
    per_field = pd.DataFrame(rows)
    per_cell.to_csv(OUT / f"study_per_cell{tag}.csv", index=False)
    per_field.to_csv(OUT / f"study_per_field{tag}.csv", index=False)

    p("\n=== DAPI bleed-through check (per-cell corr of nuclear DAPI vs nuclear green) ===")
    p(f"  median across fields = {per_field.dapi_green_corr.median():+.3f}")
    p("  high positive -> the green nuclear signal is largely DAPI leakage,")
    p("  which would inflate green_nc_ratio. Interpret that metric accordingly.")

    p("\n=== reproducibility: does rep2 agree with rep1? ===")
    p("  a metric is only reportable if the two technical replicates rank the")
    p("  six conditions the same way\n")
    p(f"  {'metric':<22}{'tier':<12}{'Spearman':>9}{'  verdict'}")
    cons = []
    for m, tier, desc in METRICS:
        piv = per_field.pivot(index="condition", columns="rep", values=m).reindex(CONDITIONS)
        rho = spearmanr(piv[1], piv[2]).statistic if piv.notna().all().all() else np.nan
        verdict = ("reproducible" if rho >= 0.7 else
                   "weak" if rho >= 0.3 else "NOT reproducible")
        cons.append({"metric": m, "tier": tier, "spearman": rho,
                     "verdict": verdict, "description": desc})
        p(f"  {m:<22}{tier:<12}{rho:>9.3f}  {verdict}")
    cons = pd.DataFrame(cons)
    cons.to_csv(OUT / f"study_reproducibility{tag}.csv", index=False)

    p("\n=== signal-to-noise: can any metric separate the conditions at all? ===")
    p("  Spearman over 6 conditions is far too noisy to lean on (rho must exceed")
    p("  ~0.83 to be significant at n=6), so measure it directly instead:")
    p("    technical SD  = field-to-field scatter, from the rep1-rep2 differences")
    p("    condition SD  = scatter of the six condition means")
    p("    SNR           = condition SD / technical SD; must exceed ~1 to discriminate")
    p("    fields needed = how many fields per well to resolve the observed spread\n")
    p(f"  {'metric':<22}{'tech SD%':>10}{'cond SD%':>10}{'SNR':>7}{'fields':>8}  verdict")
    snr_rows = []
    for m, tier, desc in METRICS:
        piv = per_field.pivot(index="condition", columns="rep", values=m).reindex(CONDITIONS)
        if piv.isna().any().any():
            continue
        mean_of = piv.mean(axis=1)
        grand = mean_of.mean()
        tech_sd = float(np.std(piv[1] - piv[2], ddof=1) / np.sqrt(2))
        cond_sd = float(mean_of.std(ddof=1))
        snr = cond_sd / tech_sd if tech_sd > 0 else np.nan
        # fields per well so that SEM of a condition mean is 1/3 of the spread
        need = int(np.ceil((3 * tech_sd / cond_sd) ** 2)) if cond_sd > 0 else 999
        verdict = ("can discriminate" if snr >= 1.5 else
                   "marginal" if snr >= 1.0 else "swamped by field noise")
        snr_rows.append({"metric": m, "tier": tier, "tech_sd_pct": 100*tech_sd/grand,
                         "cond_sd_pct": 100*cond_sd/grand, "snr": snr,
                         "fields_needed": need, "verdict": verdict, "description": desc})
        p(f"  {m:<22}{100*tech_sd/grand:>9.1f}%{100*cond_sd/grand:>9.1f}%"
          f"{snr:>7.2f}{need:>8d}  {verdict}")
    snr_df = pd.DataFrame(snr_rows)
    snr_df.to_csv(OUT / f"study_snr{tag}.csv", index=False)
    p(f"\n  -> median fields per well required: {int(snr_df.fields_needed.median())}")
    p("     (currently 2). This is the cheapest fix available: image more fields")
    p("     per well on the existing slides -- no new biology needed.")

    p("\n=== condition means, metrics that passed (Spearman >= 0.7) ===")
    good = cons[cons.spearman >= 0.7].metric.tolist()
    if not good:
        p("  none passed")
    for m in good:
        piv = per_field.pivot(index="condition", columns="rep", values=m).reindex(CONDITIONS)
        piv["mean"] = piv.mean(axis=1)
        p(f"\n  {m}  ({cons[cons.metric==m].description.iloc[0]})")
        for c in CONDITIONS:
            p(f"    {c:<18} rep1={piv.loc[c,1]:9.4f}  rep2={piv.loc[c,2]:9.4f}  "
              f"mean={piv.loc[c,'mean']:9.4f}")

    p("\n=== the confounded metric, shown as a negative control ===")
    piv = per_field.pivot(index="condition", columns="rep", values="green_cyto_mean").reindex(CONDITIONS)
    ex = per_field.pivot(index="condition", columns="rep", values="exp_G").reindex(CONDITIONS)
    for c in CONDITIONS:
        p(f"  {c:<18} rep1={piv.loc[c,1]:.4f} (exp {ex.loc[c,1]:.4f}s)   "
          f"rep2={piv.loc[c,2]:.4f} (exp {ex.loc[c,2]:.4f}s)")
    p("  health and disease swapped exposures between reps; watch the direction flip.")

    p("\n=== statistics ===")
    p("  NOT RUN. Biological n = 1 per condition (one well each, one session).")
    p("  rep1/rep2 are technical replicates; testing across them would be")
    p("  pseudoreplication. Everything above is descriptive.")

    figure(per_field, per_cell, cons, OUT / f"study_qc{tag}.png")
    (OUT / f"study_log{tag}.txt").write_text("\n".join(log), encoding="utf-8")
    p(f"\nwrote {OUT}")


def figure(per_field, per_cell, cons, path):
    show = [m for m in ["nuc_area", "nuc_circularity", "cell_nuc_area_ratio",
                        "red_cyto_mean", "red_cell_cv", "green_cyto_mean"]]
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    x = np.arange(len(CONDITIONS))
    for a, m in zip(ax.ravel(), show):
        piv = per_field.pivot(index="condition", columns="rep", values=m).reindex(CONDITIONS)
        a.plot(x, piv[1], "o-", label="rep1")
        a.plot(x, piv[2], "s-", label="rep2")
        rho = cons[cons.metric == m].spearman.iloc[0]
        tier = cons[cons.metric == m].tier.iloc[0]
        a.set_title(f"{m}\n{tier}, rep-agreement rho={rho:+.2f}",
                    color="tab:red" if tier == "CONFOUNDED" else "black")
        a.set_xticks(x)
        a.set_xticklabels([c.replace("control1_", "").replace("control2_", "")
                           .replace("treatment", "t") for c in CONDITIONS], rotation=30)
        a.grid(alpha=.3); a.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(path, dpi=100); plt.close()


if __name__ == "__main__":
    main()
