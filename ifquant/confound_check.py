"""Is the between-field intensity difference real, or an artefact of auto-exposure?

The camera metered each field and chose an exposure. cyto_mean is computed as
(linearised 8-bit) / exposure. If auto-exposure fully normalised the 8-bit
brightness, then cyto_mean carries no information beyond 1/exposure -- and
anything that changes total field brightness (cell density, focus, a bright
clump in the metering zone) propagates straight into the per-cell number.

This script quantifies how much of the observed difference survives that check.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
OUT = Path(r"C:\workspace\bio img process\results")


def main():
    d = pd.read_csv(OUT / "per_field_all.csv")
    log = []
    def p(s):
        print(s, flush=True); log.append(s)

    lo = d[d.exposure_s < 0.7]      # 0.625 s
    hi = d[d.exposure_s > 0.7]      # 0.769 s

    p("=== 1. do the two exposure levels separate cyto_mean? ===")
    p(f"  0.625s (n={len(lo)}): {sorted(np.round(lo.cyto_mean, 4))}")
    p(f"  0.769s (n={len(hi)}): {sorted(np.round(hi.cyto_mean, 4))}")
    p(f"  min(0.625s group) = {lo.cyto_mean.min():.4f}")
    p(f"  max(0.769s group) = {hi.cyto_mean.max():.4f}")
    p(f"  -> {'PERFECTLY SEPARATED' if lo.cyto_mean.min() > hi.cyto_mean.max() else 'overlapping'}")

    p("\n=== 2. but is that separation just 1/exposure? ===")
    r_meas = lo.cyto_mean.mean() / hi.cyto_mean.mean()
    r_exp = hi.exposure_s.mean() / lo.exposure_s.mean()
    p(f"  measured cyto_mean ratio      = {r_meas:.4f}")
    p(f"  pure exposure ratio           = {r_exp:.4f}")
    p(f"  residual beyond exposure      = {r_meas / r_exp:.4f}   <- 1.00 means NO extra signal")
    p(f"  raw 8-bit mean, 0.625s group  = {lo.raw8_mean.mean():.2f}")
    p(f"  raw 8-bit mean, 0.769s group  = {hi.raw8_mean.mean():.2f}")
    p(f"  -> auto-exposure equalised the 8-bit brightness to within "
      f"{100*abs(lo.raw8_mean.mean()/hi.raw8_mean.mean()-1):.1f}%")

    p("\n=== 3. cell density confound ===")
    p(f"  nuclei/field, 0.625s group = {lo.n_nuclei.mean():.1f}")
    p(f"  nuclei/field, 0.769s group = {hi.n_nuclei.mean():.1f}")
    p(f"  density ratio              = {lo.n_nuclei.mean()/hi.n_nuclei.mean():.4f}")
    p(f"  exposure-implied brightness ratio = {r_exp:.4f}")
    p("  -> a denser field is a brighter field, so auto-exposure shortens the")
    p("     exposure, which inflates every per-cell value in that field.")
    p(f"  corr(n_nuclei, cyto_mean) over all 12 fields = "
      f"{np.corrcoef(d.n_nuclei, d.cyto_mean)[0,1]:+.3f}")
    for name, g in [("0.625s", lo), ("0.769s", hi)]:
        if len(g) > 2:
            p(f"  corr within {name} group (exposure held constant) = "
              f"{np.corrcoef(g.n_nuclei, g.cyto_mean)[0,1]:+.3f}  (n={len(g)})")

    p("\n=== 4. acquisition timeline ===")
    t = pd.to_datetime(d["time"], format="%H:%M:%S")
    gaps = (t.diff().dt.total_seconds() / 60).fillna(0)
    for i, r in d.iterrows():
        mark = "   <-- %.1f min gap" % gaps[i] if gaps[i] > 3 else ""
        p(f"  {r.time}  {r.field}  exp={r.exposure_s:.4f}s  n={r.n_nuclei:3d}{mark}")
    p("  -> 766 and 769 are 21 s apart. They cannot be two different")
    p("     coverslips, so the 'control=766 / treatment=769' assumption")
    p("     inherited from the earlier chat is not supported by the timing.")

    p("\n=== 5. saturation, per field ===")
    p(f"  pct of cells containing a clipped pixel: "
      f"{d.pct_saturated_cells.min():.1f}% - {d.pct_saturated_cells.max():.1f}%")
    p(f"  fields above 40%: {', '.join(d[d.pct_saturated_cells>40].field.str[-3:])}")
    p("  -> these fields lose their brightest cells to the QC filter, which")
    p("     biases their means DOWNWARD by an unknown amount.")

    (OUT / "confound_check.txt").write_text("\n".join(log), encoding="utf-8")
    print(f"\nwrote {OUT/'confound_check.txt'}")


if __name__ == "__main__":
    main()
