"""L5 measurement: per-cell intensities in three compartments, with QC exclusions."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def measure_cells(marker: np.ndarray, cells: np.ndarray, nuclei: np.ndarray,
                  saturation: np.ndarray, offset: float = 0.0,
                  label_name: str = "", group: str = "") -> pd.DataFrame:
    """One row per cell.

    Compartments
      nuc  : the DAPI-segmented nucleus
      cell : the whole territory grown from that nucleus
      cyto : cell minus nucleus -- the ring. For a cytoplasmic adaptor like MyD88
             this is the compartment of interest, and it is also immune to any
             DAPI bleed-through into the marker channel.

    `offset` is subtracted per pixel before integrating, so both the mean
    (concentration-like) and the integrated (total-amount) columns are consistent.
    """
    ids = np.unique(cells)
    ids = ids[ids > 0]
    if len(ids) == 0:
        return pd.DataFrame()

    corrected = marker - offset
    cyto = np.where(nuclei > 0, 0, cells)

    def stats(lab, prefix):
        idx = ids
        area = np.asarray(ndi.sum(np.ones_like(lab, dtype=np.float32), lab, idx))
        total = np.asarray(ndi.sum(corrected, lab, idx))
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(area > 0, total / np.maximum(area, 1), np.nan)
        return {f"{prefix}_area": area, f"{prefix}_mean": mean, f"{prefix}_integrated": total}

    df = pd.DataFrame({"label": ids})
    for lab, prefix in [(nuclei, "nuc"), (cells, "cell"), (cyto, "cyto")]:
        df = df.assign(**stats(lab, prefix))

    # QC flags -------------------------------------------------------------
    sat_px = np.asarray(ndi.sum(saturation.astype(np.float32), cells, ids))
    df["saturated_px"] = sat_px
    df["flag_saturated"] = sat_px > 0

    a = df["cell_area"].to_numpy()
    lo, hi = np.percentile(a, [1, 99])
    df["flag_area_outlier"] = (a < lo) | (a > hi)

    # An absolute floor, not a fraction of cell area: what matters is having enough
    # pixels for a stable ring mean. A fractional rule preferentially discards cells
    # in crowded regions, which is exactly the population you must not bias against.
    df["flag_no_cyto"] = df["cyto_area"] < 2000
    df["keep"] = ~(df["flag_saturated"] | df["flag_area_outlier"] | df["flag_no_cyto"])

    df.insert(0, "group", group)
    df.insert(1, "image", label_name)
    return df


def summarise(df: pd.DataFrame, column: str = "cyto_mean") -> pd.Series:
    k = df[df["keep"]]
    return pd.Series({
        "n_cells_total": len(df),
        "n_cells_kept": len(k),
        "mean": k[column].mean(),
        "median": k[column].median(),
        "sd": k[column].std(),
        "sem": k[column].std() / np.sqrt(max(len(k), 1)),
    })


def offset_sensitivity(marker_c: np.ndarray, cells_c: np.ndarray, nuc_c: np.ndarray, sat_c: np.ndarray,
                       marker_t: np.ndarray, cells_t: np.ndarray, nuc_t: np.ndarray, sat_t: np.ndarray,
                       offsets, column: str = "cyto_mean") -> pd.DataFrame:
    """How much does the treatment/control ratio depend on the offset we cannot measure?

    Without a secondary-only control the additive offset is a free parameter, and a
    ratio of two backgrounds-subtracted means is very sensitive to it. Reporting the
    sweep is more honest than reporting one number.
    """
    rows = []
    for off in offsets:
        c = measure_cells(marker_c, cells_c, nuc_c, sat_c, offset=off)
        t = measure_cells(marker_t, cells_t, nuc_t, sat_t, offset=off)
        cm = c[c["keep"]][column].mean()
        tm = t[t["keep"]][column].mean()
        rows.append({"offset": off, "control": cm, "treatment": tm, "ratio_T_over_C": tm / cm})
    return pd.DataFrame(rows)
