"""L4b morphology, and the intensity measures that survive auto-exposure.

Why morphology leads here: shape is a geometric quantity, so it is completely
immune to the exposure / gain / gamma problems that make the green channel's
absolute intensity unusable in this dataset.

Trust levels:
  nucleus shape   -- from thresholded DAPI, genuinely biological
  cell shape      -- from the bounded seeded watershed, so it is PARTLY an
                     artefact of `max_expand`. Comparable between fields only
                     because every field used identical parameters; do not read
                     it as a true cell outline.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.measure import regionprops_table

NUC_PROPS = ("label", "area", "perimeter", "equivalent_diameter_area",
             "solidity", "eccentricity", "axis_major_length",
             "axis_minor_length", "feret_diameter_max")


def _shape_table(labels: np.ndarray, prefix: str) -> pd.DataFrame:
    t = pd.DataFrame(regionprops_table(labels, properties=NUC_PROPS))
    t = t.rename(columns={c: f"{prefix}_{c}" for c in t.columns if c != "label"})
    a, p = t[f"{prefix}_area"], t[f"{prefix}_perimeter"]
    # circularity: 1.0 = perfect circle, lower = more irregular / ramified
    t[f"{prefix}_circularity"] = np.where(p > 0, 4 * np.pi * a / np.maximum(p ** 2, 1e-9), np.nan)
    t[f"{prefix}_aspect"] = t[f"{prefix}_axis_major_length"] / np.maximum(
        t[f"{prefix}_axis_minor_length"], 1e-9)
    return t


def _intensity_by_label(img: np.ndarray, labels: np.ndarray, ids: np.ndarray):
    area = np.asarray(ndi.sum(np.ones_like(img, np.float32), labels, ids))
    total = np.asarray(ndi.sum(img, labels, ids))
    sq = np.asarray(ndi.sum(img.astype(np.float64) ** 2, labels, ids))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(area > 0, total / np.maximum(area, 1), np.nan)
        var = np.where(area > 1, sq / np.maximum(area, 1) - mean ** 2, np.nan)
        cv = np.sqrt(np.maximum(var, 0)) / np.maximum(mean, 1e-9)
    return area, total, mean, cv


def measure_field(nuclei: np.ndarray, cells: np.ndarray,
                  channels: dict[str, np.ndarray],
                  saturation: dict[str, np.ndarray],
                  condition: str, rep: int, field: str) -> pd.DataFrame:
    """One row per cell: shape, plus per-channel intensity in three compartments."""
    ids = np.unique(cells)
    ids = ids[ids > 0]
    if len(ids) == 0:
        return pd.DataFrame()

    df = _shape_table(nuclei, "nuc").merge(_shape_table(cells, "cell"), on="label", how="inner")
    df["cell_nuc_area_ratio"] = df["cell_area"] / np.maximum(df["nuc_area"], 1)
    ids = df["label"].to_numpy()

    cyto = np.where(nuclei > 0, 0, cells)
    for name, img in channels.items():
        for lab, comp in [(nuclei, "nuc"), (cells, "cell"), (cyto, "cyto")]:
            area, total, mean, cv = _intensity_by_label(img, lab, ids)
            df[f"{name}_{comp}_mean"] = mean
            df[f"{name}_{comp}_integrated"] = total
            if comp == "cyto":
                df["cyto_area"] = area
            if comp == "cell":
                df[f"{name}_cell_cv"] = cv     # scale-invariant: how punctate the marker is
        # Nuclear-to-cytoplasmic ratio. Both compartments come from the SAME
        # exposure of the SAME image, so exposure, gain and gamma cancel --
        # this is the one intensity measure here that auto-exposure cannot touch.
        df[f"{name}_nc_ratio"] = df[f"{name}_nuc_mean"] / np.maximum(df[f"{name}_cyto_mean"], 1e-9)

    for name, sat in saturation.items():
        df[f"{name}_saturated_px"] = np.asarray(
            ndi.sum(sat.astype(np.float32), cells, ids))

    # QC ------------------------------------------------------------------
    a = df["cell_area"].to_numpy()
    lo, hi = np.percentile(a, [1, 99])
    df["flag_area_outlier"] = (a < lo) | (a > hi)
    df["flag_no_cyto"] = df["cyto_area"] < 2000
    sat_cols = [c for c in df.columns if c.endswith("_saturated_px")]
    df["flag_saturated"] = df[sat_cols].sum(axis=1) > 0
    # Shape metrics stay valid even when the marker clips, so morphology keeps
    # saturated cells; only intensity columns need the stricter filter.
    df["keep_shape"] = ~(df["flag_area_outlier"] | df["flag_no_cyto"])
    df["keep_intensity"] = df["keep_shape"] & ~df["flag_saturated"]

    df.insert(0, "condition", condition)
    df.insert(1, "rep", rep)
    df.insert(2, "field", field)
    return df


def dapi_bleed_check(df: pd.DataFrame) -> float:
    """Per-cell correlation between nuclear DAPI and nuclear green, within a field.

    MyD88 is a cytoplasmic adaptor, so a strong positive correlation here points
    at DAPI emission leaking through the green filter rather than real nuclear
    MyD88 -- which would inflate every green nuclear-to-cytoplasmic ratio.
    """
    k = df[df["keep_intensity"]]
    if len(k) < 10:
        return np.nan
    return float(np.corrcoef(k["blue_nuc_mean"], k["green_nuc_mean"])[0, 1])
