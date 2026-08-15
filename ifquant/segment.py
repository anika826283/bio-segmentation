"""L4 segmentation: nuclei from DAPI, whole cells from DAPI + marker channel."""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu, sobel
from skimage.measure import label
from skimage.segmentation import watershed, clear_border
from skimage.registration import phase_cross_correlation


def check_alignment(a: np.ndarray, b: np.ndarray, ds: int = 8):
    """Sub-field shift between the two channel exposures (filter-cube change / drift)."""
    def small(x):
        H, W = x.shape
        return x[: H // ds * ds, : W // ds * ds].reshape(H // ds, ds, W // ds, ds).mean(axis=(1, 3))
    shift, error, _ = phase_cross_correlation(small(a), small(b), upsample_factor=10)
    return shift * ds, float(error)


# ------------------------------------------------------------------- nuclei
def segment_nuclei_classical(blue: np.ndarray, min_diameter: float = 110.0,
                             smooth: float = 12.0, ds: int = 4):
    """Distance-transform seeded watershed.

    Runs on a downsampled copy (nuclei are ~215 px across, so 4x costs no accuracy)
    and the labels are then upsampled. `min_diameter` sets the minimum separation
    between seeds and is the knob that trades merged clumps against over-splitting.
    """
    b = blue[: blue.shape[0] // ds * ds, : blue.shape[1] // ds * ds]
    b = b.reshape(b.shape[0] // ds, ds, b.shape[1] // ds, ds).mean(axis=(1, 3))

    sm = ndi.gaussian_filter(b, smooth / ds)
    # DAPI background here is bright and the field is 31% nuclei, so flatten the
    # slow-varying component before thresholding rather than using raw Otsu.
    detrended = sm - ndi.gaussian_filter(sm, 300 / ds)
    mask = detrended > threshold_otsu(detrended)
    mask = ndi.binary_fill_holes(mask)
    mask = ndi.binary_opening(mask, np.ones((3, 3)))

    dist = ndi.distance_transform_edt(mask)
    coords = peak_local_max(ndi.gaussian_filter(dist, 2.0), labels=mask,
                            min_distance=int(min_diameter / ds / 2), exclude_border=False)
    seeds = np.zeros(mask.shape, dtype=np.int32)
    seeds[tuple(coords.T)] = np.arange(1, len(coords) + 1)
    seeds = ndi.grey_dilation(seeds, size=(3, 3))

    lab = watershed(-dist, seeds, mask=mask)
    lab = ndi.zoom(lab, (blue.shape[0] / lab.shape[0], blue.shape[1] / lab.shape[1]), order=0)
    return _fit(lab, blue.shape)


def segment_nuclei_cellpose(blue: np.ndarray, diameter: float = 215.0):
    from cellpose import models
    m = models.CellposeModel(gpu=False, model_type="nuclei")
    img = _to_uint8(blue)
    masks, _, _ = m.eval(img, diameter=diameter, channels=[0, 0])
    return masks.astype(np.int32)


# -------------------------------------------------------------------- cells
def segment_cells_seeded(green: np.ndarray, nuclei: np.ndarray, ds: int = 4,
                         max_expand: float = 80.0):
    """Grow each nucleus into its territory along the marker-channel gradient.

    Better than Voronoi (which ignores the image entirely): the watershed lines
    settle on green-channel edges, i.e. real cell borders, wherever they exist.

    `max_expand` (full-resolution px) caps how far a territory may reach beyond its
    own nucleus. Without it, a seed sitting in a smooth region keeps flooding until
    it hits another seed's basin and can swallow a large part of the frame, which
    both corrupts that cell's mean and starves its neighbours of cytoplasm.
    """
    def down(x, order=0):
        H, W = x.shape
        x = x[: H // ds * ds, : W // ds * ds]
        if order == 0:
            return x[::ds, ::ds]
        return x.reshape(H // ds, ds, W // ds, ds).mean(axis=(1, 3))

    g = down(green, order=1)
    seeds = down(nuclei, order=0)

    reach = ndi.distance_transform_edt(seeds == 0) <= (max_expand / ds)
    mask = reach | (seeds > 0)

    elevation = sobel(ndi.gaussian_filter(g, 3.0))
    lab = watershed(elevation, seeds, mask=mask)
    lab = ndi.zoom(lab, (green.shape[0] / lab.shape[0], green.shape[1] / lab.shape[1]), order=0)
    return _fit(lab, green.shape)


def segment_cells_cellpose(green: np.ndarray, blue: np.ndarray, diameter: float = 340.0):
    from cellpose import models
    m = models.CellposeModel(gpu=False, model_type="cyto3")
    rgb = np.zeros(green.shape + (3,), dtype=np.uint8)
    rgb[..., 1] = _to_uint8(green)     # cytoplasm channel
    rgb[..., 2] = _to_uint8(blue)      # nucleus channel
    masks, _, _ = m.eval(rgb, diameter=diameter, channels=[2, 3])
    return masks.astype(np.int32)


# ------------------------------------------------------------------- helpers
def drop_edge_cells(cells: np.ndarray, nuclei: np.ndarray):
    """Remove cells whose nucleus touches the frame -- their territory is truncated."""
    keep = np.unique(clear_border(nuclei))
    keep = keep[keep > 0]
    out_c = np.where(np.isin(cells, keep), cells, 0)
    out_n = np.where(np.isin(nuclei, keep), nuclei, 0)
    return out_c, out_n


def _to_uint8(x: np.ndarray) -> np.uint8:
    lo, hi = np.percentile(x, [0.5, 99.5])
    return np.clip((x - lo) / max(hi - lo, 1e-9) * 255, 0, 255).astype(np.uint8)


def _fit(lab: np.ndarray, shape) -> np.ndarray:
    if lab.shape != shape:
        pad = [(0, max(0, shape[i] - lab.shape[i])) for i in range(2)]
        lab = np.pad(lab, pad, mode="edge")[: shape[0], : shape[1]]
    return lab.astype(np.int32)
