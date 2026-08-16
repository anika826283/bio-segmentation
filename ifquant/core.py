"""
IF quantification pipeline for RGB-JPEG widefield images (Olympus E-M5 II on scope).

L1 linearise   : inverse sRGB EOTF, then divide by EXIF exposure time and ISO gain
L2 shading     : retrospective flat-field DIVISION (multiplicative correction)
L3 offset      : additive offset subtraction (camera black + non-specific)
L4 segment     : nuclei from DAPI, whole cell from DAPI+marker
L5 measure     : per-cell nucleus / whole-cell / cytoplasmic-ring intensities

Everything downstream of L1 works in "linear photon rate" units, NOT 8-bit JPEG values.
"""
from __future__ import annotations

import numpy as np
from PIL import Image, ExifTags
from pathlib import Path
from scipy import ndimage as ndi

Image.MAX_IMAGE_PIXELS = None

CH = {"R": 0, "G": 1, "B": 2}

# The camera was left in P mode with ISO AUTO (confirmed from the Super Control
# Panel), so ISO is a free variable in principle, not a constant of the setup.
# Sensor response is linear in ISO gain, so ISO is divided out alongside exposure
# time. Normalising to the ISO the reference batch actually landed on (rather than
# to ISO 100) keeps the units identical to the numbers already reported in
# README.md -- with an all-ISO-1600 batch this factor is exactly 1.
ISO_REF = 1600.0


# ---------------------------------------------------------------- L1 linearise
def _exif(path: Path) -> dict:
    ex = Image.open(path)._getexif() or {}
    return {ExifTags.TAGS.get(k, k): v for k, v in ex.items()}


def exposure_time(path: Path) -> float:
    """EXIF exposure time in seconds."""
    return float(_exif(path)["ExposureTime"])


def iso_speed(path: Path) -> float:
    """EXIF ISO sensitivity. Pillow renames this tag across versions, and some
    writers store it as a tuple, so accept both spellings and unwrap sequences."""
    tags = _exif(path)
    for key in ("ISOSpeedRatings", "PhotographicSensitivity"):
        if key in tags:
            v = tags[key]
            return float(v[0] if isinstance(v, (tuple, list)) else v)
    raise KeyError(f"no ISO tag in EXIF of {path.name}")


def exif_audit(paths) -> "list[dict]":
    """Per-image exposure/ISO/WB/program, for checking that the batch really was
    shot under one set of camera settings.

    Worth running before trusting any cross-image comparison: P mode + ISO AUTO
    means the camera is free to change exposure AND gain between frames, and
    auto WB would additionally change the per-channel gains, which no amount of
    downstream normalisation can undo.
    """
    rows = []
    for p in map(Path, paths):
        tags = _exif(p)
        try:
            iso = iso_speed(p)
        except KeyError:
            iso = float("nan")
        rows.append({
            "file": p.name,
            "exposure_s": float(tags.get("ExposureTime", float("nan"))),
            "iso": iso,
            # EXIF spec: WhiteBalance 0 = auto, 1 = manual/preset.
            "white_balance": tags.get("WhiteBalance"),
            "exposure_program": tags.get("ExposureProgram"),
        })
    return rows


def srgb_to_linear(u8: np.ndarray) -> np.ndarray:
    """Inverse sRGB EOTF. 8-bit display values -> linear radiance in [0, 1].

    Approximation: the camera applies its own picture-mode tone curve, which is
    close to but not exactly sRGB. This is the main residual error in the rescue
    route and the reason results are semi-quantitative.
    """
    x = u8.astype(np.float32) / 255.0
    return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4).astype(np.float32)


def load_linear(path: Path, channel: str, normalise_exposure: bool = True):
    """Return (linear_rate, saturation_mask). rate is linear radiance per second
    at ISO_REF, so images shot at different auto-exposure settings are comparable."""
    raw = np.asarray(Image.open(path))[..., CH[channel]]
    sat = raw >= 254
    lin = srgb_to_linear(raw)
    if normalise_exposure:
        lin = lin / exposure_time(path) / (iso_speed(path) / ISO_REF)
    return lin, sat


# ------------------------------------------------------------------ L2 shading
def estimate_flatfield(paths, channel: str, ds: int = 8, smooth_sigma: float = 500.0):
    """Retrospective illumination estimate from a set of images of the same setup.

    Shading is multiplicative, so this is estimated in LINEAR space and applied by
    division. Each image is normalised by its own median first so that field-to-field
    brightness differences do not leak into the shading estimate; the pixel-wise
    median across images then cancels the (position-independent) sample content.

    Assumes fields were chosen without positional bias -- true for randomly picked
    fields, false if e.g. the operator always centred on a bright cluster.

    `smooth_sigma` is in FULL-RESOLUTION pixels and must be far larger than a cell
    (~350 px here) or residual cell texture leaks into the estimate and gets divided
    out of the data. Vignetting is a slow optical effect: err on the side of too smooth.
    """
    acc = []
    for p in paths:
        lin, _ = load_linear(p, channel)
        small = lin[: lin.shape[0] // ds * ds, : lin.shape[1] // ds * ds]
        small = small.reshape(small.shape[0] // ds, ds, small.shape[1] // ds, ds).mean(axis=(1, 3))
        acc.append(small / np.median(small))
    flat = np.median(np.stack(acc), axis=0)
    flat = ndi.gaussian_filter(flat, smooth_sigma / ds)
    return (flat / flat.max()).astype(np.float32)


def apply_flatfield(img: np.ndarray, flat_small: np.ndarray) -> np.ndarray:
    zoom = (img.shape[0] / flat_small.shape[0], img.shape[1] / flat_small.shape[1])
    flat = ndi.zoom(flat_small, zoom, order=1)
    flat = flat[: img.shape[0], : img.shape[1]]
    if flat.shape != img.shape:                      # guard against off-by-one from zoom
        pad = [(0, img.shape[i] - flat.shape[i]) for i in range(2)]
        flat = np.pad(flat, pad, mode="edge")
    return img / np.clip(flat, 1e-6, None)


# ------------------------------------------------------------------- L3 offset
def offset_from_negative_control(path: Path, channel: str, flat_small=None) -> float:
    """The correct offset: mean of a no-primary / secondary-only image, same settings."""
    lin, _ = load_linear(path, channel)
    if flat_small is not None:
        lin = apply_flatfield(lin, flat_small)
    return float(np.mean(lin))


def offset_from_darkest(img: np.ndarray, pct: float = 1.0) -> float:
    """Fallback when no negative control exists: low percentile of the corrected field.

    On a confluent monolayer these pixels are thin inter-cell gaps, not true
    cell-free background, so this OVER-estimates the offset and therefore
    UNDER-estimates the signal. Always report a sensitivity sweep alongside it.
    """
    return float(np.percentile(img, pct))
