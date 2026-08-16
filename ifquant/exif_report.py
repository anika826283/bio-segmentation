"""Report the acquisition settings actually recorded in the image files.

    python -m ifquant.exif_report "HMC3 activated marker IF GMy88 R INOX"

Why this exists: the camera menu shows what the body is set to *now*. The files
record what it was set to *at capture time*. Those disagreed once already (the
body is on the fixed Daylight preset, but standard EXIF WhiteBalance reads 0 =
"auto"), so settings that affect L1-L3 are read back from the files themselves.

Standard EXIF comes from Pillow. The settings that actually decide whether this
pipeline's assumptions hold -- Gradation above all -- live in Olympus's MakerNote,
which Pillow does not parse, so those need exiftool. Without exiftool the standard
half still prints and the MakerNote half is reported as unavailable rather than
guessed.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ifquant.core import exif_audit

# MakerNote fields that change how the pipeline must be interpreted.
MAKERNOTE_FIELDS = {
    "PictureMode":          "base tone curve",
    "Gradation":            "Auto = per-image adaptive curve -> L1 error varies frame to frame",
    "ColorSpace":           "confirms the L1 inverse-EOTF choice",
    "WhiteBalance2":        "the real WB preset (standard EXIF flag is unreliable here)",
    "WhiteBalanceTemperature": "WB colour temperature",
    "Sharpness":            "spatial filtering alters per-pixel statistics",
    "Contrast":             "tone curve slope",
    "Saturation":           "per-channel gain",
    "NoiseFilter":          "denoising alters per-pixel statistics",
    "NoiseReduction":       "denoising alters per-pixel statistics",
    "ShadingCompensation":  "in-camera vignetting correction -- would collide with L2 flat-field",
}


def makernote_audit(paths) -> "pd.DataFrame | None":
    """Olympus MakerNote settings per file, or None when exiftool is unavailable."""
    exe = shutil.which("exiftool")
    if exe is None:
        return None
    # No -n: values stay human-readable ("Auto", "Normal") rather than numeric codes.
    args = [exe, "-json"]
    args += [f"-{f}" for f in MAKERNOTE_FIELDS]
    args += ["-FileName", *[str(p) for p in paths]]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    return pd.DataFrame(json.loads(out))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    src = Path(sys.argv[1])
    paths = sorted(src.glob("*.JPG")) if src.is_dir() else [src]
    if not paths:
        sys.exit(f"no .JPG found in {src}")
    print(f"{len(paths)} file(s) from {src}\n")

    std = pd.DataFrame(exif_audit(paths))
    print("=== standard EXIF ===")
    for col in ("exposure_s", "iso", "white_balance", "light_source", "exposure_program"):
        vals = sorted(set(std[col].dropna()))
        flag = "" if len(vals) == 1 else "   <-- VARIES ACROSS BATCH"
        print(f"  {col:<18} {vals}{flag}")

    mn = makernote_audit(paths)
    print("\n=== Olympus MakerNote ===")
    if mn is None:
        print("  exiftool not found -- install it to read Gradation and the real WB preset:")
        print("    macOS: brew install exiftool     Ubuntu: apt install libimage-exiftool-perl")
        print("  (this is the only way to settle the Gradation question from the files)")
        return
    for field, why in MAKERNOTE_FIELDS.items():
        if field not in mn.columns:
            print(f"  {field:<24} (not recorded)")
            continue
        vals = sorted(set(mn[field].dropna().astype(str)))
        flag = "" if len(vals) == 1 else "   <-- VARIES ACROSS BATCH"
        print(f"  {field:<24} {vals}{flag}")
        print(f"  {'':<24} ^ {why}")

    _warn_adaptive(mn)


def _warn_adaptive(mn: pd.DataFrame):
    """Flag settings that make the camera's transform vary from frame to frame.

    These are the ones that break cross-image comparison, because control and
    treatment fields differ in content by construction, so a content-dependent
    transform lands differently on the two groups.
    """
    def vals(col):
        return set(mn[col].dropna().astype(str)) if col in mn.columns else set()

    if any("auto" in g.lower() for g in vals("Gradation")):
        print("\n  WARNING Gradation is Auto: a different shadow-lifting curve per frame,")
        print("  chosen from that frame's own histogram. L1's error becomes per-image")
        print("  rather than a fixed bias.")

    if any("enhance" in m.lower() for m in vals("PictureMode")):
        print("\n  WARNING Picture Mode is i-Enhance: a scene-adaptive contrast/saturation")
        print("  transform that varies with image content, is not recorded anywhere, and")
        print("  is not undone by the inverse sRGB EOTF. Control and treatment fields")
        print("  differ in content by construction, so this is a systematic confound on")
        print("  the group comparison, not just added noise.")

    if any("auto" in w.lower() for w in vals("WhiteBalance2")):
        print("\n  WARNING White balance is Auto: per-channel R/G/B gains were chosen per")
        print("  frame, so absolute G values are not directly comparable across images.")

    if any(s.lower() not in ("off", "0") for s in vals("ShadingCompensation")):
        print("\n  WARNING Shading Compensation is on: the camera already corrected")
        print("  vignetting, so the L2 retrospective flat-field would correct it twice.")


if __name__ == "__main__":
    main()
