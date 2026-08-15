"""Single source of truth for the study design.

Confirmed by the user on 2026-08-15. Each field is a consecutive triplet
blue(DAPI) / green(MyD88) / red(iNOS); the number below is the GREEN file.

Design reality, which drives everything downstream:
  - 6 conditions, ONE well each -> biological n = 1 per condition
  - rep1 and rep2 are TECHNICAL replicates (same slide, same 20-minute
    session, 17:11-17:32), not independent experiments
  - treatment1..4 are four DIFFERENT treatments, not four replicates
"""

SRC = r"C:\workspace\bio img process\HMC3 activated marker IF GMy88 R INOX"

CONDITIONS = ["control1_health", "control2_disease",
              "treatment1", "treatment2", "treatment3", "treatment4"]

# green file number per (condition, replicate)
FIELDS = {
    ("control1_health",  1): 766, ("control1_health",  2): 808,
    ("control2_disease", 1): 769, ("control2_disease", 2): 811,
    ("treatment1",       1): 796, ("treatment1",       2): 814,
    ("treatment2",       1): 799, ("treatment2",       2): 817,
    ("treatment3",       1): 802, ("treatment3",       2): 820,
    ("treatment4",       1): 805, ("treatment4",       2): 823,
}

# Exposure per channel across the 12 fields, from EXIF:
#   BLUE   0.125 - 0.200 s   (60% spread)  -> confounded
#   GREEN  0.625 - 0.769 s   (23% spread)  -> confounded
#   RED    1.000 s exactly   (0% spread)   -> CLEAN
# Auto-exposure hit its 1 s ceiling on the dim red channel and could not
# normalise it, so RED is the only channel whose absolute intensity is
# comparable between fields.
CLEAN_CHANNEL = "R"


def triplet(green_num: int):
    """(blue, green, red) file numbers for a field."""
    return green_num - 1, green_num, green_num + 1


def all_fields():
    """[(condition, rep, blue, green, red), ...] in acquisition order."""
    out = []
    for (cond, rep), g in FIELDS.items():
        b, g, r = triplet(g)
        out.append((cond, rep, b, g, r))
    return sorted(out, key=lambda x: x[3])
