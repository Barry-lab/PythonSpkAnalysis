"""The lab's standard "tint" colour map. Port of ``tintColorMap.m``."""

from matplotlib.colors import ListedColormap

# Literal 5x3 RGB triplet array from tintColorMap.m -- "comparable to tint"
# (the lab's Tint spike-sorting/ratemap-viewing software) per the MATLAB
# source's own comment.
_TINT_RGB = [
    [0.00, 0.00, 0.75],
    [0.00, 0.75, 0.90],
    [0.10, 1.00, 0.40],
    [1.00, 1.00, 0.00],
    [0.95, 0.00, 0.00],
]

TINT_CMAP = ListedColormap(_TINT_RGB, name="tint")
