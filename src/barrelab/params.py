"""Default analysis parameters.

Ported from the lab's MATLAB ``default_vars.m``. Only the parameters used by
the currently implemented analyses are kept here; the rest will be added as
their analyses are ported. Values are chosen to reproduce the historic MATLAB
behaviour exactly (see ``validation/``).
"""

# --- Position processing (postprocess_pos_data.m) ------------------------
#: Speeds above this (m/s) are treated as tracking jumps and rejected.
POS_MAX_SPEED_MS = 2.0
#: Boxcar width (seconds) used to smooth position before deriving speed/dir.
POS_SMOOTH_BOXCAR_S = 0.4

# --- Rate maps (easy_make_ratemap.m / default_vars.m) --------------------
#: Spatial bin size in cm. 2 cm is the lab standard (== old Tint 8 px bin).
RATEMAP_BIN_CM = 2.0
#: Default boxcar smoothing kernel width in bins for 2-D rate maps.
RATEMAP_SMOOTH_BINS = 5
#: Default Gaussian smoothing sigma in bins (used when smoothing="gaussian").
RATEMAP_GAUSS_SIGMA = 2.0

# --- Adaptive binning (skaggsAdaptiveBinning.m) --------------------------
#: Alpha used by Skaggs adaptive smoothing prior to spatial information.
ADAPTIVE_SMOOTH_ALPHA = 200.0

#: Value DACQ writes into the position stream for an untracked LED.
UNTRACKED_LED_FLAG = 1023
