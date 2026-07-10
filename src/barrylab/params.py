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

# --- Head direction (default_vars.m: rm.binSizeDir / rm.smthKernDir) ------
#: Directional bin width in degrees (6 -> 60 bins around the circle).
DIR_BIN_DEG = 360.0 / 60.0
#: Default circular boxcar smoothing width in bins for tuning curves.
DIR_SMOOTH_BINS = 5

# --- Adaptive binning (skaggsAdaptiveBinning.m) --------------------------
#: Alpha used by Skaggs adaptive smoothing prior to spatial information.
ADAPTIVE_SMOOTH_ALPHA = 200.0

# --- Place x direction ML model (pxd.m / wrapper_pxd.m) ------------------
#: Spatial bin size in cm for the pxd 3-D histogram. Coarser than RATEMAP_BIN_CM
#: because pxd.m's own docstring warns spatial bin count should stay <= ~30x30
#: (too many spatial bins destabilises the joint place/direction ML fit).
PXD_POS_BIN_CM = 4.0
#: Boxcar smoothing width in bins for the pxd positional (place) marginal.
PXD_POS_SMOOTH_BINS = 5
#: Circular boxcar smoothing width in bins for the pxd directional marginal.
PXD_DIR_SMOOTH_BINS = 3
#: Max EM-style iterations before giving up (pxd.m default).
PXD_MAX_ITER = 1000
#: Fractional log-likelihood change to declare convergence (pxd.m default).
PXD_ACCURACY = 1e-6
#: Floor for expected firing rate/denominators to avoid log(0)/divide-by-zero
#: (pxd.m default).
PXD_TOL = 0.1

#: Value DACQ writes into the position stream for an untracked LED.
UNTRACKED_LED_FLAG = 1023
