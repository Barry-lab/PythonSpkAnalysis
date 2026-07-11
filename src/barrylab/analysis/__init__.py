"""Format-agnostic analyses that consume a Recording (or plain arrays)."""

from ._stats import nan_corr
from .decoding import (
    decode_calc_post,
    decode_calc_bayes_post,
    decode_crunch_data,
    decode_process_post,
    decode_get_spike_counts,
)
from .direction import (
    circular_mean,
    directional_fft,
    directional_occupancy,
    directionality,
    head_direction_map,
    head_direction_stats,
    led2v1_dir_offset,
    resultant_vector,
)
from .grid import grid_properties, grid_score, gridness, spatial_autocorrelogram, square_grid_props
from .grid_ellipse import fit_grid_ellipse, regularise_elliptic_sac
from .lfp import power_spectrum
from .maps import occupancy_map, rate_map, smooth_ratemap, spike_map
from .phase_precession import circular_linear_regression
from .place_fields import field_stats, spatial_coherence, spatial_correlation
from .pxd import place_by_direction, pxd
from .spatial import adaptive_binning, spatial_information, true_ppm_bounding_box
from .spatial_fft import FFTParams, grid_fft_components, power_spectrum_2d, to_polar_fft

__all__ = [
    "rate_map",
    "occupancy_map",
    "spike_map",
    "smooth_ratemap",
    "spatial_information",
    "adaptive_binning",
    "true_ppm_bounding_box",
    "spatial_autocorrelogram",
    "gridness",
    "grid_properties",
    "grid_score",
    "square_grid_props",
    "head_direction_map",
    "head_direction_stats",
    "directional_occupancy",
    "circular_mean",
    "resultant_vector",
    "directionality",
    "directional_fft",
    "led2v1_dir_offset",
    "field_stats",
    "spatial_coherence",
    "spatial_correlation",
    "nan_corr",
    "pxd",
    "place_by_direction",
    "decode_calc_post",
    "decode_calc_bayes_post",
    "decode_get_spike_counts",
    "decode_process_post",
    "decode_crunch_data",
    "circular_linear_regression",
    "fit_grid_ellipse",
    "regularise_elliptic_sac",
    "power_spectrum",
    "FFTParams",
    "power_spectrum_2d",
    "to_polar_fft",
    "grid_fft_components",
]
