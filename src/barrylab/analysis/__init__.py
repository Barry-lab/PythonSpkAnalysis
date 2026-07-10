"""Format-agnostic analyses that consume a Recording (or plain arrays)."""

from ._stats import nan_corr
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
from .grid import grid_properties, grid_score, gridness, spatial_autocorrelogram
from .maps import occupancy_map, rate_map, smooth_ratemap, spike_map
from .place_fields import field_stats, spatial_coherence, spatial_correlation
from .pxd import place_by_direction, pxd
from .spatial import adaptive_binning, spatial_information, true_ppm_bounding_box

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
]
