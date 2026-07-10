"""Format-agnostic analyses that consume a Recording (or plain arrays)."""

from .grid import grid_properties, grid_score, gridness, spatial_autocorrelogram
from .maps import occupancy_map, rate_map, smooth_ratemap, spike_map
from .spatial import adaptive_binning, spatial_information

__all__ = [
    "rate_map",
    "occupancy_map",
    "spike_map",
    "smooth_ratemap",
    "spatial_information",
    "adaptive_binning",
    "spatial_autocorrelogram",
    "gridness",
    "grid_properties",
    "grid_score",
]
