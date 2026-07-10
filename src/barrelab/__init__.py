"""barrelab — spatial-coding analysis for tetrode and Neuropixels recordings.

Front-end loaders (DACQ today; Neuropixels/Open Ephys next) populate a common
:class:`~barrelab.core.recording.Recording`, which the analysis functions
consume. Typical use::

    import barrelab as bl

    rec = bl.load_dacq("/data/r1838_10-12-14/83812141")
    rmap = bl.rate_map(rec, electrode=2, cluster=3)
    bits_per_spike, _ = bl.spatial_information(rmap, bl.occupancy_map(rec))
"""

from .core.recording import Electrode, Recording
from .io.axona import load_dacq
from .analysis.direction import (
    circular_mean,
    directional_fft,
    directional_occupancy,
    directionality,
    head_direction_map,
    head_direction_stats,
    resultant_vector,
)
from .analysis.grid import (
    grid_properties,
    grid_score,
    gridness,
    spatial_autocorrelogram,
)
from .analysis.maps import occupancy_map, rate_map, smooth_ratemap, spike_map
from .analysis.spatial import adaptive_binning, spatial_information

__version__ = "0.1.0"

__all__ = [
    "Recording",
    "Electrode",
    "load_dacq",
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
    "head_direction_map",
    "head_direction_stats",
    "directional_occupancy",
    "circular_mean",
    "resultant_vector",
    "directionality",
    "directional_fft",
]
