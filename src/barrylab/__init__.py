"""barrylab — spatial-coding analysis for tetrode and Neuropixels recordings.

Front-end loaders (DACQ today; Neuropixels/Open Ephys next) populate a common
:class:`~barrylab.core.recording.Recording`, which the analysis functions
consume. Typical use::

    import barrylab as bl

    rec = bl.load_dacq("/data/r1838_10-12-14/83812141")
    rmap = bl.rate_map(rec, electrode=2, cluster=3)
    bits_per_spike, _ = bl.spatial_information(rmap, bl.occupancy_map(rec))
"""

from .core.recording import Electrode, Recording
from .io.axona import load_dacq
from .analysis._stats import nan_corr
from .analysis.decoding import (
    bayesian_posterior,
    bayesian_posterior_with_prior,
    crunch_decode,
    decode_maximum_likelihood,
    spike_counts_per_window,
)
from .analysis.direction import (
    circular_mean,
    directional_fft,
    directional_occupancy,
    directionality,
    head_direction_map,
    head_direction_stats,
    led2v1_dir_offset,
    resultant_vector,
)
from .analysis.grid import (
    fast_gridness,
    grid_properties,
    grid_score,
    gridness,
    spatial_autocorrelogram,
)
from .analysis.grid_ellipse import fit_grid_ellipse, regularise_elliptic_sac
from .analysis.lfp import power_spectrum
from .analysis.line_trajectory import best_line, line_residual, reget_line_band
from .analysis.maps import occupancy_map, rate_map, smooth_ratemap, spike_map
from .analysis.phase_precession import circular_linear_regression
from .analysis.place_fields import field_stats, spatial_coherence, spatial_correlation
from .analysis.pxd import place_by_direction, pxd
from .analysis.spatial import adaptive_binning, spatial_information, true_ppm_bounding_box

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
    "true_ppm_bounding_box",
    "spatial_autocorrelogram",
    "gridness",
    "grid_properties",
    "grid_score",
    "fast_gridness",
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
    "bayesian_posterior",
    "bayesian_posterior_with_prior",
    "spike_counts_per_window",
    "decode_maximum_likelihood",
    "crunch_decode",
    "circular_linear_regression",
    "fit_grid_ellipse",
    "regularise_elliptic_sac",
    "power_spectrum",
    "best_line",
    "reget_line_band",
    "line_residual",
]
