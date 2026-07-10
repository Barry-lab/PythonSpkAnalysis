"""The common in-memory representation shared by every front-end loader.

Whatever the acquisition system (Axona/DACQ today; Neuropixels via SpikeGLX or
Open Ephys next), a loader's job is to populate a :class:`Recording`. All
analysis code consumes this object and nothing else, so a rate map does not care
whether the spikes came from a 2010 tetrode or a 2026 Neuropixel.

Kept as a plain, transparent container (arrays + light metadata) so a student
can inspect every field. It intentionally has no file-format knowledge.
"""

from typing import Dict

import numpy as np


class Electrode(object):
    """Sorted spikes from one recording channel-group (a tetrode, or a probe).

    Attributes
    ----------
    identifier : int
        Tetrode number (Axona) or group id.
    spike_times : ndarray (n_spikes,)
        Spike times in seconds.
    pos_index : ndarray (n_spikes,)
        Index into the position arrays of the sample co-occurring with each
        spike (0-based). May be ``None`` if position is absent.
    clusters : ndarray (n_spikes,)
        Cluster id per spike (0 = unsorted/noise), or ``None`` if uncut.
    metadata : dict
        Raw acquisition header, kept for provenance.
    """

    def __init__(self, identifier, spike_times, pos_index=None,
                 clusters=None, metadata=None):
        self.identifier = identifier
        self.spike_times = np.asarray(spike_times)
        self.pos_index = None if pos_index is None else np.asarray(pos_index)
        self.clusters = None if clusters is None else np.asarray(clusters)
        self.metadata = metadata or {}

    def cluster_ids(self):
        """Sorted unique cluster ids present (empty if uncut)."""
        if self.clusters is None:
            return np.array([], dtype=int)
        return np.unique(self.clusters)

    def spike_mask(self, cluster):
        """Boolean mask selecting spikes belonging to ``cluster``.

        Truncated to the length of ``pos_index`` when position is present, so
        the mask always indexes valid position samples (matches the historic
        ``easy_make_ratemap`` behaviour).
        """
        if self.clusters is None:
            raise ValueError("electrode %r has no cluster (cut) data" % (self.identifier,))
        mask = self.clusters == cluster
        if self.pos_index is not None:
            mask = mask[: len(self.pos_index)]
        return mask

    def spike_pos_index(self, cluster):
        """Position-sample indices of spikes from ``cluster`` (0-based)."""
        if self.pos_index is None:
            raise ValueError("electrode %r has no position mapping" % (self.identifier,))
        return self.pos_index[self.spike_mask(cluster)]

    def __repr__(self):
        n_clusters = len(self.cluster_ids())
        return "Electrode(id=%r, n_spikes=%d, n_clusters=%d)" % (
            self.identifier, self.spike_times.size, n_clusters,
        )


class Recording(object):
    """A single recording session: position/behaviour plus sorted spikes.

    Attributes
    ----------
    name : str
        Session name (source file stem).
    xy : ndarray (n_pos, 2)
        Position in pixels, in the tracked-window frame.
    direction : ndarray (n_pos,)
        Head direction in degrees, anticlockwise from the +x axis.
    speed : ndarray (n_pos,)
        Running speed in cm/s.
    pos_times : ndarray (n_pos,)
        Time of each position sample in seconds.
    pos_sample_rate : float
        Position sampling rate in Hz.
    pixels_per_metre : float
        Spatial scale for converting pixels to cm.
    window : tuple (min_x, max_x, min_y, max_y)
        Tracked camera window, in pixels.
    electrodes : dict[int, Electrode]
        Sorted spikes keyed by tetrode/group id.
    lfp : dict[int, dict]
        Local field potential channels keyed by EEG channel number. Each
        value is a dict with ``samples`` (1-D, volts), ``sample_rate`` (Hz),
        ``channel`` and ``kind`` (``"eeg"`` = 250 Hz or ``"egf"`` = 4.8 kHz,
        both undecimated). Empty if no LFP was loaded.
    metadata : dict
        Free-form provenance (source path, format, raw headers).
    """

    def __init__(self, name, xy=None, direction=None, speed=None,
                 pos_times=None, pos_sample_rate=None, pixels_per_metre=None,
                 window=None, electrodes=None, lfp=None, metadata=None):
        self.name = name
        self.xy = xy
        self.direction = direction
        self.speed = speed
        self.pos_times = pos_times
        self.pos_sample_rate = pos_sample_rate
        self.pixels_per_metre = pixels_per_metre
        self.window = window
        self.electrodes = electrodes if electrodes is not None else {}  # type: Dict[int, Electrode]
        self.lfp = lfp if lfp is not None else {}  # type: Dict[int, dict]
        self.metadata = metadata or {}

    @property
    def n_pos(self):
        return 0 if self.xy is None else self.xy.shape[0]

    @property
    def duration(self):
        """Session duration in seconds (from the position clock)."""
        return None if self.pos_times is None else float(self.pos_times[-1])

    @property
    def window_extent(self):
        """(height, width) of the tracked window in pixels (y, x order)."""
        min_x, max_x, min_y, max_y = self.window
        return (max_y - min_y, max_x - min_x)

    def electrode(self, identifier):
        return self.electrodes[identifier]

    def __repr__(self):
        return "Recording(name=%r, n_pos=%d, electrodes=%s)" % (
            self.name, self.n_pos, sorted(self.electrodes),
        )
