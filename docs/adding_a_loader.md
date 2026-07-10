# Adding a new front-end loader

barrelab has three strict layers:

```
RAW FILES  ->  FRONT-END LOADERS (io/)  ->  COMMON OBJECT (core/Recording)  ->  ANALYSIS (analysis/)
```

A loader's only job is to read a raw acquisition format and return a fully
populated `Recording`. Every analysis consumes a `Recording` and nothing else,
so a rate map does not care whether the spikes came from a 2010 Axona tetrode or
a 2026 Neuropixel. To support a new acquisition system you write a new loader;
you touch no analysis code.

Today `load_dacq` (Axona/DACQ) is the **only implemented loader**. Everything
below describes the contract it satisfies, so the next loader can satisfy the
same one.

## The loader contract

A loader is a plain function (or class) that returns a
`barrelab.core.recording.Recording`. It fills the following fields. Reading
`src/barrelab/io/axona.py` alongside this list shows exactly how `load_dacq`
does it.

| Field | Type | Meaning | Optional? |
| --- | --- | --- | --- |
| `name` | str | Session name (source file stem) | required |
| `xy` | ndarray `(n_pos, 2)` | Position in pixels, tracked-window frame | optional* |
| `direction` | ndarray `(n_pos,)` | Head direction, degrees anticlockwise from +x | optional* |
| `speed` | ndarray `(n_pos,)` | Running speed in cm/s | optional* |
| `pos_times` | ndarray `(n_pos,)` | Time of each position sample, seconds | optional* |
| `pos_sample_rate` | float | Position sampling rate, Hz | optional* |
| `pixels_per_metre` | float | Pixels-to-cm scale | optional* |
| `window` | tuple `(min_x, max_x, min_y, max_y)` | Tracked camera window, pixels | optional* |
| `electrodes` | dict `{id: Electrode}` | Sorted spikes keyed by group id | required in practice |
| `metadata` | dict | Provenance (source path, format, raw headers) | recommended |

`*` The position block is optional as a block: `load_dacq(load_position=False)`
returns a `Recording` with all position fields left as `None`, and spatial
analyses simply cannot run on it. If you fill position at all, fill it
consistently — every position array must share the same `n_pos`, and
`window`/`pixels_per_metre` must match the frame `xy` lives in, because
`occupancy_map` and `rate_map` derive the grid from them.

Each `Electrode` (see `core/recording.py`) is:

```python
Electrode(identifier, spike_times, pos_index=None, clusters=None, metadata=None)
```

- `identifier` — tetrode number or probe-group id (the dict key).
- `spike_times` — ndarray of spike times in **seconds**.
- `pos_index` — for each spike, the 0-based index into the position arrays of
  the co-occurring sample. `None` if there is no position. `load_dacq` computes
  this from `spike_times * pos_sample_rate` and drops spikes past the end of the
  position clock; a new loader must produce the same 0-based mapping.
- `clusters` — cluster id per spike (0 = unsorted/noise), or `None` if uncut.
- `metadata` — the raw header, kept for provenance.

`spike_times`, `pos_index` and `clusters` must line up element-for-element (before
`pos_index` truncation), so downstream `spike_mask` / `spike_pos_index` select
the right samples.

## Planned: `load_neuropixels` (design sketch, not yet implemented)

The v0.3 roadmap adds a Neuropixels front end. The intended shape — **not yet
written** — reuses the same contract:

```python
def load_neuropixels(path, position_source=None):
    """PLANNED — not implemented. Assemble a Recording from a Neuropixels trial."""
    import spikeinterface.full as si

    # 1. Spikes + probe geometry via SpikeInterface (SpikeGLX or Open Ephys
    #    recording + a Kilosort4 sorting output).
    recording_ext = si.read_spikeglx(path)          # or si.read_openephys(path)
    sorting = si.read_kilosort(path / "kilosort4")

    # 2. Position / behaviour + sync via a lab-owned step that aligns the
    #    behavioural clock to the ephys clock and returns xy, direction,
    #    speed, pos_times, pos_sample_rate, window, pixels_per_metre.
    pos = load_position_and_sync(position_source)

    # 3. Group units into Electrodes (e.g. one group, or by shank/depth) and
    #    map each spike to a position sample:
    #       pos_index = searchsorted(pos["times"], spike_times) - 1
    #    Then build the Recording exactly as load_dacq does.
    rec = Recording(name=path.stem, metadata={"format": "neuropixels"})
    # ... fill rec.xy, rec.direction, ..., rec.electrodes ...
    return rec
```

SpikeInterface and probeinterface arrive with the `[neuropixels]` extra and are
never a runtime dependency of the core. The position/sync step is lab-specific
and does not exist yet.

## Checklist for a new loader

1. **Fill the contract.** Return a `Recording` with the fields above; match
   units (seconds, cm/s, pixels) and the 0-based `pos_index` convention.
2. **Add a test.** Put a small fixture and a test under `tests/` that loads it
   and asserts shapes, dtypes and a couple of known values.
3. **Validate against ground truth.** Follow the golden-master philosophy: check
   the loaded numbers against the acquisition system's own tools (e.g.
   SpikeInterface / ephysiopy for Neuropixels) rather than only internal
   consistency.
4. **Register it.** Add the function to `src/barrelab/io/__init__.py` (its
   `__all__`) and re-export it from `src/barrelab/__init__.py` so users can
   `from barrelab import load_neuropixels`.
