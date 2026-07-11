# Adding a new front-end loader

barrylab has three strict layers:

```
RAW FILES  ->  FRONT-END LOADERS (io/)  ->  COMMON OBJECT (core/Recording)  ->  ANALYSIS (analysis/)
```

A loader's only job is to read a raw acquisition format and return a fully
populated `Recording`. Every analysis consumes a `Recording` and nothing else,
so a rate map does not care whether the spikes came from a 2010 Axona tetrode or
a 2026 Neuropixel. To support a new acquisition system you write a new loader;
you touch no analysis code.

There are two implemented loaders: `load_dacq` (Axona/DACQ, the mature,
golden-master-validated one) and `load_neuropixels` (experimental). They are
worth reading side by side — `load_dacq` reads raw binary files, while
`load_neuropixels` composes *already-sorted* spikes (from SpikeInterface) with a
separately-supplied position track — because they show the same contract
satisfied two very different ways. Everything below describes that contract, so
the next loader can satisfy the same one.

## The loader contract

A loader is a plain function (or class) that returns a
`barrylab.core.recording.Recording`. It fills the following fields. Reading
`src/barrylab/io/axona.py` alongside this list shows exactly how `load_dacq`
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

## The second loader: `load_neuropixels` (a worked example)

`load_neuropixels` (`src/barrylab/io/neuropixels.py`, experimental) satisfies the
same contract a very different way, so it's the clearest example of what a
composite loader looks like. It does **not** read raw binary or spike-sort —
that's Kilosort's job upstream. It reads *already-sorted* spikes plus a
*separately-supplied* position track and fills a `Recording`:

```python
def load_neuropixels(sorting, positions, pos_sample_rate, pixels_per_metre, ...):
    # 1. Spikes: `sorting` is a Kilosort/phy output folder (read via
    #    SpikeInterface's read_phy/read_kilosort, the [neuropixels] extra), a
    #    SpikeInterface Sorting object, or any object exposing get_unit_ids() /
    #    get_unit_spike_train(unit_id) / sampling_frequency.
    # 2. Position: a Neuropixels probe records neural data only, so `positions`
    #    (n_pos, 2) comes from your behavioural source (VR / camera / DLC) on
    #    its own clock, which the caller must already have aligned to the
    #    neural clock's t=0. direction/speed are derived from it if not given.
    # 3. All units become one Electrode (the probe): clusters = sorted unit ids,
    #    pos_index = floor(spike_time * pos_sample_rate). rate_map(rec,
    #    electrode=probe_id, cluster=unit_id) then works exactly as for DACQ.
```

Contrast with `load_dacq`, which reads raw Axona binary files directly: same
`Recording` contract, completely different plumbing. SpikeInterface arrives
only with the `[neuropixels]` extra and is never a runtime dependency of the
core — if you already hold a `Sorting` object, no extra install is needed.
`load_neuropixels` is experimental and, unlike the bit-exact DACQ loader, has
no golden-master validation (see the README) — a good target to harden.

## Checklist for a new loader

1. **Fill the contract.** Return a `Recording` with the fields above; match
   units (seconds, cm/s, pixels) and the 0-based `pos_index` convention.
2. **Add a test.** Put a small fixture and a test under `tests/` that loads it
   and asserts shapes, dtypes and a couple of known values.
3. **Validate against ground truth.** Follow the golden-master philosophy: check
   the loaded numbers against the acquisition system's own tools (e.g.
   SpikeInterface / ephysiopy for Neuropixels) rather than only internal
   consistency.
4. **Register it.** Add the function to `src/barrylab/io/__init__.py` (its
   `__all__`) and re-export it from `src/barrylab/__init__.py` so users can
   `from barrylab import load_neuropixels`.
