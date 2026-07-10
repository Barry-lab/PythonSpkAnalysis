# Architecture

`barrelab` (imported as `import barrelab`; the repo is `Barry-lab/PythonSpkAnalysis`)
is a modern, validated successor to the lab's MATLAB `Universal_Matlab` toolbox. It
reads historic Axona/DACQ tetrode data today and is built to grow into Neuropixels
without rewriting a single line of analysis code. This document explains the design
so you know *where* new code goes and *why*.

## The three layers

Data flows one way through three strict layers:

```
RAW FILES  ->  FRONT-END LOADERS  ->  COMMON OBJECT  ->  ANALYSIS
 (.pos,        (io/: load_dacq,       (core/Recording   (analysis/: rate_map,
  .set,         _axona_files,          + Electrode)       spatial_information,
  .N, .cut)     _position)                                adaptive_binning)
```

- **Front-end loaders** (`src/barrelab/io/`) are the only code that knows about
  byte layouts and file extensions. `load_dacq` (in `io/axona.py`) reads a DACQ
  trial via the internal `_read_pos`, `_read_tetrode` and `_read_cut` helpers plus
  the `_position` pipeline, and returns a filled `Recording`. A loader's *only*
  job is to populate a `Recording`.
- **The common object** (`core/recording.py`) is the contract between the two
  halves of the codebase. It is a plain, transparent container of NumPy arrays and
  light metadata — no file-format knowledge at all.
- **Analysis** (`analysis/`) consumes a `Recording` and nothing else. Read
  `analysis/maps.py` and `analysis/spatial.py`: every public function takes a
  `Recording` (or plain arrays) and none of them ever open a file path.

**The invariant that makes multi-format work:** analyses depend only on the
`Recording`. Add a Neuropixels loader that fills the same object and every existing
rate map, spatial-information and adaptive-binning function works unchanged.

## The Recording / Electrode contract

A loader must populate these `Recording` fields (see `core/recording.py`):

| Field | Meaning |
|-------|---------|
| `name` | session name (source file stem) |
| `xy` | `(n_pos, 2)` position in pixels, tracked-window frame |
| `direction` | `(n_pos,)` head direction, degrees anticlockwise from +x |
| `speed` | `(n_pos,)` running speed, cm/s |
| `pos_times` | `(n_pos,)` sample times, seconds |
| `pos_sample_rate` | position rate, Hz |
| `pixels_per_metre` | pixel-to-cm scale |
| `window` | `(min_x, max_x, min_y, max_y)` camera window, pixels |
| `electrodes` | `dict[int, Electrode]` keyed by tetrode/group id |
| `metadata` | free-form provenance (source path, format, raw headers) |

Derived helpers — `n_pos`, `duration`, `window_extent` and `electrode(id)` — are
computed, not stored. Each `Electrode` carries `identifier`, `spike_times` (s),
`pos_index` (position sample co-occurring with each spike), `clusters` (cut ids),
and `metadata`, with methods `cluster_ids()`, `spike_mask(cluster)` and
`spike_pos_index(cluster)`.

What each analysis reads:

- `occupancy_map` uses `xy`, `window_extent`, `pixels_per_metre`, `pos_sample_rate`.
- `spike_map` / `rate_map` additionally use `electrode(id).spike_pos_index(cluster)`.
- `spatial_information` and `adaptive_binning` take only binned arrays (rates and
  dwell times), so they are format-agnostic by construction.

If a loader fills these fields correctly, every analysis is guaranteed to work.

## Key design decisions

**A lab-owned, lightweight `Recording` — not NWB or pynapple.** The common object
is deliberately a plain class of arrays a student can print and inspect field by
field. This keeps the core transparent, dependency-light and debuggable by rotating
students and postdocs, many arriving from MATLAB. NWB is powerful but heavyweight
and opaque for day-to-day debugging; we reserve it for export and interchange
(roadmap v0.4+), not as the in-memory representation.

**Adopt SpikeInterface/Kilosort for Neuropixels loading, but own the spatial
metrics.** When the Neuropixels front end lands (roadmap v0.3), spike sorting and
raw I/O come from mature community tools (SpikeInterface, Kilosort4, SpikeGLX and
Open Ephys) feeding into the *same* `Recording`. But the spatial-coding metrics —
rate maps, Skaggs information, adaptive binning — are the lab's validated
intellectual property, ported from `Universal_Matlab` and checked numerically
against MATLAB output (golden-master validation). We keep and own those.

**Analysis functions never touch a file path.** This is enforced by convention and
visible in the code: `maps.py` and `spatial.py` import no I/O. It is what lets a
rate map ignore whether spikes came from a 2010 tetrode or a 2026 Neuropixel.

## Where does my code go?

- New acquisition format -> a new loader in `io/` that fills a `Recording`.
- New spatial or temporal metric -> a function in `analysis/` taking a `Recording`.
- New stored quantity every analysis needs -> a field on `Recording`/`Electrode`.

If you find yourself parsing a file inside `analysis/`, or importing `io` from it,
stop: that logic belongs in a loader.
