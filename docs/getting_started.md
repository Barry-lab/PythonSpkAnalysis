# Getting started

This is a hands-on tour for a new lab member: load a real DACQ trial, look at
what came back, make a rate map, and quantify how spatial a cell is. If you have
used the MATLAB `Universal_Matlab` toolbox, this covers the same ground as
`read_DACQ` + `easy_make_ratemap` + `skaggsInfo`, in a few lines of Python.

## Install

Grab the code and install it in editable mode with the plotting extra:

```bash
git clone git@github.com:Barry-lab/PythonSpkAnalysis.git
cd PythonSpkAnalysis
uv pip install -e ".[dev,plot]"   # or: pip install -e ".[dev,plot]"
```

`[plot]` pulls in matplotlib for the figures below. For the full developer
setup (tests, linting, MATLAB validation) see [CONTRIBUTING](../CONTRIBUTING.md).

## Load a DACQ trial

A DACQ trial is a set of files that share one basename (`.pos`, `.set`, the
tetrode files `.1`, `.2`, ..., and the `.cut` cluster files). Point
`load_dacq` at *any* file of the trial, or at the shared basename with no
extension. It reads them all and returns a single `Recording`:

```python
import barrylab as bl

rec = bl.load_dacq("Universal_Matlab/Read_Bin_Data/Example_data/r1838_10-12-14/83812141")
```

The shipped example data lives in
`Universal_Matlab/Read_Bin_Data/Example_data/r1838_10-12-14/` (the trial basename
is `83812141`). In that trial, **tetrode 2, cluster 3 is a grid cell** — we use it
throughout.

## Inspect the Recording

A `Recording` is a plain, transparent container: position arrays plus a dict of
sorted `Electrode`s. Poke at it directly.

```python
print(rec)                      # Recording(name='83812141', n_pos=60050, electrodes=[1, 2, 3, 4, 5, 6, 7, 8])
print(rec.n_pos)                # number of position samples
print(rec.duration)             # session length in seconds
print(rec.pixels_per_metre)     # spatial scale
print(sorted(rec.electrodes))   # tetrode ids that were loaded

elec = rec.electrode(2)         # one tetrode
print(elec.spike_times.size)    # how many spikes it recorded
print(elec.cluster_ids())       # sorted unique cluster ids (0 = unsorted)
```

`rec.xy` is the `(n_pos, 2)` position track in pixels, `rec.speed` is running
speed in cm/s, and `rec.direction` is head direction in degrees — all sampled on
the same clock as `rec.pos_times`.

## Make a rate map

`rate_map` bins occupancy and spikes over space and divides them. With the
defaults it reproduces the lab-standard `easy_make_ratemap` (2 cm bins, 5-bin
boxcar smoothing):

```python
rmap = bl.rate_map(rec, electrode=2, cluster=3)   # boxcar default
```

The result is a `(n_y_bins, n_x_bins)` array of firing rate in Hz, with `NaN`
in bins the animal never visited. Swap the smoothing however you like:

```python
rmap_g   = bl.rate_map(rec, 2, 3, smoothing="gaussian", smooth_size=2)  # sigma in bins
rmap_raw = bl.rate_map(rec, 2, 3, smoothing=None)                       # unsmoothed
```

`smooth_size` means the boxcar width in bins for `"boxcar"`, or the Gaussian
standard deviation in bins for `"gaussian"`; leave it out for the sensible
default (5-bin box, or 2-bin sigma).

## Spatial information

The Skaggs spatial information (`skaggsInfo`) tells you how many bits each spike
carries about position. Feed it the rate map and the matching occupancy map:

```python
occ = bl.occupancy_map(rec)                         # dwell time (s) per bin
bits_per_spike, bits_per_sec = bl.spatial_information(rmap, occ)
print(f"{bits_per_spike:.3f} bits/spike, {bits_per_sec:.3f} bits/s")
```

A grid cell like t2 c3 should score well above a spatially uninformative cell.

## Plot trajectory, spikes and rate map

Finally, a quick figure that mirrors `examples/quickstart.py`: the animal's path
with spike locations on top, and the rate map beside it.

```python
import matplotlib.pyplot as plt

spk_idx = elec.spike_pos_index(3)   # position samples where cluster 3 fired

fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
ax[0].plot(rec.xy[:, 0], rec.xy[:, 1], color="0.7", lw=0.4)
ax[0].scatter(rec.xy[spk_idx, 0], rec.xy[spk_idx, 1], s=6, color="crimson")
ax[0].set_title("trajectory + spikes (t2 c3)")
ax[0].set_aspect("equal"); ax[0].invert_yaxis()

im = ax[1].imshow(rmap, origin="lower", interpolation="nearest")
ax[1].set_title("rate map")
fig.colorbar(im, ax=ax[1], label="Hz", fraction=0.046)
for a in ax:
    a.set_xticks([]); a.set_yticks([])
fig.tight_layout()
plt.show()
```

You should see the tell-tale hexagonal firing lattice of a grid cell.

## Where to go next

Everything above runs off one `Recording`, and every analysis takes only a
`Recording` — so the exact same code runs on Neuropixels data too, via
`bl.load_neuropixels(...)` (experimental — see the
[README](../README.md#neuropixels-experimental)). Run
`examples/quickstart.py <trial> 2 3` for the script version, browse the
["What's here" tour in the README](../README.md#whats-here) for every analysis
area, and see the analysis modules in `src/barrylab/analysis/`.
