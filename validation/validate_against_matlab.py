"""Validate the Python port against the MATLAB golden output.

Loads the reference ``gridCell1838_t2c3_14-12-10.mat`` shipped with the MATLAB
toolbox (the fully-loaded ``data`` struct plus the golden rate map ``smthRm``
for tetrode 2, cluster 3), runs the barrylab port over the same raw trial, and
reports the discrepancy at every stage.

What is proven bit-exact here (machine precision):
  * the .pos / tetrode / cut binary readers,
  * the full position pipeline (swap, speed filter, interpolation, smoothing)
    in ``midpoint`` mode,
  * the spike-time -> position-sample mapping and cluster selection,
  * the spatial binning (``histnd``), also cross-checked brute-force.

What is NOT bit-exact against THIS fixture (and why): the final smoothed rate
map. The shipped README states ``smthRm`` was made with
"make_smooth_ratemap.m (or its precursor)"; empirically the fixture matches a
Gaussian (sigma~2 bins, r=0.98), not the current 5-bin boxcar (r=0.94). The
fixture's smoothing step is therefore stale relative to the current code.
See ``FINDINGS.md``. To bit-validate the smoothing and the (default)
tracking-weighted head position, regenerate fixtures from the *current* MATLAB
with ``validation/generate_fixtures.m``.

Run:  python3 validation/validate_against_matlab.py
"""

import os
import sys

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from barrylab.analysis.maps import occupancy_map, rate_map, spike_map  # noqa: E402
from barrylab.analysis.spatial import spatial_information  # noqa: E402
from barrylab.io.axona import _read_pos, load_dacq  # noqa: E402

REPO = "/Users/caswell/Documents/Git/Universal_Matlab"
EX = os.path.join(REPO, "Read_Bin_Data", "Example_data")
TRIAL = os.path.join(EX, "r1838_10-12-14", "83812141")
GOLDEN = os.path.join(EX, "gridCell1838_t2c3_14-12-10.mat")

TET, CLUSTER = 2, 3
EXACT = 1e-9
results = []


def check(name, ok, detail=""):
    results.append(ok)
    print("  [%s] %-32s %s" % ("PASS" if ok else "FAIL", name, detail))


def note(name, detail):
    print("  [note] %-32s %s" % (name, detail))


def nanmaxabs(a, b):
    return float(np.nanmax(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def main():
    m = sio.loadmat(GOLDEN, squeeze_me=True, struct_as_record=False)
    g = m["data"]
    g_pos = g.pos
    g_tet = next(t for t in g.tetrode if int(t.id) == TET)

    print("\n=== Stage 1: raw .pos reader (bit-exact) ===")
    led_pos, _, _ = _read_pos(TRIAL + ".pos")
    check("led_pos shape", led_pos.shape == g_pos.led_pos.shape)
    check("led_pos NaN pattern", np.array_equal(np.isnan(led_pos), np.isnan(g_pos.led_pos)))
    finite = ~np.isnan(led_pos)
    check("led_pos values exact", np.array_equal(led_pos[finite], g_pos.led_pos[finite]))

    print("\n=== Stage 2: position pipeline, midpoint mode (bit-exact) ===")
    rec = load_dacq(TRIAL, head_pos_weighting="midpoint")
    d_xy, d_sp = nanmaxabs(rec.xy, g_pos.xy), nanmaxabs(rec.speed, g_pos.speed)
    dd = float(np.nanmax(np.abs(np.mod(rec.direction - g_pos.dir + 180, 360) - 180)))
    check("xy", d_xy < EXACT, f"max|d|={d_xy:.2e} px")
    check("direction", dd < EXACT, f"max|d|={dd:.2e} deg")
    check("speed", d_sp < EXACT, f"max|d|={d_sp:.2e} cm/s")

    print("\n=== Stage 3: spike timing -> pos_sample (bit-exact) ===")
    elec = rec.electrode(TET)
    check("pos_sample", np.array_equal(elec.pos_index + 1, np.asarray(g_tet.pos_sample)))

    print("\n=== Stage 4: cut file + spike selection (bit-exact) ===")
    check("exactCut", np.array_equal(elec.clusters, m["exactCut"]))
    check(
        "spikePosSamp (t2c3)", np.array_equal(elec.spike_pos_index(CLUSTER) + 1, m["spikePosSamp"])
    )

    print("\n=== Stage 5: spatial binning (bit-exact, brute-force cross-check) ===")
    occ = occupancy_map(rec)
    spk = spike_map(rec, elec.spike_pos_index(CLUSTER))
    check("occupancy total == duration", abs(occ.sum() - rec.duration) < 1e-6, f"{occ.sum():.4f} s")
    check(
        "spike total == n spikes",
        spk.sum() == elec.spike_pos_index(CLUSTER).size,
        "%d spikes" % int(spk.sum()),
    )
    check("binning == brute force", _binning_matches_brute(rec, occ, spk, elec))

    print("\n=== Stage 5b: smoothed rate map vs stale fixture (sanity, NOT bit-exact) ===")
    rm = rate_map(rec, TET, CLUSTER)
    gold_rm = m["smthRm"]
    check(
        "rate map shape + NaN pattern",
        rm.shape == gold_rm.shape and np.array_equal(np.isnan(rm), np.isnan(gold_rm)),
    )
    both = ~np.isnan(rm) & ~np.isnan(gold_rm)
    corr = float(np.corrcoef(rm[both], gold_rm[both])[0, 1])
    note("current boxcar vs fixture", f"r={corr:.4f} (fixture used a precursor smoother)")
    check("rate map correlated (sanity)", corr > 0.9, f"r={corr:.4f}")

    print("\n=== Stage 6: spatial information (sanity) ===")
    bps, bpsec = spatial_information(rm, occ)
    check(
        "spatial info finite & positive",
        np.isfinite(bps) and bps > 0,
        f"{bps:.4f} bits/spike, {bpsec:.4f} bits/s",
    )

    print("\n" + "=" * 60)
    n = sum(results)
    print("RESULT: %d/%d checks passed" % (n, len(results)))
    return 0 if all(results) else 1


def _binning_matches_brute(rec, occ, spk, elec):
    bs = (rec.pixels_per_metre / 100.0) * 2.0
    yext, xext = rec.window_extent
    nby, nbx = int(np.ceil(yext / bs)), int(np.ceil(xext / bs))

    def brute(points):
        yb = np.clip(
            np.ceil(np.where(points[:, 1] == 0, 1e-12, points[:, 1]) / bs).astype(int), 1, nby
        )
        xb = np.clip(
            np.ceil(np.where(points[:, 0] == 0, 1e-12, points[:, 0]) / bs).astype(int), 1, nbx
        )
        out = np.zeros((nby, nbx))
        np.add.at(out, (yb - 1, xb - 1), 1.0)
        return out

    occ_bf = brute(rec.xy) / rec.pos_sample_rate
    spk_bf = brute(rec.xy[elec.spike_pos_index(3)])
    return np.max(np.abs(occ_bf - occ)) == 0 and np.max(np.abs(spk_bf - spk)) == 0


if __name__ == "__main__":
    sys.exit(main())
