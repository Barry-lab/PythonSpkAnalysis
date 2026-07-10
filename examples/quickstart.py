"""Load a DACQ trial, make a rate map, and print spatial information.

    python3 examples/quickstart.py /path/to/r1838_10-12-14/83812141 2 3

Saves a rate-map + trajectory figure to examples/output/ if matplotlib is
installed.
"""

import os
import sys

# Make the package importable when run from a clone without installing.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np

import barrelab as bl


def main(trial, electrode, cluster):
    rec = bl.load_dacq(trial)
    print(rec)
    print("  duration: %.1f s   pixels/m: %g" % (rec.duration, rec.pixels_per_metre))

    elec = rec.electrode(electrode)
    print("  tetrode %d: %d spikes, clusters %s"
          % (electrode, elec.spike_times.size, list(elec.cluster_ids())))

    rmap = bl.rate_map(rec, electrode, cluster)
    occ = bl.occupancy_map(rec)
    bits_per_spike, bits_per_s = bl.spatial_information(rmap, occ)
    print("  cluster %d: peak %.2f Hz, %.3f bits/spike, %.3f bits/s"
          % (cluster, np.nanmax(rmap), bits_per_spike, bits_per_s))

    _try_plot(rec, elec, cluster, rmap)


def _try_plot(rec, elec, cluster, rmap):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (install matplotlib to save a figure)")
        return

    out_dir = os.path.join(os.path.dirname(__file__), "output")
    os.makedirs(out_dir, exist_ok=True)
    spk_idx = elec.spike_pos_index(cluster)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.6))
    ax[0].plot(rec.xy[:, 0], rec.xy[:, 1], color="0.7", lw=0.4)
    ax[0].scatter(rec.xy[spk_idx, 0], rec.xy[spk_idx, 1], s=6, color="crimson")
    ax[0].set_title("trajectory + spikes (t%d c%d)" % (elec.identifier, cluster))
    ax[0].set_aspect("equal"); ax[0].invert_yaxis()

    im = ax[1].imshow(rmap, origin="lower", interpolation="nearest")
    ax[1].set_title("rate map  (peak %.1f Hz)" % np.nanmax(rmap))
    fig.colorbar(im, ax=ax[1], label="Hz", fraction=0.046)
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout()
    out = os.path.join(out_dir, "ratemap_t%dc%d.png" % (elec.identifier, cluster))
    fig.savefig(out, dpi=130)
    print("  saved %s" % out)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
