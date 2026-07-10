# Validation against MATLAB

The port is validated numerically against the lab's MATLAB output, not just
unit-tested for internal consistency. Two layers:

### 1. Against the shipped example fixture (works now, no MATLAB needed)

`gridCell1838_t2c3_14-12-10.mat` (in the MATLAB repo's `Read_Bin_Data/Example_data`)
holds the fully-loaded `data` struct and a golden rate map. Run:

```bash
python3 validation/validate_against_matlab.py
```

This proves the whole loader chain + binning bit-exact (see the script header and
`../FINDINGS.md`). It does **not** bit-validate the final smoothing or the default
head-position weighting, because this particular fixture predates those code paths.

### 2. Against fresh fixtures from the current MATLAB (recommended)

To lock down *every* stage — including the current tracking-weighted head position
and the current `make_smooth_ratemap` — regenerate fixtures from the current code:

```matlab
% in MATLAB, with Universal_Matlab on the path
generate_fixtures('/path/to/r1838_10-12-14/83812141', 2, 3, 'tests/fixtures/t2c3.mat')
```

Then a regression test loads that `.mat` and asserts `np.allclose` per field
(`xy`, `direction`, `speed`, `occupancy`, `spikeMap`, `smthRm`, `bits_per_spike`, …).
This is the golden-master harness recommended in the architecture plan: build it
once and every future change is guarded against silent numerical drift from the
MATLAB ground truth.
```
