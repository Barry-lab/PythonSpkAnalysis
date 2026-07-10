# Contributing to barrylab

Thanks for helping build the lab's Python analysis toolbox. This guide covers how to
set up a development environment, run the tests and the MATLAB validation harness, and
the conventions we hold PRs to. The repo is `Barry-lab/PythonSpkAnalysis`; the import
name is `barrylab`.

## Development environment

We recommend [`uv`](https://github.com/astral-sh/uv), but plain `pip` works too.

```bash
# with uv (recommended)
uv venv
uv pip install -e ".[dev,plot]"

# or with pip
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,plot]"
```

The `[dev]` extra pulls in pytest, pytest-cov, hypothesis, ruff and mypy; `[plot]` adds
matplotlib. Two further extras exist for later phases: `[neuropixels]` (spikeinterface,
probeinterface) and `[nwb]` (pynwb, neuroconv). A DACQ-only user needs none of those, so
keep them optional.

## Running the tests

```bash
pytest          # 140+ tests (count grows as the port progresses)
```

Tests live under `tests/` and run against fixtures checked into the repo, so no MATLAB
is required.

## The MATLAB validation harness

We validate the port numerically against the lab's MATLAB output — this is the
**golden-master** philosophy. A ported function is not "done" when it looks right or
passes an internal-consistency check; it is done when it reproduces the MATLAB ground
truth to numerical tolerance. Run the shipped harness:

```bash
python3 validation/validate_against_matlab.py     # 15/15 checks
```

This loads the shipped example fixture (`gridCell1838_t2c3_14-12-10.mat`) and proves the
whole loader chain plus binning are bit-exact. It does not bit-validate the final
smoothing or the tracking-weighted head position, because that fixture predates those
code paths (see `FINDINGS.md`).

### Adding a new golden fixture

To lock down a stage the shipped fixture cannot cover, generate a fresh fixture from the
current MATLAB and add a regression test:

1. In MATLAB, with `Universal_Matlab` on the path, run `generate_fixtures`:

   ```matlab
   generate_fixtures('/path/to/r1838_10-12-14/83812141', 2, 3, 'tests/fixtures/t2c3.mat')
   ```

2. In a Python test, load the `.mat` and assert `np.allclose` per field (`xy`,
   `direction`, `speed`, `occupancy`, `spikeMap`, `smthRm`, `bits_per_spike`, ...).

Once the fixture is in place, every future change is guarded against silent numerical
drift from the MATLAB reference.

## Code conventions

- **Lint and types.** `ruff` and `mypy` are configured in `pyproject.toml` (line length
  100, target 3.9). Run `ruff check .` and `mypy src` before pushing; both must be green.
- **Docstrings.** Use numpy-style docstrings, as in `analysis/maps.py`. Document every
  parameter, its units (cm, bins, seconds, pixels) and its default.
- **The layering rule (hard).** Analysis functions never do file I/O. They consume a
  `Recording` (or plain arrays) and nothing else. A loader in `io/` fills a `Recording`;
  everything in `analysis/` takes one. This is what lets a rate map stay indifferent to
  whether spikes came from a 2010 tetrode or a 2026 Neuropixel — do not reach for a file
  path inside an analysis.
- **Tests travel with code.** Add a test with every new function. If the function ports a
  MATLAB routine, prefer a golden-master test over a hand-written expectation.

## Python 3.9 target and the local-3.6 caveat

The package targets **Python >=3.9**. However, the only interpreter available on the
original development machine was 3.6, so the runtime code deliberately avoids
3.7+-only syntax so it could be validated locally. This is why `Recording` is a plain
class rather than a `@dataclass`. When editing runtime code, keep to that constraint;
tooling and tests may of course assume 3.9.

## Branch and PR flow

- Branch off `main`; do not commit directly to it.
- Keep PRs small and reviewable — one coherent change, not a grab bag.
- Before opening a PR, make sure `pytest` passes and `ruff`/`mypy` are clean. If your
  change touches a validated code path, run `validate_against_matlab.py` too.
- Write a short description of what changed and why, and note any deviation from MATLAB
  behaviour (these are PI decisions and belong in `FINDINGS.md`).
