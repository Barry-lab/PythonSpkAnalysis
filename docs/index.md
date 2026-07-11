# barrylab documentation

`barrylab` is a Python port of the Barry lab's `Universal_Matlab` toolbox for
spatial-coding analysis. It reads Axona/DACQ tetrode data (the mature, validated
path) and, experimentally, Neuropixels data — the same analysis code runs on
both. Start with getting started, then follow the links below for architecture,
migration notes, and how to extend the package.

- [Getting started](getting_started.md) — for new users (often coming from MATLAB): install `barrylab`, load a DACQ trial, and make your first rate map.
- [Architecture](architecture.md) — for anyone extending the code: the three-layer design (raw files → loaders → `Recording` → analysis) and why one backend serves many frontends.
- [MATLAB migration](matlab_migration.md) — for `Universal_Matlab` users: a function-by-function map from the old MATLAB scripts to their `barrylab` equivalents.
- [Adding a loader](adding_a_loader.md) — for developers adding a new data source: how to fill a `Recording` so every existing analysis works unchanged (the DACQ and Neuropixels loaders are the two worked examples).
- [Contributing](../CONTRIBUTING.md) — for contributors: dev setup, coding conventions, the Python 3.9+ constraint (kept 3.6-runnable), and how to run tests, linting, and validation.
- [Findings](../FINDINGS.md) — for maintainers and the PI: the per-function validation record and every MATLAB bug/quirk found (and how it was handled) during the port.
- [Validation](../validation/README.md) — for anyone trusting the numbers: how ported functions are checked against golden-master MATLAB output and how to regenerate fixtures.
- [Project README](../README.md) — for a first look: what the package does, a quickstart, a tour of every analysis area, and the Neuropixels notes.
