# barrylab documentation

`barrylab` is a modern, validated Python successor to the Barry lab's `Universal_Matlab` toolbox: it reads historic Axona/DACQ tetrode data today and is built to grow into Neuropixels without rewriting analysis code. Start with getting started, then follow the links below for architecture, migration notes, and how to extend the package.

- [Getting started](getting_started.md) — for new users (often coming from MATLAB): install `barrylab`, load a DACQ trial, and make your first rate map.
- [Architecture](architecture.md) — for anyone extending the code: the three-layer design (raw files → loaders → `Recording` → analysis) and why one backend serves many frontends.
- [MATLAB migration](matlab_migration.md) — for `Universal_Matlab` users: a function-by-function map from the old MATLAB scripts to their `barrylab` equivalents.
- [Adding a loader](adding_a_loader.md) — for developers adding a new data source: how to fill a `Recording` so every existing analysis works unchanged.
- [Contributing](../CONTRIBUTING.md) — for contributors: dev setup, coding conventions, the Python 3.9+ constraint, and how to run tests, linting, and validation.
- [Findings](../FINDINGS.md) — for maintainers and the PI: notes from the port, two documented behavioural choices that differ from the shipped fixture, and MATLAB quirks spotted along the way.
- [Validation](../validation/README.md) — for anyone trusting the numbers: how ported functions are checked against golden-master MATLAB output and how to regenerate fixtures.
- [Project README](../README.md) — for a first look: the elevator pitch, quickstart, package layout, and roadmap from v0.1 to v1.0.
