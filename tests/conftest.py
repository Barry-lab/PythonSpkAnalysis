"""Shared test fixtures.

Golden-master tests need the example DACQ trial + reference ``.mat`` that ship
with the MATLAB ``Universal_Matlab`` repo. If they are not found (e.g. in CI
without the data) those tests skip; the synthetic invariant tests always run.

Point ``BARRELAB_EXAMPLE_DATA`` at the ``Example_data`` directory to enable them.
"""

import os

import pytest

# Ensure the src/ layout is importable without an editable install.
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

_CANDIDATES = [
    os.environ.get("BARRELAB_EXAMPLE_DATA", ""),
    "/Users/caswell/Documents/Git/Universal_Matlab/Read_Bin_Data/Example_data",
]


def _example_dir():
    for cand in _CANDIDATES:
        if cand and os.path.isdir(cand):
            return cand
    return None


@pytest.fixture(scope="session")
def example_trial():
    d = _example_dir()
    trial = None if d is None else os.path.join(d, "r1838_10-12-14", "83812141")
    if trial is None or not os.path.exists(trial + ".set"):
        pytest.skip("example DACQ trial not available")
    return trial


@pytest.fixture(scope="session")
def golden_mat():
    d = _example_dir()
    path = None if d is None else os.path.join(d, "gridCell1838_t2c3_14-12-10.mat")
    if path is None or not os.path.exists(path):
        pytest.skip("golden .mat not available")
    sio = pytest.importorskip("scipy.io")
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)
