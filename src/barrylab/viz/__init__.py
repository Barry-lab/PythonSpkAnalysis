"""Optional plotting helpers (requires the ``plot`` extra: ``pip install "barrylab[plot]"``).

Not imported by the top-level ``barrylab`` package (matching the
``analysis/field_shuffle.py``/``scikit-image`` precedent), so ``import
barrylab`` never requires ``matplotlib``. Import directly, e.g.::

    from barrylab.viz.plots import image_rm
    from barrylab.viz.colormaps import TINT_CMAP
"""
