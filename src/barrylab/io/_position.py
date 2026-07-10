"""Position post-processing for DACQ tracking data.

Ports ``postprocess_pos_data.m`` and its helpers ``led_swap_filter.m`` and
``led_speed_filter.m``. Turns the raw per-LED pixel positions read from a
``.pos`` file into cleaned ``xy`` (pixels), head ``direction`` (degrees) and
``speed`` (cm/s), reproducing the historic pipeline: swap-correct, speed-reject,
interpolate, boxcar-smooth, then derive direction and speed.
"""

import numpy as np
from scipy.ndimage import convolve1d

from ._axona_files import key_value


def _boxcar_replicate(arr, kernel):
    """Emulate MATLAB ``imfilter(arr, kernel, 'same', 'conv', 'replicate')``.

    ``arr`` is smoothed column-wise (axis 0) with the 1-D ``kernel`` using edge
    replication. The kernel is symmetric (a boxcar) so convolution and
    correlation coincide, and scipy's default centring (``origin=0``) matches
    MATLAB's for both odd and even kernels — verified bit-exact against the
    lab's golden position output for the 20-tap (0.4 s) kernel.
    """
    return convolve1d(arr, kernel, axis=0, mode="nearest")


def led_swap_filter(led_pos, led_pix, thresh=5.0):
    """Detect samples where the two tracked LEDs have swapped identity.

    Returns the (0-based) row indices whose two LEDs should be exchanged.
    Direct port of ``led_swap_filter.m`` (euclidean-distance variant).
    """
    mean_npix = np.nanmean(led_pix, axis=0)
    std_npix = np.nanstd(led_pix, axis=0, ddof=1)  # MATLAB nanstd uses N-1

    cur = slice(1, led_pix.shape[0])
    prev = slice(0, led_pix.shape[0] - 1)

    def dist(a, b):
        diff = led_pos[cur, a, :] - led_pos[prev, b, :]
        return np.sqrt(np.nansum(diff ** 2, axis=1))

    dist12, dist11 = dist(0, 1), dist(0, 0)
    dist21, dist22 = dist(1, 0), dist(1, 1)

    switched = (dist12 < dist11 - thresh) & (
        np.isnan(led_pos[cur, 1, 0]) | (dist21 < dist22 - thresh)
    )
    z11 = (mean_npix[0] - led_pix[cur, 0]) / std_npix[0]
    z12 = (led_pix[cur, 0] - mean_npix[1]) / std_npix[1]
    with np.errstate(invalid="ignore"):  # NaN pix -> NaN compare -> False, as in MATLAB
        shrunk = z11 > z12
    return np.where(switched & shrunk)[0] + 1


def led_speed_filter(led_pos, max_pix_per_sample, led):
    """NaN-out points reached impossibly fast (tracker grabbing a stray point).

    Mutates ``led_pos`` in place for the given LED. Sequential port of
    ``led_speed_filter.m``: the previous *accepted* point is the reference.
    """
    ok = np.where(~np.isnan(led_pos[:, led, 0]))[0]
    if len(ok) < 2:
        return
    max_sq = max_pix_per_sample ** 2
    prev = ok[0]
    for idx in range(1, len(ok)):
        pos = ok[idx]
        step_sq = (
            (led_pos[pos, led, 0] - led_pos[prev, led, 0]) ** 2
            + (led_pos[pos, led, 1] - led_pos[prev, led, 1]) ** 2
        ) / (pos - prev) ** 2
        if step_sq > max_sq:
            led_pos[pos, led, :] = np.nan
        else:
            prev = pos


def _interp_nans(led_pos, led):
    """Linearly interpolate NaN gaps for one LED, holding the endpoints."""
    n = led_pos.shape[0]
    idx = np.arange(n)
    missing = np.isnan(led_pos[:, led, 0]) | np.isnan(led_pos[:, led, 1])
    ok = ~missing
    if not ok.any():
        return
    for k in range(2):
        # np.interp holds the endpoint value outside the sample range, which
        # matches the explicit endpoint handling in postprocess_pos_data.m.
        led_pos[missing, led, k] = np.interp(idx[missing], idx[ok], led_pos[ok, led, k])


def postprocess_pos_data(led_pos, led_pix, pos_header, set_header,
                         max_speed=2.0, box_car=0.4, head_pos_weighting="tracking"):
    """Clean raw LED positions into xy / direction / speed.

    Parameters
    ----------
    led_pos : ndarray (n_pos, n_leds, 2)
        Raw LED positions in pixels (``x`` then ``y``), with untracked samples
        already flagged as NaN by the ``.pos`` reader.
    led_pix : ndarray (n_pos, n_leds)
        Pixel count per LED (used for swap detection).
    pos_header, set_header
        Parsed ``.pos`` and ``.set`` headers.
    max_speed : float
        Jump-rejection threshold in m/s.
    box_car : float
        Position smoothing boxcar width in seconds.
    head_pos_weighting : {"tracking", "midpoint"}
        How to combine the two LEDs into a head position (2-LED trials only).
        ``"tracking"`` (default) reproduces the current MATLAB
        ``postprocess_pos_data.m``: a mean weighted by how often each LED was
        tracked, which is robust to a poorly-tracked light. ``"midpoint"`` is
        the plain geometric midpoint of the two lights (the legacy convention;
        this is what the shipped example ``.mat`` was generated with).

    Returns
    -------
    dict with keys ``xy`` (n,2), ``direction`` (n,), ``speed`` (n,),
    ``times`` (n,), ``n_leds``.
    """
    if head_pos_weighting not in ("tracking", "midpoint"):
        raise ValueError("head_pos_weighting must be 'tracking' or 'midpoint'")
    led_pos = np.array(led_pos, dtype=float, copy=True)
    led_pix = None if led_pix is None else np.array(led_pix, dtype=float, copy=True)
    n_pos, n_leds = led_pos.shape[0], led_pos.shape[1]

    sample_rate = key_value(pos_header, "sample_rate", "num")
    ppm = key_value(pos_header, "pixels_per_metre", "num")

    # Blend weights: fraction of samples each LED was actually tracked for
    # (computed before any swapping, matching the MATLAB order).
    weights = np.sum(~np.isnan(led_pos[:, :, 0]), axis=0) / float(n_pos)

    if n_leds == 2 and led_pix is not None and led_pix.size:
        swap = led_swap_filter(led_pos, led_pix)
        led_pos[swap, 0, :], led_pos[swap, 1, :] = (
            led_pos[swap, 1, :].copy(),
            led_pos[swap, 0, :].copy(),
        )
        led_pix[swap, 0], led_pix[swap, 1] = (
            led_pix[swap, 1].copy(),
            led_pix[swap, 0].copy(),
        )

    max_pix_per_sample = max_speed * ppm / sample_rate
    for led in range(n_leds):
        led_speed_filter(led_pos, max_pix_per_sample, led)
    for led in range(n_leds):
        _interp_nans(led_pos, led)

    if box_car > 0:
        width = int(np.ceil(box_car * sample_rate))
        kernel = np.ones(width) / width
    else:
        kernel = np.array([1.0])

    times = np.arange(1, n_pos + 1) / sample_rate

    if n_leds == 1:
        xy = _boxcar_replicate(led_pos[:, 0, :], kernel)
        direction = _heading_from_displacement(xy, sample_rate)
    else:
        sm_front = _boxcar_replicate(led_pos[:, 0, :], kernel)
        sm_back = _boxcar_replicate(led_pos[:, 1, :], kernel)
        # 0 deg = straight ahead; the large (front) light's set-file bearing
        # corrects for its offset from the animal's head axis.
        correction = key_value(set_header, "lightBearing_1", "num")
        dy = -sm_front[:, 1] + sm_back[:, 1]
        dx = sm_front[:, 0] - sm_back[:, 0]
        direction = np.mod(np.degrees(np.arctan2(dy, dx)) - correction, 360.0)
        if head_pos_weighting == "tracking":
            xy = (weights[0] * sm_front + weights[1] * sm_back) / np.sum(weights)
        else:  # "midpoint"
            xy = (sm_front + sm_back) / 2.0

    speed = np.empty(n_pos)
    step = np.sqrt(np.diff(xy[:, 0]) ** 2 + np.diff(xy[:, 1]) ** 2)
    speed[: n_pos - 1] = step
    speed[n_pos - 1] = speed[n_pos - 2] if n_pos > 1 else 0.0
    speed = speed * (100.0 * sample_rate / ppm)  # pixels/sample -> cm/s

    return {
        "xy": xy,
        "direction": np.asarray(direction).reshape(-1),
        "speed": speed,
        "times": times,
        "n_leds": n_leds,
    }


def _heading_from_displacement(xy, sample_rate):
    """Direction (deg, anticlockwise from +x) from frame-to-frame displacement."""
    n = xy.shape[0]
    direction = np.empty(n)
    dy = -xy[1:, 1] + xy[:-1, 1]
    dx = xy[1:, 0] - xy[:-1, 0]
    direction[: n - 1] = np.mod(np.degrees(np.arctan2(dy, dx)), 360.0)
    direction[n - 1] = direction[n - 2] if n > 1 else 0.0
    return direction
