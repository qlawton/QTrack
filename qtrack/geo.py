"""
Periodic-longitude helpers for QTrack.

QTrack ingests data on a longitude axis that is physically periodic, but the
original tracker treated it as a finite line running from the first array column
to the last. With ``prep_data`` producing a -180...180 axis, that array seam sits
exactly on the international dateline, so waves crossing 180 were lost.

Everything in this module is pure NumPy and works on scalars or arrays. The two
ideas that make the rest of the package wrap-safe are:

* ``lon_delta`` -- the *signed shortest* separation between two longitudes, which
  replaces every raw ``lon_a - lon_b`` / ``lon_a <= lon_b`` comparison; and
* ``pad_lon`` / ``pad_lon_coord`` -- circular padding of a field and its
  longitude coordinate, so that neighbourhood operations (gradients, radial
  averages, peak finding) see a complete neighbourhood at the seam.

Longitudes are returned in the [-180, 180) convention used by
``qtrack.core.prep_data``.
"""

import numpy as np

__all__ = [
    "wrap180",
    "wrap360",
    "lon_delta",
    "is_westward",
    "lon_in_range",
    "lon_in_any_range",
    "unwrap_lon",
    "on_unwrapped",
    "is_global_lon",
    "grid_resolution",
    "find_nearest_lon",
    "pad_lon",
    "unpad_lon",
    "pad_lon_coord",
    "wrap_index_delta",
]


def _scalar_in(*args):
    """True when every input is a scalar, so the result should be a scalar too."""
    return all(np.ndim(a) == 0 for a in args)


def _as_float(value, scalar):
    return float(value) if scalar else value


def _as_bool(value, scalar):
    return bool(value) if scalar else value


def wrap180(lon):
    """Wrap longitude(s) into [-180, 180).

    Matches the convention produced by ``prep_data`` (``(lon + 180) % 360 - 180``),
    so wrapped values can be compared directly against a prepped longitude axis.
    """
    out = (np.asarray(lon, dtype=float) + 180.0) % 360.0 - 180.0
    return _as_float(out, _scalar_in(lon))


def wrap360(lon):
    """Wrap longitude(s) into [0, 360)."""
    out = np.asarray(lon, dtype=float) % 360.0
    return _as_float(out, _scalar_in(lon))


def lon_delta(lon_b, lon_a):
    """Signed shortest separation ``lon_b - lon_a``, in [-180, 180).

    Negative means ``lon_b`` lies to the west of ``lon_a`` (the direction easterly
    waves travel). NaN in either argument propagates.

    >>> lon_delta(-179.0, 179.0)
    2.0
    >>> lon_delta(179.0, -179.0)
    -2.0
    """
    out = (np.asarray(lon_b, dtype=float) - np.asarray(lon_a, dtype=float) + 180.0) % 360.0 - 180.0
    return _as_float(out, _scalar_in(lon_a, lon_b))


def is_westward(lon_new, lon_old):
    """True when ``lon_new`` is at or west of ``lon_old`` by the shortest path.

    Drop-in replacement for the old ``lon_new <= lon_old`` test. A stationary
    point counts as westward, matching the original behaviour, and NaN gives
    False (also matching, since NaN comparisons are False).
    """
    delta = lon_delta(lon_new, lon_old)
    with np.errstate(invalid="ignore"):
        out = np.asarray(delta) <= 0.0
    return _as_bool(out, _scalar_in(lon_new, lon_old))


def lon_in_range(lon, west, east, inclusive="both"):
    """Modular test for ``lon`` lying in the window running east from ``west`` to ``east``.

    The window is always traversed *eastward* from ``west``, so ``(170, -170)``
    is the 20-degree window straddling the dateline, while ``(-170, 170)`` is the
    340-degree window that does not. ``west == east`` means the whole globe.

    ``inclusive`` is one of ``"both"``, ``"neither"``, ``"west"`` or ``"east"``
    and controls the endpoints, so the caller can reproduce the exact ``<``/``<=``
    semantics of the code being replaced.
    """
    scalar = _scalar_in(lon)
    lon = np.asarray(lon, dtype=float)
    width = (float(east) - float(west)) % 360.0
    if width == 0.0:
        width = 360.0
    offset = (lon - float(west)) % 360.0

    # A point exactly on ``west`` has offset 0.0; one exactly on ``east`` has
    # offset ``width`` (except for the full-globe case, where they coincide).
    with np.errstate(invalid="ignore"):
        if width == 360.0:
            inside = np.isfinite(lon)
        elif inclusive == "both":
            inside = offset <= width
        elif inclusive == "neither":
            inside = (offset > 0.0) & (offset < width)
        elif inclusive == "west":
            inside = offset < width
        elif inclusive == "east":
            inside = (offset > 0.0) & (offset <= width)
        else:
            raise ValueError(f"inclusive must be both/neither/west/east, got {inclusive!r}")
        inside = inside & np.isfinite(lon)
    return _as_bool(inside, scalar)


def lon_in_any_range(lon, ranges, inclusive="both"):
    """True where ``lon`` falls inside any of a list of ``(west, east)`` windows.

    ``ranges`` of None means "no restriction" and returns True everywhere
    (except at NaN). A single ``(west, east)`` tuple is also accepted.
    """
    scalar = _scalar_in(lon)
    lon_arr = np.asarray(lon, dtype=float)
    if ranges is None:
        return _as_bool(np.isfinite(lon_arr), scalar)
    if len(ranges) == 2 and np.ndim(ranges[0]) == 0:
        ranges = [ranges]
    out = np.zeros(lon_arr.shape, dtype=bool)
    for west, east in ranges:
        out = out | np.asarray(lon_in_range(lon_arr, west, east, inclusive=inclusive))
    return _as_bool(out, scalar)


def unwrap_lon(lon_series):
    """Return a continuous longitude series, removing 360-degree jumps.

    NaN-safe: gaps are skipped rather than terminating the series, so a track
    that loses a timestep still unwraps across the gap. The first finite value is
    kept as-is, so ``wrap180(unwrap_lon(x)) == wrap180(x)``.

    This is what makes ``np.interp``, ``savgol_filter`` and ``LinearRegression``
    usable on a track that crosses the seam.
    """
    lon = np.asarray(lon_series, dtype=float)
    if lon.ndim != 1:
        raise ValueError("unwrap_lon expects a 1-D series; apply it row by row")
    out = np.full(lon.shape, np.nan)
    valid = np.isfinite(lon)
    if not valid.any():
        return out
    values = lon[valid]
    if values.size == 1:
        out[valid] = values
        return out
    steps = wrap180(np.diff(values))
    out[valid] = values[0] + np.concatenate(([0.0], np.cumsum(steps)))
    return out


def on_unwrapped(func, lon_series, *args, **kwargs):
    """Run ``func`` on an unwrapped copy of ``lon_series`` and re-wrap the result.

    Used for every operation that treats a track's longitude history as a smooth
    1-D signal (gap filling, Savitzky-Golay smoothing, linear extrapolation).
    """
    return wrap180(func(unwrap_lon(lon_series), *args, **kwargs))


def grid_resolution(lon):
    """Median absolute spacing of a longitude axis, in degrees."""
    lon = np.asarray(lon, dtype=float)
    if lon.size < 2:
        return 0.0
    return float(np.median(np.abs(wrap180(np.diff(lon)))))


def is_global_lon(lon, res=None):
    """True when a longitude axis spans the globe and can be treated as periodic.

    The axis need not start at any particular longitude; what matters is that the
    span plus one grid spacing closes the circle.
    """
    lon = np.asarray(lon, dtype=float)
    if lon.size < 2:
        return False
    if res is None:
        res = grid_resolution(lon)
    if res <= 0:
        return False
    span = float(lon.max() - lon.min()) + res
    return span >= 360.0 - 0.5 * res


def find_nearest_lon(lon_axis, value):
    """Wrap-aware nearest-gridpoint lookup on a longitude axis.

    Same ``(index, value)`` return as the tracker's ``find_nearest``, but a query
    of -181 finds +179 instead of clamping to the first column.
    """
    lon_axis = np.asarray(lon_axis, dtype=float)
    idx = int(np.abs(lon_delta(lon_axis, value)).argmin())
    return idx, float(lon_axis[idx])


def pad_lon(arr, n, axis=-1):
    """Circularly pad an array by ``n`` columns on each side of the longitude axis."""
    if n <= 0:
        return np.asarray(arr)
    arr = np.asarray(arr)
    pad_width = [(0, 0)] * arr.ndim
    pad_width[axis] = (n, n)
    return np.pad(arr, pad_width, mode="wrap")


def unpad_lon(arr, n, axis=-1):
    """Inverse of :func:`pad_lon` -- strip ``n`` columns from each end."""
    if n <= 0:
        return np.asarray(arr)
    arr = np.asarray(arr)
    slicer = [slice(None)] * arr.ndim
    slicer[axis] = slice(n, arr.shape[axis] - n)
    return arr[tuple(slicer)]


def pad_lon_coord(lon, n):
    """Circularly pad a longitude *coordinate* with continuous values.

    The coordinate must stay monotonic across the pad (..., 178, 179, 180, 181,
    ...), not wrap back to -180. Padding it with raw wrapped values would put a
    -360 spike into ``np.gradient(lon)`` and therefore into every ``dx`` derived
    from it.
    """
    lon = np.asarray(lon, dtype=float)
    if n <= 0:
        return lon.copy()
    return np.concatenate((lon[-n:] - 360.0, lon, lon[:n] + 360.0))


def wrap_index_delta(delta, n):
    """Signed shortest separation between two indices on a periodic axis of length ``n``.

    ``wrap_index_delta(1 - 359, 360) == 2``: index 1 is two gridpoints east of
    index 359, not 358 gridpoints west.
    """
    scalar = _scalar_in(delta)
    out = (np.asarray(delta) + n // 2) % n - n // 2
    return int(out) if scalar else out
