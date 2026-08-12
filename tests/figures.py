"""Opt-in diagnostic figures from the test suite.

Set ``QTRACK_TEST_FIGDIR`` to a directory and the tests write PNGs there:

    QTRACK_TEST_FIGDIR=./figs pytest tests/test_dateline.py tests/test_curvvort_seam.py

Unset (the default, and in CI) nothing is plotted and the tests are unchanged.
The figures are for eyeballing the seam behaviour; the assertions, not the
pictures, are what actually pass or fail.

Plotting conventions used throughout:

* Curvature vorticity is signed, so it gets a diverging colormap with a neutral
  midpoint, symmetric limits, and never a rainbow.
* Track overlays are drawn in ink (black / white / grey) with a contrasting halo
  rather than in hue, so they stay legible on top of the diverging field and
  carry no colour-only meaning.
* Where several tracks share a panel, they are distinguished by weight and a
  direct label, not by twenty-odd cycled hues.
"""

import os
from pathlib import Path

import numpy as np

ENV_VAR = "QTRACK_TEST_FIGDIR"

# Ink used for overlays on top of the diverging field.
TRACK_INK = "#111111"
TRACK_HALO = "#ffffff"
TRUTH_INK = "#ffffff"
TRUTH_HALO = "#111111"
MUTED_INK = "#5a5a5a"


def figure_dir():
    """Directory to write figures to, or None when figure output is off."""
    target = os.environ.get(ENV_VAR)
    if not target:
        return None
    path = Path(target)
    path.mkdir(parents=True, exist_ok=True)
    return path


def pyplot():
    """Import pyplot with a headless backend, only when actually plotting."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def save(fig, name):
    """Write ``fig`` into the figure directory and close it."""
    path = figure_dir() / f"{name}.png"
    fig.savefig(path, dpi=140, bbox_inches="tight", facecolor="white")
    pyplot().close(fig)
    print(f"[figure] {path}")
    return path


def symmetric_limits(field, percentile=99.5):
    """Colour limits centred on zero, robust to a few extreme gridpoints."""
    scale = np.nanpercentile(np.abs(field), percentile)
    if not np.isfinite(scale) or scale == 0:
        scale = 1.0
    return -scale, scale


def tidy(ax):
    """Recessive axes: no top/right spines, a faint grid behind the data."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#bbbbbb")
    ax.tick_params(colors="#555555", labelsize=9)
    ax.grid(True, color="#000000", alpha=0.08, linewidth=0.6)
    ax.set_axisbelow(True)


def mark_dateline(ax, label=True):
    """Draw the array seam so it is obvious where the old code fell apart."""
    ax.axvline(180, color="#111111", linewidth=0.8, linestyle=(0, (4, 3)), alpha=0.7)
    ax.axvline(-180, color="#111111", linewidth=0.8, linestyle=(0, (4, 3)), alpha=0.7)
    if label:
        ax.text(178, ax.get_ylim()[1], " dateline", ha="right", va="top", rotation=90, fontsize=8, color="#333333")


def shade_old_halo(ax, halo_deg=12.0):
    """Shade the longitudes GetBG used to leave at exactly zero on a global grid."""
    for x0, x1 in ((-180, -180 + halo_deg), (180 - halo_deg, 180)):
        ax.axvspan(x0, x1, color="#111111", alpha=0.10, linewidth=0, zorder=3)


def recentre_lon(lon, centre):
    """Express longitudes in the 360-degree window centred on ``centre``.

    With ``centre=180`` a track running 179 -> -179 becomes 179 -> 181, i.e.
    continuous, which is what a map centred on the dateline needs. With
    ``centre=0`` it stays 179 -> -179 and :func:`plot_track` breaks it at the
    frame edge instead of streaking it across the map.
    """
    lon = np.asarray(lon, dtype=float)
    return centre + ((lon - centre + 180.0) % 360.0 - 180.0)


def plot_track(ax, lon, y, halo=True, centre=None, **kwargs):
    """Overlay a track, splitting it where it wraps so no false line crosses the panel."""
    style = {"color": TRACK_INK, "linewidth": 1.8, "zorder": 6}
    style.update(kwargs)
    lon = np.asarray(lon, dtype=float)
    if centre is not None:
        lon = recentre_lon(lon, centre)
    y = np.asarray(y, dtype=float)
    # Break the line wherever consecutive points jump more than half the globe:
    # that is a wrap, not motion, and drawing through it would be a lie.
    # Anything the caller passed that is not a line-appearance override -- notably
    # a cartopy `transform` -- has to reach the halo as well, or the halo is drawn
    # in a different coordinate system from the line it is supposed to sit under.
    passthrough = {k: v for k, v in style.items() if k not in ("color", "linewidth", "zorder", "linestyle", "label")}
    breaks = np.flatnonzero(np.abs(np.diff(lon)) > 180.0) + 1
    for seg_lon, seg_y in zip(np.split(lon, breaks), np.split(y, breaks)):
        if np.isfinite(seg_lon).sum() < 1:
            continue
        if halo:
            ax.plot(seg_lon, seg_y, color=TRACK_HALO, linewidth=style["linewidth"] + 1.6, solid_capstyle="round", zorder=style["zorder"] - 1, **passthrough)
        ax.plot(seg_lon, seg_y, **style)


def global_map(fig, subplot_spec, lon, lat, field, extent, centre=0.0, cmap="RdBu_r", limits=None):
    """Map of a curvature vorticity field with coastlines, on a plate carree grid.

    ``centre`` is the central longitude of the projection: 0 for the conventional
    Africa-and-Atlantic view, 180 to put the dateline in the middle of the frame.
    Returns ``(ax, mesh)``.
    """
    import cartopy.crs as ccrs

    proj = ccrs.PlateCarree(central_longitude=centre)
    spec = subplot_spec if isinstance(subplot_spec, tuple) else (subplot_spec,)
    ax = fig.add_subplot(*spec, projection=proj)
    vmin, vmax = limits if limits is not None else symmetric_limits(field)
    mesh = ax.pcolormesh(lon, lat, field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True, transform=ccrs.PlateCarree())
    ax.coastlines("110m", linewidth=0.6, color="#3a3a3a")
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    grid = ax.gridlines(draw_labels=True, linewidth=0.5, color="#000000", alpha=0.12)
    grid.top_labels = False
    grid.right_labels = False
    grid.xlabel_style = {"size": 8, "color": "#555555"}
    grid.ylabel_style = {"size": 8, "color": "#555555"}
    return ax, mesh


def mark_detection_bands(ax, lat_bounds=(5, 20)):
    """Show the latitude band the meridional averaging step actually searches."""
    import cartopy.crs as ccrs

    for value in lat_bounds:
        ax.plot([-180, 180], [value, value], color="#111111", linewidth=0.7, linestyle=(0, (5, 4)), alpha=0.55, transform=ccrs.PlateCarree(), zorder=4)


def hovmoller(ax, lon, field, cmap="RdBu_r"):
    """Longitude-time shading of the meridionally averaged curvature vorticity."""
    vmin, vmax = symmetric_limits(field)
    mesh = ax.pcolormesh(lon, np.arange(field.shape[0]), field, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
    ax.set_xlabel("longitude (deg)")
    ax.set_ylabel("timestep (6-hourly)")
    tidy(ax)
    return mesh
