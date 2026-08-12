"""The radially averaged curvature vorticity must not know where the array seam is.

``GetBG`` used to leave its longitude halo at exactly 0.0. On a global grid that
is a band of roughly 11 columns of fake zero curvature vorticity straddling the
dateline, so any wave reaching 180 was erased before the tracker ever saw it.
"""

import figures
import numpy as np
import pytest
import xarray as xr
from synthetic import westward_vortex_dataset

import qtrack
from qtrack.curvvort import compute_curvvort


def _curv_vort_for_vortex_at(tmp_path, lon_c, tag):
    """Run prep + curvature vorticity for a stationary vortex centred on ``lon_c``."""
    raw = tmp_path / f"raw_{tag}.nc"
    prepped = tmp_path / f"prepped_{tag}.nc"
    cv = tmp_path / f"cv_{tag}.nc"
    westward_vortex_dataset(str(raw), n_times=1, start_lon=lon_c, speed_ms=0.0)
    qtrack.prep_data(str(raw), data_out=str(prepped))
    compute_curvvort(str(prepped), str(cv), njobs_in=1)
    with xr.open_dataset(cv) as ds:
        return ds["curv_vort"].values[0], ds["longitude"].values, ds["latitude"].values


@pytest.fixture(scope="module")
def seam_and_centre(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("cvseam")
    on_seam = _curv_vort_for_vortex_at(tmp_path, -180.0, "seam")
    at_centre = _curv_vort_for_vortex_at(tmp_path, 0.0, "centre")
    return on_seam, at_centre


def test_no_all_zero_longitude_columns(seam_and_centre):
    (field, _, _), _ = seam_and_centre
    zero_cols = np.flatnonzero((field == 0.0).all(axis=0))
    assert zero_cols.size == 0, f"{zero_cols.size} longitude columns are identically zero"


def test_latitude_halo_is_still_unfilled(seam_and_centre):
    # Only longitude was made periodic; the latitude halo behaves as before, and
    # sits far outside the 5-20N tracking bands.
    (field, _, lat), _ = seam_and_centre
    zero_rows = np.flatnonzero((field == 0.0).all(axis=1))
    assert zero_rows.size > 0
    # ...and outside the 5-20N band the tracker averages over.
    assert np.all((lat[zero_rows] < 5.0) | (lat[zero_rows] > 20.0))


def test_vortex_on_the_seam_is_detected(seam_and_centre):
    (field, lon, lat), _ = seam_and_centre
    i, j = np.unravel_index(np.nanargmax(field), field.shape)
    assert abs(lon[j] - -180.0) <= 2.0 or abs(lon[j] - 179.0) <= 2.0
    assert abs(lat[i] - 12.0) <= 2.0


def test_result_is_independent_of_seam_placement(seam_and_centre):
    """A vortex on the dateline and the same vortex at 0E must look identical.

    This is the invariance the old code broke: physics cannot depend on which
    column numpy happens to call zero.
    """
    (seam_field, seam_lon, _), (centre_field, centre_lon, _) = seam_and_centre

    # Roll the on-seam result so its vortex sits where the 0E vortex does.
    shift = int(round(180.0 / (centre_lon[1] - centre_lon[0])))
    rolled = np.roll(seam_field, shift, axis=1)

    valid = (rolled != 0.0) & (centre_field != 0.0)
    assert valid.any()
    np.testing.assert_allclose(rolled[valid], centre_field[valid], rtol=1e-6, atol=1e-12)

    assert np.nanmax(seam_field) == pytest.approx(np.nanmax(centre_field), rel=1e-6)


def test_figure_seam_invariance(seam_and_centre):
    """Not an assertion -- draws the evidence for the invariance tested above."""
    if figures.figure_dir() is None:
        pytest.skip(f"set {figures.ENV_VAR} to write diagnostic figures")
    (seam_field, lon, lat), (centre_field, _, _) = seam_and_centre
    plt = figures.pyplot()

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 10.5), gridspec_kw={"height_ratios": [1, 1, 0.85], "hspace": 0.42})
    vmin, vmax = figures.symmetric_limits(centre_field)

    for ax, field, title in (
        (axes[0], seam_field * 1e6, "Vortex centred on the dateline (180)"),
        (axes[1], centre_field * 1e6, "The same vortex centred on 0E"),
    ):
        mesh = ax.pcolormesh(lon, lat, field, cmap="RdBu_r", vmin=vmin * 1e6, vmax=vmax * 1e6, shading="auto", rasterized=True)
        figures.shade_old_halo(ax)
        figures.mark_dateline(ax, label=False)
        ax.set_title(title, fontsize=11, loc="left")
        ax.set_xlabel("longitude (deg)")
        ax.set_ylabel("latitude (deg)")
        ax.set_xlim(-180, 180)
        figures.tidy(ax)
        fig.colorbar(mesh, ax=ax, pad=0.015, aspect=28, label="curv. vort. (1e-6 s-1)")

    axes[0].text(
        -179,
        lat.min() + 1.5,
        "shaded bands: the ~12 columns GetBG used to leave at exactly 0.0",
        fontsize=8.5,
        color="#333333",
        va="bottom",
    )

    # Profiles through the vortex latitude, each re-centred on its own vortex, so
    # the two curves are directly comparable.
    row = int(np.argmin(np.abs(lat - 12.0)))
    offset = np.arange(lon.size) - lon.size // 2
    seam_profile = np.roll(seam_field[row], lon.size // 2) * 1e6
    centre_profile = np.roll(centre_field[row], lon.size // 2 - int(np.argmin(np.abs(lon - 0.0)))) * 1e6

    ax = axes[2]
    ax.plot(offset, seam_profile, color=figures.TRACK_INK, linewidth=2.0, label="vortex on 180")
    ax.plot(offset, centre_profile, color="#c1440e", linewidth=2.0, linestyle=(0, (5, 3)), label="vortex on 0E")
    ax.set_xlim(-40, 40)
    ax.set_xlabel("degrees longitude from the vortex centre")
    ax.set_ylabel("curv. vort. (1e-6 s-1)")
    ax.set_title(f"Zonal profile at {lat[row]:.0f}N -- max |difference| = {np.nanmax(np.abs(seam_profile - centre_profile)):.2e}", fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9)
    figures.tidy(ax)

    fig.suptitle("Radially averaged curvature vorticity does not depend on where the array seam is", fontsize=13, y=0.945)
    figures.save(fig, "01_curvvort_seam_invariance")
