"""End-to-end proof that a wave is tracked continuously across the dateline.

A single cyclonic vortex is translated due west at 7 m/s through 180 for ten
days and pushed through the whole pipeline. Before this branch the track died in
the band of zero curvature vorticity around the seam; even with that filled it
broke into two pieces, because the tracker's direction test read the crossing as
a 359 degree jump east.

Everything here is synthetic, so the test needs no network access.
"""

import figures
import numpy as np
import pytest
import xarray as xr
from synthetic import westward_vortex_dataset

import qtrack
from qtrack.curvvort import compute_curvvort
from qtrack.tracking import run_postprocessing, run_tracking

N_TIMES = 40
START_LON = -150.0
CENTER_LAT = 12.0


@pytest.fixture(scope="module")
def crossing(tmp_path_factory):
    """Run the full pipeline once and hand back the tracks plus the true positions."""
    tmp_path = tmp_path_factory.mktemp("dateline")
    raw = tmp_path / "winds.nc"
    prepped = tmp_path / "prepped.nc"
    cv = tmp_path / "cv.nc"
    tracks = tmp_path / "tracks_raw.nc"
    post = tmp_path / "tracks_post.nc"

    true_lon = westward_vortex_dataset(str(raw), n_times=N_TIMES, start_lon=START_LON, center_lat=CENTER_LAT)

    qtrack.prep_data(str(raw), data_out=str(prepped))
    compute_curvvort(str(prepped), str(cv), njobs_in=1)
    # initiation_bounds=None because this wave is born in the Pacific, not over Africa.
    run_tracking(input_file=str(cv), save_file=str(tracks), initiation_bounds=None)
    run_postprocessing(
        input_file=str(tracks),
        curv_data_file=str(cv),
        real_year_used=2010,
        save_obj_file=str(tmp_path / "tracks.pkl"),
        save_nc_file=str(post),
        hovmoller_save=False,
    )
    with xr.open_dataset(post) as ds, xr.open_dataset(tracks) as raw_ds:
        return ds.load(), raw_ds.load(), true_lon


def _only_track(ds):
    lon = ds["AEW_lon"].values
    keep = ~np.isnan(lon).all(axis=1)
    assert keep.sum() == 1, f"expected one track, got {keep.sum()}"
    row = int(np.flatnonzero(keep)[0])
    return lon[row], ds["AEW_lat"].values[row], ds["AEW_lon_unwrapped"].values[row]


def test_single_unbroken_track(crossing):
    ds, _, _ = crossing
    lon, _, _ = _only_track(ds)
    valid = np.flatnonzero(~np.isnan(lon))
    assert valid.size == N_TIMES, f"track covers {valid.size} of {N_TIMES} timesteps"
    assert np.array_equal(valid, np.arange(valid[0], valid[-1] + 1)), "track has a gap"


def test_no_gap_at_the_seam(crossing):
    """The old failure mode was losing the wave within about 20 degrees of 180."""
    ds, _, _ = crossing
    lon, _, _ = _only_track(ds)
    near_seam = np.abs(np.abs(lon) - 180.0) <= 20.0
    assert near_seam.sum() >= 4, "the wave never got near the dateline; test is not exercising the seam"
    assert not np.isnan(lon[near_seam]).any()


def test_unwrapped_longitude_is_monotonically_westward(crossing):
    ds, _, _ = crossing
    _, _, unwrapped = _only_track(ds)
    finite = unwrapped[~np.isnan(unwrapped)]
    assert np.all(np.diff(finite) < 0), "an easterly wave must move west at every step"
    # It really did cross: the unwrapped series leaves the -180..180 interval.
    assert finite.min() < -180.0


def test_track_follows_the_true_position(crossing):
    ds, _, true_lon = crossing
    lon, lat, _ = _only_track(ds)
    error = np.abs(qtrack.geo.lon_delta(lon, true_lon))
    assert np.nanmax(error) <= 2.0, f"max longitude error {np.nanmax(error):.2f} deg"
    assert np.nanmax(np.abs(lat - CENTER_LAT)) <= 2.0


def test_total_westward_travel_matches_the_prescribed_speed(crossing):
    ds, _, true_lon = crossing
    _, _, unwrapped = _only_track(ds)
    finite = unwrapped[~np.isnan(unwrapped)]
    travelled = finite[0] - finite[-1]
    expected = true_lon[0] - true_lon[-1]
    assert travelled == pytest.approx(expected, abs=2.0)


def test_short_tracks_are_still_removed(crossing):
    """run_tracking must drop tracks shorter than AEW_day_remove.

    The check reads the track length back out of the row after the gap fill has
    restored the leading and trailing NaNs. An earlier version of the wrap-safe
    gap fill rebound the working array, so the length was measured on the
    fully-interpolated series and nothing was ever short enough to drop -- five
    spurious one-point tracks survived into the output.
    """
    _, raw_ds, _ = crossing
    lengths = (~np.isnan(raw_ds["AEW_lon"].values)).sum(axis=1)
    cutoff = 2 * 24 // 6  # AEW_day_remove=2 days at 6-hourly resolution
    assert lengths.min() > cutoff, f"raw output kept tracks of length {sorted(lengths.tolist())}"


def test_figure_dateline_crossing(crossing):
    """Not an assertion -- draws the crossing the assertions above check."""
    if figures.figure_dir() is None:
        pytest.skip(f"set {figures.ENV_VAR} to write diagnostic figures")
    ds, raw_ds, true_lon = crossing
    lon, lat, unwrapped = _only_track(ds)
    plt = figures.pyplot()

    cv_lon = raw_ds["longitude"].values
    cv = raw_ds["curv_data_mean"].values * 1e6
    steps = np.arange(len(lon))
    true_wrapped = (true_lon + 180.0) % 360.0 - 180.0

    fig = plt.figure(figsize=(12.5, 8.6))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.55, 1.0], hspace=0.38, wspace=0.24)

    # -- wrapped Hovmoller: the track leaves one edge and re-enters the other ----
    ax = fig.add_subplot(grid[0, 0])
    mesh = figures.hovmoller(ax, cv_lon, cv)
    figures.shade_old_halo(ax)
    figures.plot_track(ax, true_wrapped, steps, color=figures.TRUTH_INK, linewidth=2.6, halo=False)
    figures.plot_track(ax, true_wrapped, steps, color=figures.TRUTH_HALO, linewidth=1.0, linestyle=(0, (4, 3)), halo=False)
    figures.plot_track(ax, lon, steps, linewidth=1.9)
    ax.set_xlim(-180, 180)
    ax.invert_yaxis()
    ax.set_title("As stored: -180 to 180", fontsize=11, loc="left")
    ax.text(-176, 1, "shaded: the old zero halo", fontsize=8.5, color="#333333", va="top")
    fig.colorbar(mesh, ax=ax, pad=0.015, aspect=30, label="5-20N mean curv. vort. (1e-6 s-1)")

    # -- unwrapped: one straight line ------------------------------------------
    ax = fig.add_subplot(grid[0, 1])
    figures.plot_track(ax, unwrapped, steps, linewidth=1.9)
    ax.plot(true_lon, steps, color="#c1440e", linewidth=1.1, linestyle=(0, (5, 3)), zorder=7)
    ax.axvline(-180, color="#111111", linewidth=0.8, linestyle=(0, (4, 3)), alpha=0.7)
    ax.text(-181, 1, " crosses 180 here", fontsize=8.5, color="#333333", rotation=90, va="top", ha="right")
    ax.set_xlabel("unwrapped longitude (deg)")
    ax.set_ylabel("timestep (6-hourly)")
    ax.invert_yaxis()
    ax.set_title("AEW_lon_unwrapped: continuous", fontsize=11, loc="left")
    figures.tidy(ax)

    # -- latitude and error ------------------------------------------------------
    ax = fig.add_subplot(grid[1, 0])
    ax.plot(steps, lat, color=figures.TRACK_INK, linewidth=1.8, label="tracked")
    ax.axhline(CENTER_LAT, color="#c1440e", linewidth=1.1, linestyle=(0, (5, 3)), label="prescribed")
    ax.set_xlabel("timestep (6-hourly)")
    ax.set_ylabel("latitude (deg)")
    ax.set_title("Latitude", fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9)
    figures.tidy(ax)

    ax = fig.add_subplot(grid[1, 1])
    error = np.abs(qtrack.geo.lon_delta(lon, true_lon))
    crossing_step = int(np.argmax(np.abs(np.diff(lon)) > 180.0)) + 1
    ax.plot(steps, error, color=figures.TRACK_INK, linewidth=1.8)
    ax.axvline(crossing_step, color="#c1440e", linewidth=1.1, linestyle=(0, (5, 3)))
    ax.text(crossing_step + 0.5, ax.get_ylim()[1] * 0.92, "seam crossing", fontsize=8.5, color="#c1440e")
    ax.set_xlabel("timestep (6-hourly)")
    ax.set_ylabel("|longitude error| (deg)")
    ax.set_title(f"Error against the prescribed track -- max {np.nanmax(error):.2f} deg, no jump at the seam", fontsize=11, loc="left")
    figures.tidy(ax)

    fig.suptitle("Synthetic vortex translating west at 7 m/s through the dateline, tracked end to end", fontsize=13, y=0.965)
    figures.save(fig, "02_dateline_crossing")
