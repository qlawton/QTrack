"""End-to-end proof that a wave is tracked continuously across the dateline.

A single cyclonic vortex is translated due west at 7 m/s through 180 for ten
days and pushed through the whole pipeline. Before this branch the track died in
the band of zero curvature vorticity around the seam; even with that filled it
broke into two pieces, because the tracker's direction test read the crossing as
a 359 degree jump east.

Everything here is synthetic, so the test needs no network access.
"""

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
    with xr.open_dataset(post) as ds:
        return ds.load(), true_lon


def _only_track(ds):
    lon = ds["AEW_lon"].values
    keep = ~np.isnan(lon).all(axis=1)
    assert keep.sum() == 1, f"expected one track, got {keep.sum()}"
    row = int(np.flatnonzero(keep)[0])
    return lon[row], ds["AEW_lat"].values[row], ds["AEW_lon_unwrapped"].values[row]


def test_single_unbroken_track(crossing):
    ds, _ = crossing
    lon, _, _ = _only_track(ds)
    valid = np.flatnonzero(~np.isnan(lon))
    assert valid.size == N_TIMES, f"track covers {valid.size} of {N_TIMES} timesteps"
    assert np.array_equal(valid, np.arange(valid[0], valid[-1] + 1)), "track has a gap"


def test_no_gap_at_the_seam(crossing):
    """The old failure mode was losing the wave within about 20 degrees of 180."""
    ds, _ = crossing
    lon, _, _ = _only_track(ds)
    near_seam = np.abs(np.abs(lon) - 180.0) <= 20.0
    assert near_seam.sum() >= 4, "the wave never got near the dateline; test is not exercising the seam"
    assert not np.isnan(lon[near_seam]).any()


def test_unwrapped_longitude_is_monotonically_westward(crossing):
    ds, _ = crossing
    _, _, unwrapped = _only_track(ds)
    finite = unwrapped[~np.isnan(unwrapped)]
    assert np.all(np.diff(finite) < 0), "an easterly wave must move west at every step"
    # It really did cross: the unwrapped series leaves the -180..180 interval.
    assert finite.min() < -180.0


def test_track_follows_the_true_position(crossing):
    ds, true_lon = crossing
    lon, lat, _ = _only_track(ds)
    error = np.abs(qtrack.geo.lon_delta(lon, true_lon))
    assert np.nanmax(error) <= 2.0, f"max longitude error {np.nanmax(error):.2f} deg"
    assert np.nanmax(np.abs(lat - CENTER_LAT)) <= 2.0


def test_total_westward_travel_matches_the_prescribed_speed(crossing):
    ds, true_lon = crossing
    _, _, unwrapped = _only_track(ds)
    finite = unwrapped[~np.isnan(unwrapped)]
    travelled = finite[0] - finite[-1]
    expected = true_lon[0] - true_lon[-1]
    assert travelled == pytest.approx(expected, abs=2.0)
