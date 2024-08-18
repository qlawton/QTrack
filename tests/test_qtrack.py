import os

import pytest

import qtrack

# from qtrack.curvvort import compute_curvvort
# from qtrack.tracking import run_postprocessing, run_tracking


def test_bad_example_id_raises():
    with pytest.raises(ValueError):
        qtrack.download_examples("bad_example_id")


def test_small_tracking(tmp_path):
    # Move to temp dir and download example data
    os.chdir(tmp_path)
    example_id = "era5_2010_10day"
    example_fn = "era5_700_wind_global_2010_10day.nc"
    qtrack.download_examples(example_id)
    assert os.path.exists(example_fn)

    # Prep
    prepped_fn = "prepped_for_tracking.nc"
    qtrack.prep_data(data_in=example_fn, data_out=prepped_fn)

    # # Compute curvature vorticity
    # cv_fn = "cv.nc"
    # compute_curvvort(data_in=prepped_fn, data_out=cv_fn, njobs_in=-1)

    # # Track
    # raw_fn = "tracks_raw.nc"
    # run_tracking(input_file=cv_fn, save_file=raw_fn)

    # # AEW post-processing
    # tracks_fn = "tracks.nc"
    # tracks_obj_fn = "tracks.pkl"
    # run_postprocessing(
    #     input_file=raw_fn,
    #     real_year_used=2010,
    #     curv_data_file=cv_fn,
    #     save_obj_file=tracks_obj_fn,
    #     save_nc_file=tracks_fn,
    # )
