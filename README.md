# General Information
This repo contains code to run the AEW Tracker described by Lawton et al. (2022), which can be found here: https://doi.org/10.1175/MWR-D-21-0321.1.

The AEW tracker was developed by Quinton Lawton, currently affiliated with the NSF National Center for Atmospheric Research (NSF NCAR). This python module was developed with the support of Zachary Moon (NOAA ARL/Texas A&M University) and Kelly Núñez Ocasio (Texas A&M University).  

Existing AEW tracks for ERA5 data from 1979 - 2022 are already produced and available here: https://osf.io/j4hpq/

## How the tracker works
The AEW tracker uses an input file containing 700hPa zonal (u) and meridional (v) wind and computes a curvature vorticity (CV) term that is smoothed using a 600km average around each gridpoint. Local maxima in this CV are then tracked in two different ways. Over the African continent, meridional averages of CV are used to find longitudinal local maxima, which are used as a first guess for a centroid location. Set thresholds of CV value define when a new AEW is initialized and whether an existing one continues being tracked. Over the Atlantic, the first guess for an existing AEW's location is estimated by the AEW's propagation speed, and it is later found again using a centroid. A post processing script cleans up some of the tracks.

More detailed information on the tracking, not provided in the Lawton et al. (2022) paper or here, can be found in an online technical guide here: https://osf.io/6hqy5

## Important Note on Non-Divergent Wind Step
Due to inconsistencies in the windspharm package, which computes the non-divergent wind using spherical harmonics, the non-divergent wind step is currently not included in the package. This step was largely superfluous and thus not including this step is is not anticipated to have any major negative impacts. Nevertheless, it could result in AEW tracks slightly differing from those produced with this step included, including the AEW tracks in the Lawton et al. (2022) AEW databases. We hope to include a non-divergent wind step in a future release.

## Input files
**You need u and v winds on a 1x1 degree grid with 6-hourly temporal outputs**. This is the optimal resolution tested in Lawton et al. (2022), but it can be adjusted manually in the code if necessary to change. Furthermore, it is highly recommended that at least 10 days of data are included. The tracking needs a bit of spinup, so if you are running this on model output, it is recommended you append a few days to 1 week of analysis/reanalysis data prior to the first model timestep. Input data should geographically cover at least a portion of the African continent and/or Atlantic Ocean. Note that by default, the tracker does not initiate new waves west of 35W, but this can be adjusted in the tracking function.

A sample model data/reanalysis combination script is provided under `tools/wind_combine.py`

## Output AEW Track files
Two different kinds of output files are produced, a netCDF file and a python pickle file (with AEWs saved as objects). It is recommended that you use netCDF if possible. However, if you need to use the pickle files for whatever reason, you will need to have the `AEW_module.py` script available locally when you open them.

Details on the AEW objects and AEW_module.py are available here: https://osf.io/jnv5u

## Citing
If you use this AEW tracker in your research, we ask that you acknowledge your use of this package and provide a citation to the original paper documenting the tracker.

Lawton, Q. A., S. J. Majumdar, K. Dotterer, C. Thorncroft, and C. J. Schreck, 2022: The Influence of Convectively Coupled Kelvin Waves on African Easterly Waves in a Wave-Following Framework. Monthly Weather Review, 150(8), 2055-2072,  https://doi.org/10.1175/MWR-D-21-0321.1.

# Getting Started

## Recommended Workflow Steps
- Adjust input data to 1x1 degree resolution and 6 hourly temporal resolution (must be done prior to running this code)
- Prepare input data using `prep_data` function
- Compute curvature vorticity fields using `curvvort` function
- Run the AEW tracking code using `run_tracking` function
- Run post processing of AEWs using `run_postprocessing` function

## Example Workflow
It is highly recommended that users view the `example_tracking.ipynb` or `example_tracking.py` files to see an example of AEW tracking.

## Example Files
The `track.download_examples` function can be used to download datasets for users to test out AEW tracking on. Possible datasets to download include:
- `mpas_2021092400`: Downsampled data from a 30km run of the MPAS-A system for 9-24-2021 at 00Z. This data originates from research by Lawton et al. (2024): https://doi.org/10.1029/2023MS004187
- `era5_2010`: ERA5 data for July through September of 2010
- `era5_2010_10day`: A compressed 10 day period in September of 2010 (for a quicker running example)
- `gfs_2024062612`: GFS run from June 26th, 2024 at 12Z. Analysis data was appended before model start to allow for "spin up" time.

## How to import required functions
~~~~
import qtrack
from qtrack.curvvort import compute_curvvort
from qtrack.tracking import run_postprocessing, run_tracking
~~~~

# Details on essential included functions

## **Prep Data** `qtrack.prep_data(data_in, cut_lev_val=700, data_out="prepped_data_for_tracking.nc")`
This function ensures the AEW tracker runs properly on the included data. Data should already be at 1x1 degree resolution and 6 hourly temporal resolution. This script will flip the lon/lat coordinates to their property directions, check that all necessary components are included, and slice out the desired height-level if it's 3D data.
- data_in: The input data. Should be 2- or 3-dimensional data with U-wind component labeled as "u" and V-wind component labelled as "v". The data file should contain longitudes named "longitude", "lon", or "lons", and latitudes labeled "latitude", "lat", or "lats'.
- cut_lev_val (Default: 700): If 3D data, the pressure level to cut out for AEW tracking. Recommended to be near 700hPa.
- data_out (Default: "prepped_data_for_tracking.nc"): The name of the output data.

## **Compute Curvature Vorticity and Averaging** `qtrack.curvvort.compute_curvvort(data_in, data_out="radial_avg_curv_vort.nc", radius_of_avg=600, data_resolution=1, njobs_in=1, nondiv_wind=False, run_animation=False, gif_dir_in="")`
This step computes the CV and takes the gridpoint averages. It is slow. **However, there is an option within the script to use multiprocessing, which is highly recommended.**

## **AEW Tracking** `qtrack.tracking.run_tracking(
    input_file="radial_avg_curv_vort.nc", save_file="AEW_tracks_raw.nc", initiation_bounds=(-35, 40), radius_used=600, threshold_initial=2e-6, threshold_continue=1e-7, threshold_continue_extrap=1e-6, extrap_day_limit=3, extrap_dist=700, extrap_dist_carib=500, extrap_latitude_max=50,
    extrap_latitude_min=5,
    extrap_longitude_start=-20,
    extrap_latitude_start=20,
    carib_longitude_start=-60,
    AEW_day_remove=2,
    centroid_radius=600,
    spatial_res=1,
    temporal_res=6,
    run_animation=True,
    speed_limit_in=True,
)`
This step runs the AEW tracker on the computed CV output from the previous steps. It is fairly quick to run.

## **Post-Processing of AEW Tracks** `qtrack.tracking.run_postprocessing()`
This step computes the netCDF4 files and saves the data there. It also tries to clean up the tracked AEW data. Importantly, there is a setting that eliminates AEWs that are not at least 2 days long. This can be adjusted within the script if necessary. There is also a feature to identify developing AEWs using HURDAT data (ONLY use this for reanalysis inputs), but this is buggy at the moment. This is why it is required that a year be input into this call, as this allows the code to select the correct year for comparision.

# Details on optional functions

## **Download Example Data** `qtrack.download_examples(input_str, output_dir="")`
Downloads example datasets for AEW tracking. Available files are included above.
- input_str: Input string (from list above) indicating example dateset to be downloaded.
- output_dir (Default: ""): directory path to save downloaded data.
