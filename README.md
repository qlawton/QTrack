# QTrack: An African Easterly Wave Tracker
**Version: 0.0.1**

The AEW tracker was developed by **Quinton Lawton**, currently affiliated with the NSF National Center for Atmospheric Research (NSF NCAR). This python module was developed with the support of **Zachary Moon** (NOAA ARL/Texas A&M University) and **Kelly Núñez Ocasio** (Texas A&M University).  

This is similar to the AEW Tracker described by Lawton et al. (2022), which can be found here: https://doi.org/10.1175/MWR-D-21-0321.1.

Existing AEW tracks for ERA5 data from 1979 - 2022 are already produced and available here: https://osf.io/j4hpq/

## module installation
`pip install qtrack`

## How the tracker works
The AEW tracker uses an input file containing 700hPa zonal (u) and meridional (v) wind and computes a curvature vorticity (CV) term that is smoothed using a 600km average around each gridpoint. Local maxima in this CV are then tracked in two different ways. Over the African continent, meridional averages of CV are used to find longitudinal local maxima, which are used as a first guess for a centroid location. Set thresholds of CV value define when a new AEW is initialized and whether an existing one continues being tracked. Over the Atlantic, the first guess for an existing AEW's location is estimated by the AEW's propagation speed, and it is later found again using a centroid. A post processing script cleans up some of the tracks.

More detailed information on the tracking, not provided in the Lawton et al. (2022) paper or here, can be found in an online technical guide here: https://osf.io/6hqy5

## Citing
If you use this AEW tracker in your research, we ask that you acknowledge your use of this package and provide a citation to the original paper documenting the tracker.

> Lawton, Q. A., S. J. Majumdar, K. Dotterer, C. Thorncroft, and C. J. Schreck, 2022: The Influence of Convectively Coupled Kelvin Waves on African Easterly Waves in a Wave-Following Framework. Monthly Weather Review, 150(8), 2055-2072,  https://doi.org/10.1175/MWR-D-21-0321.1.

## Important Note on Non-Divergent Wind Step
Due to inconsistencies in the windspharm package, which computes the non-divergent wind using spherical harmonics, the non-divergent wind step is currently not included in the package. This step was largely superfluous and thus not including this step is is not anticipated to have any major negative impacts. Nevertheless, it could result in AEW tracks slightly differing from those produced with this step included, including the AEW tracks in the Lawton et al. (2022) AEW databases. We hope to include a non-divergent wind step in a future release.

# Getting Started

## Recommended Workflow Steps
- Adjust input U/V wind data (preferably at or near pressure level 700hPa) to 1x1 degree resolution and 6 hourly temporal resolution. This must be done prior to running this code.
- Prepare input data using `prep_data` function
- Compute curvature vorticity fields using `curvvort` function
- Run the AEW tracking code using `run_tracking` function
- Run post processing of AEWs using `run_postprocessing` function

## Example Workflow
It is highly recommended that users view the [`example_tracking.ipynb`](https://github.com/qlawton/QTrack/blob/main/example_tracking.ipynb) or [`example_tracking.py`](https://github.com/qlawton/QTrack/blob/main/example_tracking.py) files located [on the github page](https://github.ocm/qlawton/QTrack) to see an example of AEW tracking.

## Example Data Files
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

## Input Wind files
**You need u and v winds on a 1x1 degree grid with 6-hourly temporal outputs**. This is the optimal resolution tested in Lawton et al. (2022), but it can be adjusted manually in the code if necessary to change. Furthermore, it is highly recommended that at least 10 days of data are included. The tracking needs a bit of spinup, so if you are running this on model output, it is recommended you append a few days to 1 week of analysis/reanalysis data prior to the first model timestep. Input data should geographically cover at least a portion of the African continent and/or Atlantic Ocean. Note that by default, the tracker does not initiate new waves west of 35W, but this can be adjusted in the tracking function.

A sample model data/reanalysis combination script is provided under `tools/wind_combine.py`

## Output AEW Track files
Two different kinds of output files are produced, a netCDF file and a python pickle file (with AEWs saved as objects). It is recommended that you use netCDF if possible. However, if you need to use the pickle files for whatever reason, you will need to have the `AEW_module.py` script available locally when you open them.

Details on the AEW objects and AEW_module.py are available here: https://osf.io/jnv5u

# Essential Function Details

## **Prep Data** `qtrack.prep_data`

This function ensures the AEW tracker runs properly on the included data. Data should already be at 1x1 degree resolution and 6 hourly temporal resolution. This script will flip the lon/lat coordinates to their property directions, check that all necessary components are included, and slice out the desired height-level if it's 3D data.

*`qtrack.prep_data(data_in, cut_lev_val=700, data_out="prepped_data_for_tracking.nc")`*

- **data_in**: The input data. Should be 2- or 3-dimensional data with U-wind component labeled as "u" and V-wind component labelled as "v". The data file should contain longitudes named "longitude", "lon", or "lons", and latitudes labeled "latitude", "lat", or "lats'.
- **cut_lev_val (Default: 700)** If 3D data, the pressure level to cut out for AEW tracking. Recommended to be near 700hPa.
- **data_out (Default: "prepped_data_for_tracking.nc")**: The name of the output data.

## **Compute Curvature Vorticity and Averaging** `qtrack.curvvort.compute_curvvort`

This step computes the CV and takes the gridpoint averages. It is slow. **However, there is an option within the script to use multiprocessing, which is highly recommended.**

*`qtrack.curvvort.compute_curvvort(data_in, data_out="radial_avg_curv_vort.nc", radius_of_avg=600, data_resolution=1, njobs_in=1, nondiv_wind=False, run_animation=False, gif_dir_in="")`*

- **data_in**: Prepped input U and V data used for computation of CV. Recommended that `prep_data` function is run before computing CV.
- **data_out (Default: "radial_avg_curv_vort.nc')**: Name of CV file to be saved out.
- **radius_of_avg (Default: 600)**: Radius (km) of CV averaging around each grid point.
- **data_resolution (Default: 1)**: Spatial resolution (degrees) of input dataset.
- **njobs_in (Default: 1)**: For multiprocessing (through joblib), the number of CPU jobs to run at a given time. Must be equal to or below the available number of CPUs. -1 will utilize all available CPUs. The more you use, the faster this step will run.
- **nondiv_wind (Default: False)**: DEPRECATED, KEEP AT TRUE.
- **run_animation (Default: False)**: Boolean switch to output an optional animation of CV, saved at the specified gif_dir_in. This will slow down the computation.
- **gif_dir_in (Default: "")**: If animation is output, the file directory to output the animation.

## **AEW Tracking** `qtrack.tracking.run_tracking`

This step runs the AEW tracker on the computed CV output from the previous steps. It is fairly quick to run.

*`qtrack.tracking.run_tracking(input_file="radial_avg_curv_vort.nc", save_file="AEW_tracks_raw.nc", initiation_bounds=(-35, 40), radius_used=600, threshold_initial=2e-6, threshold_continue=1e-7, threshold_continue_extrap=1e-6, extrap_day_limit=3, extrap_dist=700, extrap_dist_carib=500, extrap_latitude_max=50, extrap_latitude_min=5, extrap_longitude_start=-20, extrap_latitude_start=20, carib_longitude_start=-60, AEW_day_remove=2, centroid_radius=600, spatial_res=1, temporal_res=6, run_animation=True, speed_limit_in=True)`*

- **input_file (Default: "radial_avg_curv_vort.nc")**: name of input curvature vorticity file.
- **save_file (Default: "AEW_tracks_raw.nc")**: name of raw AEW output file to be saved.
- **initiation_bounds (Default: (-35, 40)): Longitudes for which the tracker will allow new AEWs to be initiated, from west to east.
- **radius used (Default: 600)**: Averaging radius used in previous CV calculation step.
- **threshold_initial (Default: 2e-6)**: CV threshold (s-1) for initiating new AEW event.
- **threshold_continue (Default: 1e-7)**: CV threshold (s-1) for continued tracking of AEWs over land.
- **threshold_continue_extrap (Default: 1e-6)**: CV threshold (s-1) for continued tracking of AEWs when the extrapolation step has been started. Purposely higher over land than to limit extrapolation errors.
- **extrap_day_limit (Default: 3)**: How long (days) an AEW must exist before the tracker switches to extrapolation from meridional averaging method. AEW must also be west of "extra_longitude_start" for extrapolation to occur.
- **extrap_dist (Default: 700)**: Maximum distance (km) extrapolated center can be from last timestep to be included in new wave track.
- **extrap_dist_carib (Default: 500)**: Maximum distance (km) extrapolated center can be from lat timestep to be included in new wave track, but for AEWs that exist at or west of "Caribbean" longitude threshold (defined in carib_longitude_start).
- **extrap_latitude_max (Default: 50)**: Maximum latitude extrapolation will allow AEW tracks to be extended to.
- **extrap_latitude_min (Default: 5)**: Minimum latitude extrapolation will allow AEWs to be extended to.
- **extrap_longitude_start (Default: -20)**: Longitude where extrapolation can occur, at or west of this point.
- **extrap_latitude_start (Default: 20)**: Extrapolation will automatically be started, no matter the longitude, if the wave center is at or above this latitude.
- **carib_longitude_start (Default: -60)**: Longitude threshold of the "Caribbean". This is the point for which the extrapolation distance rules change to "extrap_dist_carib". Purposely smaller than over open ocean due to the influence of non-AEW CV centers that can create erroneous AEW tracks.
- **AEW_day_remove (Default: 2)**: DEPRECATED, CHANGE IN POST-PROCESSING STEP NOT HERE.
- **centroid_radius (Default: 600)**: Radius (km) of weighted centroid used to find AEW centers.
- **spatial_res (Default: 1)**: Spatial resolution (degrees) of the input data.
- **temporal_res (Default: 6)**: Temporal resolution (hours) of the input data.
- **run_animation (Default: True)**: DEPRECATED, NO ANIMATION IS GENERATED.
- **spped_limit_in (Default: True)**: Boolean, determines if a AEW "speed limit" is incorporated during tracking to prevent erroneously fast, slow, or backwards tracks.

## **Post-Processing of AEW Tracks** `qtrack.tracking.run_postprocessing`

This step computes the netCDF4 files and saves the data there. It also cleans up the tracked AEW data (remove duplicates, combine similar tracks, etc.). Importantly, there is a setting that eliminates AEWs that do not exist for the specified period. This can be adjusted within the script if necessary. There is also a feature to identify developing AEWs using HURDAT data (ONLY use this for reanalysis inputs), but this is buggy at the moment. This HURDAT step is why it is required that a year be input into this call, as this allows the code to select the correct year for TC tracks.

*`qtrack.tracking.run_postprocessing(input_file="AEW_tracks_raw.nc", curv_data_file="radial_avg_curv_vort.nc", radius_used=600, AEW_day_remove=2, real_year_used="None", AEW_merge_dist=500, AEW_forward_connect_dist=700, AEW_backward_connect_dist=200, TC_merge_dist=500, TC_pairing=False, TC_pair_lat_max=25, remove_duplicates=True, hovmoller_save=True, object_data_save=True, netcdf_data_save=True, save_obj_file="AEW_tracks_post_processed.pkl", save_nc_file="AEW_tracks_post_processed.nc", hov_save_file="final_hovmoller.png", hov_name_prefix="", hov_AEW_lat_lim=25,hov_over_africa_color=True)`*

- **input_file (Default: "AEW_tracks_raw.nc")**: path to input raw AEW file to be processed.
- **curv_data_file (Default: "radial_avg_curv_vort.nc")**: path to CV data file from previous step, necessary to compute the AEW strength.
- **radius_used (Default: 600)**: Radius (km) of averaging used in previous CV computation step.
- **AEW_day_remove (Default: 2)**: Number of days a AEW needs to have a track to be included in final dataset. Helps eliminate short tracks that are likely associated with MCSs over Africa or erroneous.
- **real_year_used (Default: "none")**: Year of data, used to compare AEW tracks to TC tracks for the separation of developing and non-developing AEWs. If you are using model data or data without a real-world analogy, put "none".
- **AEW_merge_dist (Default: 500)**: Maximum distance (km) apart two AEW tracks have to be in order to be merged. Used to determine existence of duplicate tracks.
- **AEW_forward_connect_dist (Default: 700)**: Maximum distance (km) in the forward (west) direction AEWs can be reconnected across time steps. This helps fix AEW tracks that may have erroneously gotten separated into two tracks.
- **AEW_backward_connect_dist (Default: 200)**: Maximum distance (km) in the backwards (east) direction AEWs can be reconnected across time steps.
- **TC_merge_dist (Default: 500)**: Distance (km) threshold used to associated AEWs with TCs. If a AEW is at or closer to this distance from a TC's genesis point (HURDAT), the AEW will be considered a developing wave and relevant information will be saved in the output file(s).
- **TC_pairing (Default: False)**: Boolean indicating whether AEWs should be paired with TC data. Only use if your input data corresponds to observatons/reanalysis. Will not work if year is not specified in "real_year_used".
- **TC_pair_lat_max (Default: 25)**: For optional hovmoller output, the maximum latitude of TC tracks to be displayed on hovmoller diagram.
- **remove_duplicates (Default: True)**: Boolean indicating whether the code should try to remove duplicate AEW tracks.
- **hovmoller_save (Default: True)**: Boolean indicating whether an output hovmoller should be generated after post-processing.
- **object_data_save (Default: True)**: Boolean indicating whether AEW object files should be saved after post-processing.
- **netcdf_data_save (Default: True)**: Boolean indicating whether netCDF files should be saved after post-processing.
- **save_obj_file (Default: "AEW_tracks_post_processed.pkl")**: Name of output object file, if saved.
- **save_nc_file (Default: "AEW_tracks_post_processed.nc")**: Name of output netCDF file, if saved.
- **hov_save_file (Default: "final_hovmoller.png")**: Name of output hovmoller, if saved.
- **hov_name_prefix (Default: "")**: Name to be displayed in title of hovmoller plot.
- **hov_AEW_lat_lim (Default: 25)**: For optional hovmoller output, the maximum latitude of AEW tracks to be displayed on hovmoller diagram.
- **hov_over_africa_color (Default: True)**: If true and hovmoller is produced, use a different color to indicate AEWs tracks that originate over Africa versus those that originate over the ocean.  

# Optional Function Details

## **Download Example Data** `qtrack.download_examples`

Downloads example datasets for AEW tracking. Available files are included above.

*`qtrack.download_examples(input_str, output_dir="")`*

- **input_str**: Input string (from list above) indicating example dateset to be downloaded.
- **output_dir (Default: "")**: directory path to save downloaded data.
