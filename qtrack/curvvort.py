def compute_curvvort(data_in, data_out="radial_avg_curv_vort.nc", radius_of_avg=600, data_resolution=1, njobs_in=1, nondiv_wind=False, run_animation=False, gif_dir_in=""):
    """
    Compute the Curvature Vorticity averaged within a radius of each gridpoint. This is computationally expensive, and can be run in parallel.


    """
    # from COMPUTE_SAVE_CURV_VORT_NON_DIV_UPDATE_FIX_PARALLEL import COMPUTE_CURV_VORT_NON_DIV_UPDATE

    from qtrack.core import COMPUTE_CURV_VORT_NON_DIV_UPDATE

    #     if len(sys.argv) == 1: #Basically, no inputs
    #         data_in = 'CURV_VORT/HELMHOLTZ/wind_700_helmholtz.nc'
    #         data_out = 'CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc'
    #     elif len(sys.argv) == 2: #Only input the data-in
    #         data_in = sys.argv[1]
    #         data_out = 'CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc'
    #     elif len(sys.argv) == 3:
    #         data_in = sys.argv[1]
    #         data_out = sys.argv[2]
    #     else:
    #         raise Exception("Too many arguments provided.")

    # SETTINGS #
    SAVE_OUTPUT = True
    RUN_ANIMATION = run_animation
    res = data_resolution  # Resolution of input data (Default: 1 [1x1 data]. Not recommended one deviates from this)
    rad = radius_of_avg  # Radius of averaging (km) used here. (Default: 600).
    nondiv = nondiv_wind  # Non-divergent component of wind (Default: True. Only true if you have global data and can run the 'non_divergent_wind.py' script from before.

    ###RUNNING CURV VORT IN PARALLEL USING JOBLIB
    njobs = njobs_in  # Set to 1 for non-parallel, -1 for all available CPUs, or to the number of CPUs requested. (Default = 1)

    COMPUTE_CURV_VORT_NON_DIV_UPDATE(data_in, data_out, res, rad, njobs, nondiv, RUN_ANIMATION, SAVE_OUTPUT)
