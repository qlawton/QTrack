# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 16:06:58 2020

@author: Quin7
"""
import numpy as np
from COMPUTE_SAVE_CURV_VORT_NON_DIV_UPDATE_FIX_PARALLEL import COMPUTE_CURV_VORT_NON_DIV_UPDATE
import sys

if len(sys.argv) == 1: #Basically, no inputs
    data_in = 'CURV_VORT/HELMHOLTZ/wind_700_helmholtz.nc'
    data_out = 'CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc'
elif len(sys.argv) == 2: #Only input the data-in
    data_in = sys.argv[1]
    data_out = 'CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc'
elif len(sys.argv) == 3:
    data_in = sys.argv[1]
    data_out = sys.argv[2]
else:
    raise Exception("Too many arguments provided.")

# SETTINGS #
SAVE_OUTPUT = True
RUN_ANIMATION = False
res = 1 #Resolution of input data (Default: 1 [1x1 data]. Not recommended one deviates from this)
rad = 600 #Radius of averaging (km) used here. (Default: 600).
nondiv = True #Non-divergent component of wind (Default: True. Only true if you have global data and can run the 'non_divergent_wind.py' script from before.

###RUNNING CURV VORT IN PARALLEL USING JOBLIB
njobs = 1 #Set to 1 for non-parallel, -1 for all available CPUs, or to the number of CPUs requested. (Default = 1)

COMPUTE_CURV_VORT_NON_DIV_UPDATE(data_in, data_out, res, rad, njobs, nondiv, RUN_ANIMATION, SAVE_OUTPUT)
