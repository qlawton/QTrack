# -*- coding: utf-8 -*-
"""
Main African Easterly wave Tracking script. Requires an input of radially averaged curvature vorticity at the 700 hPa level. A description of this methodology can be found in the associated Monthly Weather Review publication: https://doi.org/10.1175/MWR-D-21-0321.1

More technical details on the algorithm can be found here: https://osf.io/6hqy5

This tracking algoritm was written by Quinton Lawton while at the University of Miami.
"""

def run_tracking():
    ##### ALL SETTINGS #####
    # -----------------------------------------------------------------------------
    ##### INPUT SETTINGS  #####
    #The following variables control the data being ingested, the name of the system
    #(and corresponding datafiles) being focused on, and information on the radius
    #of curvature vorticity averaging we want to pull from.
    import numpy as np
    import sys
    np.warnings.filterwarnings('ignore')

    #### INFO OF DATA TO RUN ON #####

    rad_used = 600 #km, (Default: 600km)
    out_attr = 'preliminary_tracks' #NAME OF RUN
    extra_info = out_attr
    temp_res = 6 #IMPORTANT -- set temporal resolution of data (hour)
    res = 1 #IMPORTANT -- set spatial resolution of data (deg)

    #### BINARY SETTINGS/OPTIONS -- Decide whether to run certain parts of code
    run_extender = False #Turn linear extending of smooth tracks at begin/end on and off (Default: False)
    raw_ani = True #True will plot "raw" centers, false will plot smoothed. (Default: True)
    wide_ani = True #Include up to gulf with true, otherwise stick over Atlantic with false (Default: True)
    wind_overlay = True #IF ANIMATION CREATED, decide if wind vectors overlain or not (Default: True)
    duplicate_removal = False #Determines if we remove duplciate tracks here (Default: True)
    save_data = True #Turn on and off the saving of output data (netCDF4) (Default: True)
    merge_wave = False #Merge wave during run (Default: False) (NOTE: we end up doing this in postprocessing now)
    speed_limit = True #Turn on or off forward speed limit, prevening backward motion/stalling (Default: True)

    # -----------------------------------------------------
    #### BANDING (ELLESS AND TORN INSPIRED) SETTINGS ####
    # -----------------------------------------------------
    # -- Over land, the tracker uses a similar method to that of Elless and Torn (2018) to add wave points. This
    #is also used as a backup method over the ocean (ONLY EAST OF TRANSITION LONGITUDE) if extrapolation does not
    #extend the wave, typically due to meridional shifts in wave position. 


    banding_t = True # If True, will look at 6 bands (5-15, 6-16,...10-20) (default: True)

    #THRESHOLD SETTINGS (for Curvature Vorticity CV)
    thres_init = 2e-6 #threshold to initialize a AEW event
    thres_cont = 1e-7 #threshold for continued tracking of AEWs

    # TURN FEATURES ON AND OFF FOR THIS PETHOD (default: True)
    force_forward = True #Allow "backward" (eastward) points to be considered, but when appended longitude will be equal to previous AEW track point
    force_bump = True #Preference "forward" (westward) points to be appended to existing AEW tracks over backwards/nearby

    # DISTANCE TESTING
    #NOTE -- IGNORE 3/6 hr denotations. In reality, these refer to 1 or 2 temporal time steps backwards (for 6hr data, that would be 6hr and 12hr not 3hr and 6hr)
    step_3hr = 700 #How close new AEW point must be to existing AEW, looking at 1 timestep before (KM) (default: 700km)
    step_6hr = 1000 #Same but looking at track data for 2 timesteps before (KM) (default: 1000km)
    bump_num = 500./step_3hr #(NOTE: NOT RELEVANT IN FINAL VERSION OF CODE) When "bumping" points (if turned on), only do so for points within this RATIO of the forward step (i.e. 5/8*800 = 500KM) 
    back_cutoff_land = 100 #Only consider "backward" points within this radius of existing AEW (KM), here over land (default: 100km)
    back_cutoff_ocean = 300 #Same setting as before, but now over teh ocean (default: 300km)
    back_cutoff_long_land = 100 #"Backward" point setting but for extending waves with a missing timestep, land (default: 100km)
    back_cutoff_long_ocean = 300 #"Backward" point setting but for extending waves with a missing timestep, ocean (default: 300km)
    long_edit = True #Allow tracker to use the banding method to try to extend waves with a missing timestep (default: True)
    land_lat_limit = 3 #RESTRICTION over land, this is the degrees latitude a wave can shift meridionally for each timestep (default: 2)
    #Note that for the missing timestep portion, the latitude limit becomes 3*land_lat_limit

    #Merge settings
    merge_distance = 500 #Merge distance, obsolete now as it's handed in the post processing #(Default: 500km)

    # WAVE INITATION SETTINGS
    lon_west = -35# degrees -- AEWs can only be initiated east of this longitude
    lon_east = 40 #degrees -- AEWs can only be initated west of this longitude
    step_3hr_init = step_3hr #To initiate new point, it must not be within this radius of an existing AEW 1 timestep ago (KM)
    step_6hr_init = step_6hr #Same as before but for two timesteps ago (KM)
    stuck_thresh_day =1 #Number of days can go by with wave stuck in same location before cutting track
    stuck_lon = 2

    # -----------------------------------------------------
    #### EXTRAPOLATION SETTINGS ####
    # -----------------------------------------------------
    # -- Over the ocean, this code will first attempt to extend a waves track by predicting it's next
    # position and then running a mass-centering script with said position as it's first guess. This
    # will only be added to the existing wave if it's within a certain distance of the last wave timestep,
    # and not too weak. This is to prevent runaway track extension and to try to keep the tracker from following
    # unrelated blobs of curvature vorticity.

    ##### OUTPUT 
    EXTRAPOLATE = True #Turn on or off the extrapolation entirely (Default, True)
    EXTRAPOLATE_SKIP = True #Allow extrapolation code to run again 2-3 timesteps back in case wave was only temporarily lost (Default: True)
    EXTRAPOLATE_FORWARD = True #Allow thresholds to prevent too much extrapolation in wrong direction, with parameters set below (Default: True)
    extrapolate_thresh = 1e-6 #Strength final centroid has to be for centroid to be appended to wave (default: 1e-6)
    speed_curv_thresh = 8e-6  #EXCEPTION, allowing free movement of wave centers at or above this threshold w/o restriction (Default: 8e-6)
    extrapolate_time = 24/temp_res*3 #How long a wave must exist before extrapolation can occur (Default, 24/temp_res*3 or 3 days)
    extra_distance = 700 #Distance (km) extrapolated center can be from last timestep, otherwise excluded (default: 700km)
    extra_dist2 = 700 #If running on skipped point, distance (km) extrapolated center can be from timestep 2 back or otherwise excluded (default: 1000km)
    extra_dist3 = 700 #If running on skipped point, distance (km) extrapolated center can be from timestep 3 back or otherwise excluded (default: 1000km)

    extra_distance_60 = 500 #Distance (km) extrapolated center can be from last timestep, at or west of transition longitude (default, 500km)
    extra_dist2_60 = 700 #Same as before, but west of transition longitude (default: 1000km)
    extra_dist3_60 = 700 #Same as before, but west of transition longitude (default: 1000km)

    lat_max = 50 #The maximum latitude tracker will run extrapolation to (default: 50)
    extra_lat_bottom = 5 #The minimum latitude tracker will run extrapolation on, otherwise defaults to other method (default: 5)
    extra_lon = -20 #The longitude at and west of that extrapolation begin (default: -20)
    extra_lat_start = 20 #Was at 18, perhaps now do 20 to account for weird action just N of lesser antilles (default: 20)
    extra_transition = -60 #Longitude at and west of that extrapolation settings become more restrictive at (default: -60)
    extra_it_20 = 5 #Number of timesteps used for linear extrapolation east of transition threshold (default: 5)
    extra_it_60 = 5 #Number of timesteps used for linear extrapolation west of transition threshold (default: 5)
    extra_backcut = 100 #Maximum distance (km) "backwards" (to east) extrapolation can put a new wave point, unless there is exceptions (default: 100)
    #Right now, 12 hour land_lat_limit is three times that of the 6_hour limit above

    # -----------------------------------------------------
    #### MASS CENTROID SETTINGS ####
    # -----------------------------------------------------
    # A curvature voticity centroid is run on local maxima before distance testing and either AEW initation or adding of points.
    centroid_rad = 600 #Radius to weight centroid 
    centroid_rad_extra = 600 #Radius to weight extrapolation centroid
    centroid_rad_extra_60 = 600 #Same as above, but west of transition longitude
    centroid_it = 0 #Number of iterations to run centroid code for banding
    centroid_it_extra = 0 #Same as above, but for the extrapolation code

    # -----------------------------------------------------
    #### CLEANUP SETTINGS ####
    # -----------------------------------------------------

    cleanup = True #Turn cleanup of points on/off 
    upper_limit = False #old, obsolete feature. (default: False)
    cut_stuck = False# old, obsolute feature. (default: False)
    days_remove = 2 #Number of days a AEW is required to last to be put in data
    deg_sep = 2 #Degrees of separations initiation points must be from each other.
    arb_thresh = 0.5e-6 
    smooth_len = 7


    # -----------------------------------------------------------------------------
    ##### END SETTINGS #####

    final_merge_list = np.array([])

    stuck_thresh = int(stuck_thresh_day*24/temp_res)

    ##### IMPORT STATEMENTS AND FUNCTIONS #####
    # -----------------------------------------------------------------------------

    # IMPORT STATEMENTS
    import numpy as np
    import datetime
    from netCDF4 import Dataset, num2date
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import matplotlib
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import metpy.calc as mpcalc
    import time as tm
    from scipy import signal
    from geopy.distance import geodesic
    from sklearn.linear_model import LinearRegression
    from scipy.signal import savgol_filter

    # USER-DEFINED FUNCTIONS
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        #print(idx)
        return idx, array[idx]


    def haversine(lon1, lat1, lon2, lat2):
        from math import radians, cos, sin, asin, sqrt
        """
        Calculate the great circle distance between two points 
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians 
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 
        # Radius of earth in kilometers is 6371
        km = 6371* c
        return km

    def rad_mask(i, j, dx, dy, radius):
        '''This is a computationally efficient (at least, in python-world) way of calculating a mask of gridpoints around a
        center point. This generates an array of points with values corresponding to the distance from some center point. 
        Then the code masks out any values above a certain radius. 

        This uses the assumption of a spherical Earth, not accounting for the equatorial "bulge" in real life.'''
        start = tm.time()

        earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
        lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
        lat_km = lat_met/1000
        res = 1
        buffer = int(np.ceil(radius/lat_km/res))
        boolean_array = np.zeros(np.shape(dx), dtype=bool)
        #i_array = np.zeros(np.shape(dy))
        #j_array = np.zeros(np.shape(dx))
        dy = -dy

        # Before doing this, we recognize something -- there is a max radius in the x and y direction that will
        # be computed. To save on computational time, we need to slice out the extraneous points that we already
        # know from this NOT to be in the circle. How? We know the absolute highest distance between gridspaces
        # will be the maximum value of latitude on the elliptical earth... ~ 111 km. So slice out roughly a 112km "box" 
        # plus one gridbox for buffer. 
        i_st = (i-(buffer+1))
        i_end = (i+(buffer+1))
        j_st = (j-(buffer+1))
        j_end = (j+(buffer+1))

        new_i = i-i_st
        new_j = j-j_st

        dy_slc = dy[i_st:i_end,j_st:j_end]
        dx_slc = dx[i_st:i_end,j_st:j_end]

        i_array_sub = np.zeros(np.shape(dy_slc))
        j_array_sub = np.zeros(np.shape(dx_slc))

        i_array_sub[new_i,:] = 0
        i_array_sub[(new_i+1):,:] = np.add.accumulate(dy_slc[(new_i+1):,:])
        i_array_sub[(new_i-1)::-1,:] = np.add.accumulate(dy_slc[(new_i-1)::-1,:])

        j_array_sub[:,new_j] = 0
        j_array_sub[:,(new_j+1):] = np.add.accumulate(dx_slc[:,(new_j+1):], axis=1)
        j_array_sub[:,(new_j-1)::-1] = np.add.accumulate(dx_slc[:,(new_j-1)::-1], axis=1)

        radial_array = (np.sqrt(np.square(i_array_sub)+np.square(j_array_sub))/1000) < radius # in km
        boolean_array[i_st:i_end,j_st:j_end] = radial_array
        end = tm.time()
        #print(end - start)
        return boolean_array

    def general_centroid(x_in, y_in, PV_in, radius, lati, loni, it, exclude=True):
        '''Similar to that used for calculating PV centers, this computes a centroid point for a given field, using a prescribed
        radius as specified in the function input. Can iterate a given number of times: set in=0 for no iterations.

        Uses the more computationally efficient (rad_mask) function to test if a point is within a radial distance or not.'''
        #start= time.time()

        # What does our lat/lon refer to?
        guess_lat = y_in[lati]
        guess_lon = x_in[loni]
        #Other important stuff
        res = 1
        #box_num = 6/res; #First number is the degrees you want to cut out from the dataset for analysis

        def get_dist_meters(lon, lat):
            earth_circ = 6371*2*np.pi*1000 #earth's circumference in meters
            lat_met = earth_circ/360 #get the number of meters in a degree latitude (ignoring "bulge")
            lat_dist = np.gradient(lat, axis=0)*lat_met
            lon_dist = np.gradient(lon, axis=1)*np.cos(np.deg2rad(lat))*lat_met
            return lon_dist, lat_dist 


        ## FIRST ITERATION ##
        #use the last bt_lat, bt_lon as the "genesis point" and commence slicing to make data sizes more managable
        lati, latval = find_nearest(y_in, guess_lat)
        loni, lonval = find_nearest(x_in, guess_lon)

        #Commence slicing out the correct data
        #lat_slice = slice(int(lati-box_num),int(lati+box_num))
        #lon_slice = slice(int(loni-box_num),int(loni+box_num))
        lon_new = x_in#[lon_slice]
        lat_new = y_in#[lat_slice]
        PV_slice = PV_in#[lat_slice,lon_slice]
        #Get the distance array for later use    

        #Make a meshed grid
        LON_X, LAT_Y = np.meshgrid(lon_new, lat_new)
        dx, dy = get_dist_meters(LON_X, LAT_Y)

        #Finally, output a 1/0 filter for PV within the given radius
        #pv_filt = inner_filter(lat_new, lon_new, guess_lat, guess_lon, radius)
        pv_filt = rad_mask(lati, loni, dx, dy, radius)

        PV_mask = PV_slice.copy()
        lon_mask = LON_X.copy()
        lat_mask = LAT_Y.copy()
        PV_mask[pv_filt == 0] = np.NaN
        lon_mask[pv_filt == 0] = np.NaN
        PV_mask[pv_filt== 0] = np.NaN

            #Mask out (for all three fields) the filtered PV locations
        #PV_mask = np.ma.masked_where(pv_filt==0, PV_slice).copy()
        #lon_mask = np.ma.masked_where(pv_filt==0, LON_X).copy()
        #lat_mask = np.ma.masked_where(pv_filt==0, LAT_Y).copy()

        if exclude == True:
            PV_mask[PV_mask<0] = np.NaN
            lon_mask[PV_mask<0] = np.NaN
            lat_mask[PV_mask<0] = np.NaN


        #Calculate the center of mass, finally
        x_cent = np.nansum(lon_mask*PV_mask)/np.nansum(PV_mask)
        y_cent = np.nansum(lat_mask*PV_mask)/np.nansum(PV_mask)

        ## LOOPING THE ITERATIONS ##
        #Now we will iterate over this it-1 times to converge on the true weighted mass. Each time, we will filter based on 
        #the previous answer's "guess/calculation" of where the PV weighted center is.

        if it == 0: #This means we are running to convergence
            dist_check = 100
            loop_it = 0
            while dist_check >= 50: #While the new point is at least 50km away from old
                if loop_it >= 10:
                    break
                old_x = x_cent
                old_y = y_cent
                x_centi, n = find_nearest(lon_new, x_cent)      
                y_centi, n = find_nearest(lat_new, y_cent)  

                try:
                    pv_filt = rad_mask(y_centi, x_centi, dx, dy, radius)
                except:
                    #print('Exception: Edge of boundary. Skipping to next iteration.')
                    break
                #Mask out (for all three fields) the filtered PV locations
                PV_mask = PV_slice.copy()
                lon_mask = LON_X.copy()
                lat_mask = LAT_Y.copy()
                PV_mask[pv_filt == 0] = np.NaN
                lon_mask[pv_filt == 0] = np.NaN
                PV_mask[pv_filt== 0] = np.NaN

                if exclude == True:
                    PV_mask[PV_mask<0] = np.NaN
                    lon_mask[PV_mask<0] = np.NaN
                    lat_mask[PV_mask<0] = np.NaN

                x_cent = np.nansum(lon_mask*PV_mask)/np.nansum(PV_mask)
                y_cent = np.nansum(lat_mask*PV_mask)/np.nansum(PV_mask)

                dist_check = haversine(x_cent, y_cent, old_x, old_y)
                loop_it = loop_it + 1
        else:
            for step in np.arange(it-1):
                #print(x_cent, y_cent)
                x_centi, n = find_nearest(lon_new, x_cent)      
                y_centi, n = find_nearest(lat_new, y_cent)
                try:
                    pv_filt = rad_mask(y_centi, x_centi, dx, dy, radius)
                except:
                    #print('Exception: Edge of boundary. Skipping to next iteration.')
                    continue
                #Mask out (for all three fields) the filtered PV locations
                PV_mask = PV_slice.copy()
                lon_mask = LON_X.copy()
                lat_mask = LAT_Y.copy()
                PV_mask[pv_filt == 0] = np.NaN
                lon_mask[pv_filt == 0] = np.NaN
                PV_mask[pv_filt== 0] = np.NaN

                if exclude == True:
                    PV_mask[PV_mask<0] = np.NaN
                    lon_mask[PV_mask<0] = np.NaN
                    lat_mask[PV_mask<0] = np.NaN

                x_cent = np.nansum(lon_mask*PV_mask)/np.nansum(PV_mask)
                y_cent = np.nansum(lat_mask*PV_mask)/np.nansum(PV_mask)


        #end = time.time()
        #print(str(end-start)+' secs')
        return lon_new, lat_new, pv_filt, PV_mask, x_cent, y_cent

    def find_maxima(data_in, lon, lat, lon_west, lon_east, thres_init, thres_cont, separate_bands = banding_t, exclude = True, smooth = True):
        '''This finds local maxima for a 5-20N AVERAGE VALUE AT EACH LONGITUDE POINT. Two thresholds are considered:
        -- "Init" is the threshold for initially defining a AEW, and here is ONLY defined if east of the "lon_land" value
        -- "Cont" is the threshold for continuing the tracking/propagation of an existing AEW, and is defined everywhere.'''
        avg_list = []    

        def find_nearest(array, value): #This way, we can set the lat range from 5 to 20N
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            #print(idx)
            return idx, array[idx]
        band_len = 6
    #     if extra_bands == True:
    #         band_len = 8
    #     else:
    #         band_len = 6

        ## ----- TAKE ELLESS'S FEEDBACK INTO CONSIDERATION, CONSIDER 5 SEPARATE "BANDS" -----
        if separate_bands == True:
            lat_extra_end, n = find_nearest(lat, 15)
            lat_extra_st, n = find_nearest(lat,5)
            #print(lat_extra_st, lat_extra_end)
            if lat_extra_end > lat_extra_st:
                band_vals = np.zeros((band_len,np.shape(data_in[lat_extra_st:lat_extra_end, :])[1]))
                test_vals = np.zeros((band_len,np.shape(data_in[lat_extra_st:lat_extra_end, :])[1]))
            else:
                band_vals = np.zeros((band_len,np.shape(data_in[lat_extra_end:lat_extra_st, :])[1]))
                test_vals = np.zeros((band_len,np.shape(data_in[lat_extra_end:lat_extra_st, :])[1]))    
            for i in range(band_len): #We will have six bands
                lat_avg_end, n = find_nearest(lat, 15+i) #We only want to average the cells in a 5-20N range
                lat_avg_st, n = find_nearest(lat, 5+i)

                if lat_avg_end > lat_avg_st:
                    data_slice = data_in[lat_avg_st:lat_avg_end, :]
                else:
                    data_slice = data_in[lat_avg_end:lat_avg_st, :]
                #print(lat_avg_st, la)
                if exclude == True: #Are we excluding negative values or not?
                    banded_noneg = data_slice.copy()
                    banded_noneg[banded_noneg<0] = np.NaN
                else:
                    banded_noneg = data_slice.copy()

                data_mean_int = np.nanmean(banded_noneg, axis = 0)
                data_mean_both = np.nanmean(data_slice, axis = 0)
                band_vals[i,:] = data_mean_int[:]
                test_vals[i,:] = data_mean_both[:]

            data_mean_noneg = np.nanmax(band_vals, axis = 0)
            data_mean = np.nanmax(test_vals, axis = 0)
            lat_center = np.argmax(band_vals, axis = 0)+10

            #Need to get final list of latitudes that is the same length as the init/cont out. Actually right now only
            #init since we only initate new waves with the latitude banding... use previous latitude for existing wavese
            #print(lat_center)


        else:
            lat_avg_st, n = find_nearest(lat, 15) #We only want to average the cells in a 5-20N range
            lat_avg_end, n = find_nearest(lat, 5)

            data_slice = data_in[lat_avg_st:lat_avg_end, :]
            data_mean = np.nanmean(data_slice, axis = 0)
            if exclude == True:
                data_slice_noneg = data_slice.copy()
                data_slice_noneg[data_slice_noneg<0] = np.NaN
                data_mean_noneg = np.nanmean(data_slice_noneg, axis = 0)
            else:
                data_mean_noneg = data_mean

        if smooth == True:
            data_mean_noneg = savgol_filter(data_mean_noneg, 3, 2)
            data_mean = savgol_filter(data_mean, 3, 2)

        init_i = signal.argrelextrema(data_mean, np.greater)[0]
        cont_i = signal.argrelextrema(data_mean, np.greater)[0]

        init_max = data_mean_noneg[init_i].astype(float)
        cont_max = data_mean_noneg[cont_i].astype(float)

        if (np.isnan(init_max)).any():
            init_max[np.isnan(init_max)] = 0
        if (np.isnan(cont_max)).any():
            cont_max[np.isnan(cont_max)] = 0

        #Filter out maximas such that they have to both be positive and greater than a given threshold
        init_max = init_max[lon[init_i]<40]
        init_i = init_i[lon[init_i] <40]
        init_max = init_max[lon[init_i]>-100]
        init_i = init_i[lon[init_i]>-100]
        init_i = init_i[init_max >= thres_init] 

        cont_max = cont_max[lon[cont_i] < 40]
        cont_i = cont_i[lon[cont_i] < 40]
        cont_max = cont_max[lon[cont_i] > -100]
        cont_i = cont_i[lon[cont_i] > -100]
        cont_i = cont_i[cont_max >= thres_cont] 

        #Now filter "initial" AEWs such that they have to have a longitude greater than 17W (the Africa land cutoff)
        lon_init = lon[init_i]
        init_i = init_i[lon_init >= lon_west]
        lon_init = lon_init[lon_init >= lon_west]
        init_i = init_i[lon_init <= lon_east]
        #lat_init = lat_center[init_i]

        return data_mean, init_i, cont_i, lat_center

    def remove_nearby(data_list, data_values, arb_thresh, lon_thresh, data_res):
        '''Using an extrema threshold is likely to result in duplicate data points nearby each other. In this
        function, we will test if subsequent maxima points are within a certain distance of one another. Everytime they are not,
        the "counter" will reset and a new list will be generated. Once there is a list of points and the chain is broken, two 
        things will happen:
            - Will find the maximum value, and any points arb_thresh or greater less than it is automatically removed
            - For all remaining points, a centroid is calculated and a final point is output'''
        counter = []
        counter_data = []
        rm_list = []
        adj_point = []
        grid_thresh = lon_thresh/data_res
        #print(grid_thresh)

        for i in range(len(data_list)):
            if (i == 0): #If first element 
                continue #Skip this loop
            if (data_list[i] - data_list[i-1]) <= grid_thresh: #If closer together than supposed to be
                counter.append(data_list[i-1]) #Append previous point
                counter.append(data_list[i]) #Append current point

                #Do same for counter_data points
                counter_data.append(data_values[i-1])
                counter_data.append(data_values[i])

                rm_list.append(data_list[i-1])
                rm_list.append(data_list[i])

            elif len(counter) ==0: #If counter doesn't have anything, go ahead and continue to next loop
                continue

            elif len(counter) !=0: #But if counter does have something, then do our removal process
                count_max = np.max(counter_data)
                cond = (count_max - counter_data)
                counter = np.array(counter)
                counter_data = np.array(counter_data)

                counter = counter[cond < arb_thresh] #Only keep things wihtin arb_thresh of the max
                counter_data = counter_data[cond < arb_thresh] #Do same for actual values

                #take weighted mean of the whole thing

                weighted_point = np.nansum(counter*counter_data)/np.nansum(counter_data)
                weighted_i = int(np.round(weighted_point))
                adj_point.append(weighted_i)
                counter = []
                counter_data = []
            if (i == (len(data_list)-1)): #If it's the last interation
                if len(counter) ==0: #If counter doesn't have anything, go ahead and continue to next loop
                    continue
                elif len(counter) !=0: #But if counter does have something, then do our removal process
                    count_max = np.max(counter_data)
                    cond = (count_max - counter_data)
                    counter = np.array(counter)
                    counter_data = np.array(counter_data)

                    counter = counter[cond < arb_thresh] #Only keep things wihtin arb_thresh of the max
                    counter_data = counter_data[cond < arb_thresh] #Do same for actual values

                    #take weighted mean of the whole thing

                    weighted_point = np.nansum(counter*counter_data)/np.nansum(counter_data)
                    weighted_i = int(np.round(weighted_point))
                    adj_point.append(weighted_i)
                    counter = []
                    counter_data = []



        return adj_point, rm_list

    def execute_cleanup1(cont_max, init_max, data_mean, arb_thresh, lon_thresh, data_res):
        '''This just executes the remove_nearby script and cleans up output to generate the initial list of lon points.'''
        #adj_cont, rm_cont = remove_nearby(cont_max, data_mean[cont_max], arb_thresh, lon_thresh, data_res)
        adj_init, rm_init = remove_nearby(init_max, data_mean[init_max], arb_thresh, lon_thresh, data_res)
        final_max_init = []
        final_max_cont = []

        # ---------- DO IT FOR CONT --------------
    #     for i in range(len(cont_max)):
    #         if cont_max[i] not in rm_cont:
    #             final_max_cont.append(cont_max[i])

    #     for i in range(len(adj_cont)):
    #         final_max_cont.append(adj_cont[i])
    #     final_max_cont.sort()

        # ---------- DO IT FOR INIT ---------------
        for i in range(len(init_max)):
            if init_max[i] not in rm_init:
                final_max_init.append(init_max[i])

        for i in range(len(adj_init)):
            final_max_cont.append(adj_init[i])
        final_max_cont.sort()

        #REMOVE DUPLICATES (Need to check... why does this happen??)
    #     new_cont = []
        new_init = []

    #     [new_cont.append(x) for x in final_max_cont if x not in new_cont] 
        [new_init.append(x) for x in final_max_init if x not in new_init] 

    #     final_max_cont = new_cont
        final_max_init = new_init
        return False, final_max_init

    def execute_cleanup2(final_max_cont, final_max_init, data_mean, wide_thresh, data_res):
        '''This just executes the wide_wave script and cleans up output to generate the final list of lon points.'''
        wide_array, rm_wide, adj_wide = wide_wave_tester(final_max_cont, data_mean, wide_thresh, data_res)
        wide_array_init, rm_wide_init, adj_wide_init = wide_wave_tester(final_max_init, data_mean, wide_thresh, data_res)
        #print(adj_wide)
        new_final_cont = []
        new_final_init = []
        # ---------- DO IT FOR CONT --------------
        for i in range(len(final_max_cont)):
            if final_max_cont[i] not in rm_wide:
                new_final_cont.append(final_max_cont[i])

        for i in range(len(adj_wide)):
            new_final_cont.append(adj_wide[i])
        new_final_cont.sort()

        # ---------- DO IT FOR INIT ---------------
        for i in range(len(final_max_init)):
            if final_max_init[i] not in rm_wide_init:
                new_final_init.append(final_max_init[i])

        for i in range(len(adj_wide_init)):
            new_final_init.append(adj_wide_init[i])
        new_final_init.sort()


        #print(lon_cont)
        #print(lon_init)

        #print(wide_array, rm_wide, adj_wide)
        lon_init = lon[new_final_init]
        lon_cont = lon[new_final_cont]
        return lon_init, lon_cont, new_final_cont, new_final_init

    def within_distance_direct(lon, lat, lon_in,lat_in, lon_old, lat_old, step, forward_weight = True, forward_scale = 1/5):
        '''Right now, we want centers to be within a distance of the time step. 

        Forward bias prevents extreme cases of "back_building" by limiting back_lon search to half the forward distance.'''
        if lon_in <= lon_old: #IF true, then the wave is moving the correct direction
            fwd = True
        else:
            fwd = False
        distance = haversine(lon_in, lat_in, lon_old, lat_old)
        #print(distance)
        if distance <= step:
            return True, distance, fwd
        else: 
            return False, distance, fwd
    # lon1 = lon[25]
    # lon2 = lon[25]
    # lat1 = lat[25]
    # #lat2 = lat1
    # lat2 = lat[30]

    #print(lon1, lon2)
    #haversine(lon1, lat1, lon2, lat2)
    #print(lon1, lon2, lat1, lat2)
    #print(within_distance_direct(lon, lat, lon1, lat1, lon2, lat2, 300))
    #win_3, win_6 = within_distance(-100, -98, 300, 600, 1)
    #print(win_3, win_6)

    def nan_helper(y):
        """Helper to handle indices and logical indices of NaNs. CREDIT TO: user "eat" on StackOverflow

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]

    def cleanup_AEW(AEW_lon_in, AEW_lat_in, temporal_res, days_remove):
        del_list = []
        '''Oftentimes, a wave is "lost" to the tracker for 3-9 hours due to not meeting a threshold requirement, but then
        re-acheiving this requirement later on. This leaves annoying "nan" values in the missing places. Thus, this script
        will fill these holes with linear interpolated values.'''
        cutoff_len = int(np.round(days_remove*24/temporal_res))
        for row in range(np.shape(AEW_lon_in)[0]):
            AEW_lon_slc = AEW_lon_in[row,:]
            AEW_lat_slc = AEW_lat_in[row,:]
            # First, find the indices of the first 'non-nan' value and the last real value
            real_val = np.argwhere(~np.isnan(AEW_lon_slc))
            first = real_val[0][0]
            last = real_val[-1][0]

            #Linear interpolation of "intermediate" points
            nans, x = nan_helper(AEW_lon_slc)
            AEW_lon_slc[nans]= np.interp(x(nans), x(~nans), AEW_lon_slc[~nans])

            nans, x = nan_helper(AEW_lat_slc)
            AEW_lat_slc[nans]= np.interp(x(nans), x(~nans), AEW_lat_slc[~nans])

            AEW_lon_in[row,:] = AEW_lon_slc
            AEW_lat_in[row,:] = AEW_lat_slc

            #But this resulted in interpolations before and after the data -- thus, need to use "first, last" to restore to Nans
            AEW_lon_in[row, 0:first] = np.NaN
            AEW_lon_in[row, (last+1):] = np.NaN
            AEW_lat_in[row, 0:first] = np.NaN
            AEW_lat_in[row, (last+1):] = np.NaN

            #Remove any track that has less than 48 hours worth of data
            data_len = len(AEW_lon_slc[~np.isnan(AEW_lon_slc)])
            #print(data_len)
            if data_len <= cutoff_len:
                del_list.append(row)
        AEW_lon_in = np.delete(AEW_lon_in, del_list, axis = 0)
        AEW_lat_in = np.delete(AEW_lat_in, del_list, axis = 0)


        return AEW_lon_in, AEW_lat_in

    def extend_AEW(data_stuff, smooth_len): 
        #Test to make sure that we can add 13 points... otherwise, add other points
        data_in = data_stuff.copy()
        real_val = np.argwhere(~np.isnan(data_in))
        first = real_val[0][0]
        last = real_val[-1][0]
        data_len = len(data_in)

        if ((data_len-last)>=smooth_len) and ((first+1)>=smooth_len):
            append_len = smooth_len
        elif ((data_len-1)==last) or (first==0):
            print('Zero detected, will not attempt extension') #Don't do anything else
            return 0, 0, data_in
        else:
            append_len = np.min([data_len-last-1,first+1])-1
            print('Shorter Append = '+str(append_len))

        #We want to perform linear regression on first 10 points, then predict into the past
        #Then we want to perform linear regression for last 10 points,and predict into the future

        # Get linear fit for the first data points
        data_1 = data_in[~np.isnan(data_in)][0:append_len].reshape(-1,1)
        x1 = np.arange(len(data_1)).reshape(-1,1)
        x1_new = (np.arange(append_len) - append_len)
        x1_new = x1_new.reshape(-1,1)

        model1 = LinearRegression().fit(x1,data_1)
        reg_1 = model1.predict(x1_new)

        #print(append_len)
        #Get linear fit for last 5 data points
        data_2 = data_in[~np.isnan(data_in)][-append_len::].reshape(-1,1)
        x2 = np.arange(len(data_2)).reshape(-1,1)
        x2_new = (np.arange(append_len) + append_len)
        x2_new = x2_new.reshape(-1,1)

        model2 = LinearRegression().fit(x2,data_2)
        reg_2 = model2.predict(x2_new)

        #Finally, want to append the data before and after)
        #AEw_lon_in[first-]
        data_in = data_in.reshape(-1,1)
        data_in[(first-append_len):first] = reg_1
        data_in[(last+1):(last+append_len+1)] = reg_2
        return reg_1, reg_2, data_in.reshape(1,-1)

    #reg1, reg2, data_in = extend_AEW(data, smooth_len)

    def static_background(zoom = False, wider = False):
        #PROJECTION SETTINGS
        dataproj = ccrs.PlateCarree()

        #Define the extent of our static atlantic plot
        if zoom == True:
            lon1 = -70
            lon2 = 30
            lat1 = 5
            lat2 = 30
        if wider == True:
            lon1 = -110
            lon2 = 40
            lat1 = 0
            lat2 = 25
        else: #Modify to region of interest
            lon1 = -60
            lon2 = 0
            lat1 = 0
            lat2 = 40

        lon_delta = (lon2-lon1)/5
        lat_delta = (lat2-lat1)/5

        #Actually plot this
        fig=plt.figure(figsize=(15, 5))
        ax=plt.subplot(111, projection=dataproj)
        ax.set_extent([lon1, lon2, lat1, lat2],ccrs.PlateCarree())
        ax.coastlines('50m', linewidth=1.5)
        ax.add_feature(cfeature.STATES, linewidth=1.0)
        ax.add_feature(cfeature.BORDERS, linewidth=1.0)
        gl = ax.gridlines(color='gray',alpha=0.5,draw_labels=True) 
        gl.xlabels_top, gl.ylabels_right = False, False
        gl.xlabel_style, gl.ylabel_style = {'fontsize': 16}, {'fontsize': 16}
        gl.xlocator = mticker.FixedLocator([lon1,lon1+lon_delta, lon1+2*lon_delta, lon1+3*lon_delta, lon1+4*lon_delta, lon1+5*lon_delta])
        gl.ylocator = mticker.FixedLocator([lat1,lat1+lat_delta, lat1+2*lat_delta, lat1+3*lat_delta, lat1+4*lat_delta, lat1+5*lat_delta])
        gl.xformatter = LongitudeFormatter(zero_direction_label=True)
        gl.yformatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(gl.xformatter)
        ax.yaxis.set_major_formatter(gl.yformatter)
        ax.stock_img()

        return fig, ax

    def AEW_duplicate(lon_data, num_cutoff, close_cutoff, perc):
        lon_in = lon_data.copy()
        rm_list = []
        merge_list = []
        for row in range(np.shape(lon_in)[0]):
            print(row)
            if row in rm_list:
                continue #This means this row has already been flagged for removal so don't bother
            else:
                sample_data = lon_in[row,:] #Sample data for row currently in
                for new_row in range(row+1,np.shape(lon_in)[0]): #Look at the data for upcoming rows
                    new_data = lon_in[new_row,:]
                    equal_num = np.abs(new_data-sample_data)<=close_cutoff
                    equal_num = len(equal_num[equal_num == True])
                    #print(old_equal_num, equal_num)
                    if equal_num !=0:
                        first_i = np.where(np.abs(new_data-sample_data)<=close_cutoff)[0][0]
                        first_i2 = np.argwhere(new_data-sample_data<=close_cutoff)[0][0]
                        print(first_i, first_i2)
                        beg_new = new_data[:(first_i+1)]
                        beg_old = sample_data[:(first_i+1)]
                        beg_new = beg_new[~np.isnan(beg_new)]
                        beg_old = beg_old[~np.isnan(beg_old)]
                        prior_num_new = len(beg_new[~np.isnan(beg_new)])
                        prior_num_old = len(beg_old[~np.isnan(beg_old)])

                        # NEW EDITED SECTION.... TEST OUT
                        if prior_num_new != 0 and prior_num_old != 0:
                            if prior_num_new<=prior_num_old and (np.abs(beg_new[-1] - beg_new[0])<=5): #If the "new" data is greater or equal to 20, continue
                                print(beg_new, beg_old)
                                rm_list.append(new_row)
                            elif prior_num_old<=prior_num_new and (np.abs(beg_old[-1] - beg_old[0])<=5): #If the old data also fits these requirements along with new data...
                                print(beg_new, beg_old)
                                rm_list.append(row)

        return rm_list, merge_list
    # -----------------------------------------------------------------------------
    ##### END OF IMPORT STATEMENTS AND FUNCTIONS #####

    ##### MAIN CODE STARTS HERE ##### 
    # -----------------------------------------------------------------------------

    if len(sys.argv) == 1: #Basically, no inputs
        datafile = 'CURV_VORT/RADIAL_AVG/radial_avg_curv_vort.nc'
        #windfile = 'wind_for_tracking.nc'
    elif len(sys.argv) == 2: #Only input the data-in
        datafile = sys.argv[1]
        #windfile = 'wind_for_tracking.nc'
    elif len(sys.argv) == 3:
        data_in = sys.argv[1]
        data_out = sys.argv[2]
    else:
        raise Exception("Too many arguments provided.")

    # First, import the data
    nc_file = Dataset(datafile, 'r')
    time = nc_file.variables['time'][:]
    time_units = nc_file.variables['time'].units
    nclat = nc_file.variables['latitude'][:]
    nclon = nc_file.variables['longitude'][:]

    #Find, slice out only the area we are interested in (to reduce file size 
    #and prevent memory overuse/dumps!)
    # lon_st, n = find_nearest(nclon, -65)
    # lon_end, n = find_nearest(nclon, 60)
    # lat_end, n = find_nearest(nclat, 50)
    # lat_st, n = find_nearest(nclat, -20)

    # if lat_st> lat_end:
    #     lat = nclat[lat_end:lat_st]
    #     lon = nclon[lon_st:lon_end]
    # else:
    #     lat = nclat[lat_st:lat_end]
    #     lon = nclon[lon_st:lon_end]
    #print(lat, lon)
    lon = nclon
    lat = nclat
    curv_vort = nc_file.variables['curv_vort']#[:,lat_st:lat_end,lon_st:lon_end]

    #wind_file = Dataset(windfile, 'r')

    #wnd_lon = wind_file.variables['longitude'][:]
    #wnd_lat = wind_file.variables['latitude'][:]
    # wlon_st, n = find_nearest(wnd_lon, -60)
    # wlon_end, n = find_nearest(wnd_lon, 60)
    # wlat_end, n = find_nearest(wnd_lat, 50)
    # wlat_st, n = find_nearest(wnd_lat, -20)

    # if wlat_st> wlat_end:
    #     wlat = wnd_lat[wlat_end:wlat_st]
    #     wlon = wnd_lon[wlon_st:wlon_end]
    # else:
    #     wlat = wnd_lat[wlat_st:wlat_end]
    #     wlon = wnd_lon[wlon_st:wlon_end]
    #wlon = wnd_lon
    #wlat = wnd_lat

    #uwnd = wind_file.variables['u']#[:,wlat_st:wlat_end, wlon_st:wlon_end]
    #vwnd = wind_file.variables['v']#[:,wlat_st:wlat_end, wlon_st:wlon_end]

    #uwnd = wind_file.variables['u'][:,wlat_st:wlat_end, wlon_st:wlon_end]
    #vwnd = wind_file.variables['v'][:,wlat_st:wlat_end, wlon_st:wlon_end]

    dont_extra_list = []


    # To make sure this works properly, need to delete the variables ----------
    try:
        pass#del AEW_time, AEW_lon, AEW_lat
    except:
        print('Exception: Nothing to Delete')
    # ----------------- INITIALIZE IMPORTANT VARIABLES ------------------------
    mean_array = np.zeros((np.shape(time)[0],np.shape(lon)[0]))
    empty_stack = np.ones((1, np.shape(time)[0]))*np.NaN
    ##### MAIN ITERATION LOOP #####
    # -------------------------------------------------------------------------
    # ITERATE OVER EACH TIMESTEP IN THE GIVEN DATASET WE ARE TRACKING AEWS IN
    for slc_num in range(len(time)):
        print(str(slc_num+1)+' out of '+str(len(time)))

        #FIRST -- TAKE A SLICE OUT OF DATA FOR THE CURRENT FRAME
        curv_vort_data = curv_vort[slc_num, :, :]
        #PART I  ---- TEST FOR TWO IMPORTANT EXCEPTIONS -----------------------
        # -- FIRST ITERATION or
        # -- NO INITIATED WAVES EXIST
        # -- For these, we only care about trying to designate "new" waves
        if (slc_num == 0) or ('AEW_lon' not in locals()): #This only runs if an array has not already been created
            data_mean, init_max, cont_max, lat_center = find_maxima(curv_vort_data, lon, lat, lon_west, lon_east, thres_init, thres_cont, separate_bands = banding_t)
            mean_array[slc_num,:] = data_mean[:]
            #print(data_mean)
            #print(data_mean)
            if cleanup == True: #Cleanup is optional, relict of early tests. This only runs for initiation points
                new_cont_out = cont_max
                cont_out, init_out = execute_cleanup1(cont_max, init_max, data_mean, arb_thresh, deg_sep, res)
                new_init_out = init_out
                lon_init = lon[new_init_out]
                lon_cont = lon[new_cont_out]
            else:
                new_cont_out = cont_max
                new_init_out = init_max
                lon_init = lon[new_init_out]
                lon_cont = lon[new_cont_out]

            # GET ROUGH ESTIMATE OF LATITUDE POINTS -- we ONLY initate these from a center point (center of selected band)
            #This will find the latitude first guess as the center of the initating band (band with max value)
            lat_first_guess = []
            for idx in new_init_out:
                lati_append, n = find_nearest(lat, lat_center[idx])
                lat_first_guess.append(lati_append)

            # ------ Now, run centroid weighting to get the final list of cent_lat and cent_lons -----------                                    
            cent_lat = []
            cent_lon = []

            for i in range(len(lat_first_guess)): #Iterate over the total # of initial points we have
                lati = lat_first_guess[i] #Latitude first guess
                loni = new_init_out[i] #Longitude first guess
                lont, latt, pvout, pvmask, x_out, y_out = general_centroid(lon, lat, curv_vort_data, centroid_rad, lati, loni, centroid_it)
                cent_lat.append(y_out)
                cent_lon.append(x_out)

            # ------ Now -- do we have any initated AEWs? If so, let's build the main array ----------
            if cent_lon: #Will only run if something is in cent_lon... otherwise, move to next loop
                AEW_lon = np.ones((len(new_init_out), np.size(time)))*np.NaN
                AEW_lat = np.ones((len(new_init_out), np.size(time)))*np.NaN
                AEW_time = np.ones((len(new_init_out), np.size(time)))*np.NaN

                for row in range(len(cent_lon)): #Now, iterate each new AEW wave to put them in their own row
                    AEW_lon[row, slc_num] = cent_lon[row]
                    AEW_lat[row, slc_num] = cent_lat[row]
                    AEW_time[row, slc_num] = time[slc_num]
                    print(AEW_lon)
        else: 
            # IF WE ALREADY CREATED THE ARRAY... GOOD! But now we have to deal with a few other things
            #First, we need to get the new initation points. However, we only declare a new AEW if
            #the detected wave is NOT within a certain radius of a previous point!

            ##### EXTRAPOLATION SECTION: add in Alan Brammer's suggested extrapolation metric 

            data_mean, init_max, cont_max, lat_center = find_maxima(curv_vort_data, lon, lat, lon_west, lon_east, thres_init, thres_cont, separate_bands = banding_t)
            mean_array[slc_num,:] = data_mean[:]

            if EXTRAPOLATE == True:
                for existing in range(np.shape(AEW_lon)[0]): #First iterate over all the waves
                    AEW_lon_slc = AEW_lon[existing,:]
                    AEW_lat_slc = AEW_lat[existing,:]

                    if ~final_merge_list.any(): #Check to make sure we are not working with a merged wave
                        pass
                    elif existing in final_merge_list[:,1]:
                        continue
                    if AEW_lat[existing, (slc_num-1)]<=extra_lat_bottom:
                        continue

                    if len(AEW_lon_slc[~np.isnan(AEW_lon_slc)]) >= extrapolate_time and AEW_lon[existing, (slc_num-1)] <= extra_lon and AEW_lat[existing, (slc_num-1)]<=lat_max: #If track is long enough to do extrapolation, over ocean
                            #print(AEW_lon_slc[slc_num-1], AEW_lon[existing, slc_num-1])
                            # DO LINEAR EXTRAPOLATION BASED ON THE LINEAR FIT FOR LAST 5 TIMESTEPS

                            if AEW_lon[existing, (slc_num-1)] <= extra_transition:
                                extra_it = extra_it_60
                            else:
                                extra_it = extra_it_20

                            data_lon = AEW_lon_slc[~np.isnan(AEW_lon_slc)][-extra_it::].reshape(-1,1)
                            x = np.arange(len(data_lon)).reshape(-1,1)
                            x_new = np.array([len(data_lon)]) 
                            x_new = x_new.reshape(-1,1)

                            model = LinearRegression().fit(x,data_lon)
                            lon_guess_extra = data_lon[-1]+model.coef_

                            data_lat = AEW_lat_slc[~np.isnan(AEW_lon_slc)][-extra_it::].reshape(-1,1)
                            x2 = np.arange(len(data_lat)).reshape(-1,1)
                            x2_new = np.array([len(data_lat)])
                            x2_new = x2_new.reshape(-1,1)

                            model2 = LinearRegression().fit(x2,data_lat)
                            lat_guess_extra = data_lat[-1]+model2.coef_

                            lon_guess_i, n = find_nearest(lon, lon_guess_extra)
                            lat_guess_i, n = find_nearest(lat, lat_guess_extra)
                            #print(lon_guess_extra, lat_guess_extra)

                            # RUN A CENTROID FINDER FOR THE NEXT TIME STEP
                            if AEW_lon[existing, (slc_num-1)] <= extra_transition:
                                extra_dist = extra_distance_60
                                centroid_weight = centroid_rad_extra_60
                            else:
                                extra_dist = extra_distance
                                centroid_weight = centroid_rad_extra

                            try:
                                lont_off, latt_off, pvout_off, pvmask_off, lon_extra_out, lat_extra_out = general_centroid(lon, lat, curv_vort_data, centroid_weight, lat_guess_i, lon_guess_i, centroid_it_extra)
                            except:
                                print('Something went wrong with centroid, possibly out of bounds')
                                continue
                            # FIND THE CURVATURE VORTICITY VALUE AT THIS NEXT CENTER
                            lon_extra_posi, n = find_nearest(lon, lon_extra_out)
                            lat_extra_posi, n = find_nearest(lat, lat_extra_out)

                            curv_extra = curv_vort_data[lat_extra_posi, lon_extra_posi]


                            #FIND DISTANCE
                            win3_extra, distance3_extra, fwd_extra = within_distance_direct(lon, lat, lon_extra_out,lat_extra_out,
                                                                               AEW_lon[existing, (slc_num-1)], AEW_lat[existing, (slc_num-1)],
                                                                               extra_dist)
                            #replace with "guess_lon_extra", "guess_lat_extra" for difference

                            if win3_extra == True and curv_extra >=extrapolate_thresh:
                                if EXTRAPOLATE_FORWARD == True:
                                    backward_dist = haversine(lon_extra_out, lat_extra_out, AEW_lon[existing, (slc_num-1)], lat_extra_out) 
                                    curv_lon_extra,n = find_nearest(lon, AEW_lon[existing, slc_num-1])
                                    curv_lat_extra,n = find_nearest(lat, AEW_lat[existing, slc_num-1]) 
                                    curv_val_extra = curv_vort_data[curv_lat_extra, curv_lon_extra]
                                    if lon_extra_out>= AEW_lon[existing, (slc_num-1)] and backward_dist >= extra_backcut and AEW_lat[existing, (slc_num-1)]<=extra_lat_start and curv_val_extra< speed_curv_thresh:
                                        continue
                                AEW_lon[existing, slc_num] = lon_extra_out
                                AEW_lat[existing, slc_num] = lat_extra_out
                    elif EXTRAPOLATE_SKIP == True:
                        if len(AEW_lon_slc[~np.isnan(AEW_lon_slc)]) >= extrapolate_time and AEW_lon[existing, (slc_num-2)] <= extra_lon and AEW_lat[existing, (slc_num-2)]<=lat_max: #If track is long enough to do extrapolation, over ocean
                            if AEW_lon[existing, (slc_num-2)] <= extra_transition:
                                extra_it = extra_it_60
                            else:
                                extra_it = extra_it_20

                            data_lon = AEW_lon_slc[~np.isnan(AEW_lon_slc)][-(extra_it+1)::].reshape(-1,1)
                            x = np.arange(len(data_lon)).reshape(-1,1)
                            x_new = np.array([len(data_lon)]) 
                            x_new = x_new.reshape(-1,1)

                            model = LinearRegression().fit(x,data_lon)
                            lon_guess_extra = data_lon[-1]+2*model.coef_

                            data_lat = AEW_lat_slc[~np.isnan(AEW_lon_slc)][-(extra_it+1)::].reshape(-1,1)
                            x2 = np.arange(len(data_lat)).reshape(-1,1)
                            x2_new = np.array([len(data_lat)])
                            x2_new = x2_new.reshape(-1,1)

                            model2 = LinearRegression().fit(x2,data_lat)
                            lat_guess_extra = data_lat[-1]+2*model2.coef_

                            lon_guess_i, n = find_nearest(lon, lon_guess_extra)
                            lat_guess_i, n = find_nearest(lat, lat_guess_extra)
                            #print(lon_guess_extra, lat_guess_extra)

                            # RUN A CENTROID FINDER FOR THE NEXT TIME STEP
                            if AEW_lon[existing, (slc_num-2)] <= extra_transition:
                                extra_dist = extra_distance_60
                                centroid_weight = centroid_rad_extra_60
                            else:
                                extra_dist = extra_distance
                                centroid_weight = centroid_rad_extra
                            try:
                                lont_off, latt_off, pvout_off, pvmask_off, lon_extra_out, lat_extra_out = general_centroid(lon, lat, curv_vort_data, centroid_weight, lat_guess_i, lon_guess_i, centroid_it_extra)
                            except:
                                print('Something went wrong with centroid, possibly out of bounds')
                                continue
                            # FIND THE CURVATURE VORTICITY VALUE AT THIS NEXT CENTER
                            lon_extra_posi, n = find_nearest(lon, lon_extra_out)
                            lat_extra_posi, n = find_nearest(lat, lat_extra_out)

                            curv_extra = curv_vort_data[lat_extra_posi, lon_extra_posi]

                            if lon_guess_extra<=-60:
                                extra_dist2_final = extra_dist2_60
                            else:
                                extra_dist2_final = extra_dist2
                            #FIND DISTANCE
                            win3_extra, distance3_extra, fwd_extra = within_distance_direct(lon, lat, lon_extra_out,lat_extra_out,
                                                                               lon_guess_extra, lat_guess_extra,
                                                                               extra_dist2_final)
                            #replace with "guess_lon_extra", "guess_lat_extra" for difference
                            center_distance = haversine(lon_extra_out, lat_extra_out, AEW_lon[existing, (slc_num-2)], AEW_lat[existing, (slc_num - 2)])

                            if win3_extra == True and curv_extra >=extrapolate_thresh and center_distance<=extra_dist2_final and lon_extra_out<=AEW_lon[existing, (slc_num-3)]:
                            #if curv_extra >=extrapolate_thresh:
                                #print('Boom')
                                AEW_lon[existing, slc_num] = lon_extra_out
                                AEW_lat[existing, slc_num] = lat_extra_out

                        elif len(AEW_lon_slc[~np.isnan(AEW_lon_slc)]) >= extrapolate_time and AEW_lon[existing, (slc_num-3)] <= extra_lon and AEW_lat[existing, (slc_num-3)]<=lat_max: #If track is long enough to do extrapolation, over ocean
                            if AEW_lon[existing, (slc_num-2)] <= extra_transition:
                                extra_it = extra_it_60
                            else:
                                extra_it = extra_it_20

                            data_lon = AEW_lon_slc[~np.isnan(AEW_lon_slc)][-(extra_it+1)::].reshape(-1,1)
                            x = np.arange(len(data_lon)).reshape(-1,1)
                            x_new = np.array([len(data_lon)]) 
                            x_new = x_new.reshape(-1,1)

                            model = LinearRegression().fit(x,data_lon)
                            lon_guess_extra = data_lon[-1]+3*model.coef_

                            data_lat = AEW_lat_slc[~np.isnan(AEW_lon_slc)][-(extra_it+1)::].reshape(-1,1)
                            x2 = np.arange(len(data_lat)).reshape(-1,1)
                            x2_new = np.array([len(data_lat)])
                            x2_new = x2_new.reshape(-1,1)

                            model2 = LinearRegression().fit(x2,data_lat)
                            lat_guess_extra = data_lat[-1]+3*model2.coef_

                            lon_guess_i, n = find_nearest(lon, lon_guess_extra)
                            lat_guess_i, n = find_nearest(lat, lat_guess_extra)
                            #print(lon_guess_extra, lat_guess_extra)

                            if AEW_lon[existing, (slc_num-1)] <= extra_transition:
                                extra_dist = extra_distance_60
                                centroid_weight = centroid_rad_extra_60
                            else:
                                extra_dist = extra_distance
                                centroid_weight = centroid_rad_extra
                            # RUN A CENTROID FINDER FOR THE NEXT TIME STEP
                            try:
                                lont_off, latt_off, pvout_off, pvmask_off, lon_extra_out, lat_extra_out = general_centroid(lon, lat, curv_vort_data, centroid_weight, lat_guess_i, lon_guess_i, centroid_it_extra)
                            except:
                                print('Something went wrong with centroid, possibly out of bounds')
                                continue
                            # FIND THE CURVATURE VORTICITY VALUE AT THIS NEXT CENTER
                            lon_extra_posi, n = find_nearest(lon, lon_extra_out)
                            lat_extra_posi, n = find_nearest(lat, lat_extra_out)

                            curv_extra = curv_vort_data[lat_extra_posi, lon_extra_posi]

                            if lon_guess_extra<=-60:
                                extra_dist3_final = extra_dist3_60
                            else:
                                extra_dist3_final = extra_dist3

                            #FIND DISTANCE
                            win3_extra, distance3_extra, fwd_extra = within_distance_direct(lon, lat, lon_extra_out,lat_extra_out,
                                                                               lon_guess_extra, lat_guess_extra,
                                                                               extra_dist3_final)
                            #replace with "guess_lon_extra", "guess_lat_extra" for difference
                            center_distance = haversine(lon_extra_out, lat_extra_out, AEW_lon[existing, (slc_num-3)], AEW_lat[existing, (slc_num - 3)])

                            if win3_extra == True and curv_extra >=extrapolate_thresh and center_distance<=extra_dist3_final and lon_extra_out<=AEW_lon[existing, (slc_num-3)]:
                            #if curv_extra >=extrapolate_thresh:
                                #print('Boom')
                                AEW_lon[existing, slc_num] = lon_extra_out
                                AEW_lat[existing, slc_num] = lat_extra_out
                    #Now if it still doesn't append... prevent the code from actually doing the other method
                    if ~np.isnan(AEW_lon[existing, slc_num]) and existing not in dont_extra_list and AEW_lon[existing, slc_num]<-60:
                        dont_extra_list.append(existing)



            if cleanup == True: #Cleanup is optional, relict of early tests. This only runs for initiation points.
                new_cont_out = cont_max
                cont_out, init_out = execute_cleanup1(cont_max, init_max, data_mean, arb_thresh, deg_sep, res)
                new_init_out = init_out
                lon_init = lon[new_init_out]
                lon_cont = lon[new_cont_out]
            else:
                new_cont_out = cont_max
                new_init_out = init_max
                lon_init = lon[new_init_out]
                lon_cont = lon[new_cont_out]

            # GET ROUGH ESTIMATE OF LATITUDE POINTS -- we ONLY initate these from a center point (center of selected band)
            lat_first_guess = []
            for idx in new_init_out:
                lati_append, n = find_nearest(lat, lat_center[idx])
                lat_first_guess.append(lati_append)

            # ------ Now, run centroid testing to get the final list of cent_lat and cent_lons -----------                                    
            cent_lat = []
            cent_lon = []
            for i in range(len(lat_first_guess)): #Iterate over the total # of initial points we have
                lati = lat_first_guess[i]
                loni = new_init_out[i]
                lont, latt, pvout, pvmask, x_out, y_out = general_centroid(lon, lat, curv_vort_data, centroid_rad, lati, loni, centroid_it)
                cent_lat.append(y_out)
                cent_lon.append(x_out)

            # ----- DISTANCE TESTING FOR INITATED POINTS -----
            #NOW DO WE ADD THESE NEW POINTS TO THE LIST? DEPENDS ON IF THEY ARE WITHIN A RADIUS OF THE POINTS IN QUESTION
            for init_point in range(len(cent_lon)): #For each potential initation point  
                win3_list = []
                win6_list = []
                for existing in range(np.shape(AEW_lon)[0]): #Check with each existing AEW track to see if it is within a distance of them
                    win3, distance3, fwd_init = within_distance_direct(lon, lat, cent_lon[init_point],cent_lat[init_point],
                                                                       AEW_lon[existing, (slc_num-1)],AEW_lat[existing, (slc_num-1)],
                                                                       step_3hr_init)
                    win3_list.append(win3)
                    try:
                        win6, distance6, fwd_init = within_distance_direct(lon, lat, cent_lon[init_point],cent_lat[init_point],
                                                                           AEW_lon[existing, (slc_num-2)],AEW_lat[existing, (slc_num-2)],
                                                                           step_6hr_init)
                        win6_list.append(win6)
                    except:
                        if slc_num < 3:
                            print('Exception: Likely normal, not enough timesteps to run code')
                        else:
                            print('Exception: Likely error in code. Check.')
                #Now, if there are ANY "True" values in the win3 array, we don't want to include point. Otherwise, we do.
                #if not any(win6_list) and not any(win3_list) and not any(win9_list) and not np.isnan(cent_lon[init_point]):
                if not any(win3_list) and not any(win6_list) and not np.isnan(cent_lon[init_point]):
                    AEW_lon = np.vstack((AEW_lon, empty_stack)) #Add another row for the new AEW
                    AEW_lat = np.vstack((AEW_lat, empty_stack)) #Add antoher row for the new AEW

                    AEW_lon[-1, slc_num] = cent_lon[init_point]
                    AEW_lat[-1, slc_num] = cent_lat[init_point]
            #NEXT -- ADDING POINTS TO EXISTING AEW POINTS
            #Now we have to run centroid finding and append the new lat/lons to another list for further testing
            lat_point_array = []
            lon_point_array = []
            lat_first_guess_cont = []
            for idx in new_cont_out:
                lati_append, n = find_nearest(lat, lat_center[idx])
                lat_first_guess_cont.append(lati_append)
            for i in range(len(new_cont_out)): 
                loni = new_cont_out[i]
                lati = lat_first_guess_cont[i]

                try:
                    lont, latt, pvout, pvmask, lon_new_pt, lat_new_pt = general_centroid(lon, lat, curv_vort_data, centroid_rad, lati, loni, centroid_it)
                    if upper_limit == True: #Here we make all points that are centered north of 20N and limit them to just 20
                        if lat_new_pt > 20 and lat_new_pt <=25: #Only allow points that are somewhat over the threshold, don't get too crazy
                            lat_new_pt = 20 #Set it equal to 20N now
                            print("Upper Limit Engaged")
                    lat_point_array.append(lat_new_pt)
                    lon_point_array.append(lon_new_pt)


                except: 
                    print('Exception -- something major went wrong and not sure why.')
                    continue


        #NOW DO WE ADD THESE NEW POINTS TO THE LIST? DEPENDS ON IF THEY ARE WITHIN A RADIUS OF THE POINTS IN QUESTION
            for existing in range(np.shape(AEW_lon)[0]): #Check with each 
                if ~final_merge_list.any(): #Check to make sure we are not working with a merged wave
                    pass
                elif existing in final_merge_list[:,1]:
                    continue
                if existing in dont_extra_list: #Check to see if wave has been extrapolated... if so, don't append
                    continue

                if np.isnan(AEW_lon[existing,slc_num]):
                    temp_dist_list3 = []
                    temp_dist_list6 = []
                    temp_fwd_list3 = []
                    temp_fwd_list6 = []

                    temp_lon_arr3 = []
                    temp_lat_arr3 = []
                    temp_lon_arr6 = []
                    temp_lat_arr6 = []

                    for cont_i in range(len(lon_point_array)): #For each potential initation point 
                        temp_bool3, temp_dist3, fwd3 = within_distance_direct(lon, lat, lon_point_array[cont_i], lat_point_array[cont_i],
                                                                     AEW_lon[existing, (slc_num-1)], AEW_lat[existing, (slc_num-1)],
                                                                       step_3hr)
                        temp_bool6, temp_dist6, fwd6 = within_distance_direct(lon, lat, lon_point_array[cont_i], lat_point_array[cont_i],
                                                                     AEW_lon[existing, (slc_num-2)], AEW_lat[existing, (slc_num-2)],
                                                                       step_6hr)



                        if temp_bool3 == True: #If it is within the specified range 
                            temp_dist_list3.append(temp_dist3)
                            temp_fwd_list3.append(fwd3)
                            temp_lon_arr3.append(lon_point_array[cont_i])
                            temp_lat_arr3.append(lat_point_array[cont_i])
                        elif temp_bool6 == True: #Same for timesteps 6 hours ago
                            temp_dist_list6.append(temp_dist6)
                            temp_fwd_list6.append(fwd6)
                            temp_lon_arr6.append(lon_point_array[cont_i])
                            temp_lat_arr6.append(lat_point_array[cont_i])


                    #Now we have a list of points that satisfy the distance requirements for this given AEW. What now?
                    #We determine the point the closest to the AEW from one timestep ago. If they don't exist, we try the same for 
                    #6 hourly.

                    #Caveat -- we don't necessarily prevent points from being fixed to two waves. Will see how this performs and
                    #will adjust accordingly. Could be okay to identify wave mergers. 
                    temp_dist_list3 = np.array(temp_dist_list3)
                    temp_fwd_list3 = np.array(temp_fwd_list3)
                    temp_lat_list3 = np.array(temp_lat_arr3)
                    temp_lon_list3 = np.array(temp_lon_arr3)

                    if any(temp_dist_list3): #If there are any in the 1 timestep list
                        min_i = np.argmin(temp_dist_list3)   
                        #Test to see if latitudinal jump occurs
                        if np.abs(temp_lat_arr3[min_i]-AEW_lat[existing, (slc_num-1)])>= land_lat_limit:
                            continue

                        #Get some curvature value characteristics for later
                        curv_lon_back,n = find_nearest(lon, temp_lon_list3[min_i])
                        curv_lat_back,n = find_nearest(lat, temp_lat_list3[min_i]) 

                        curv_val_back = curv_vort_data[curv_lat_back, curv_lon_back]

                        if force_bump == True:
                            if temp_fwd_list3[min_i] == False and any(temp_fwd_list3): #If the shortest distance is backwards and there exist points that are forward within the threshold..
                                new_temp_fwd_list = temp_fwd_list3[temp_fwd_list3 == True]
                                new_temp_lon_list = temp_lon_list3[temp_fwd_list3 == True]
                                new_temp_lat_list = temp_lat_list3[temp_fwd_list3 == True]
                                new_temp_dist_list = temp_dist_list3[temp_fwd_list3 == True]

                                min_adj_i = np.argmin(new_temp_dist_list)
                                if new_temp_dist_list[min_adj_i] <= bump_num*step_3hr:
                                    AEW_lon[existing, slc_num] = new_temp_lon_list[min_adj_i]
                                    AEW_lat[existing, slc_num] = new_temp_lat_list[min_adj_i]
                                    continue
                        ##Different backward cutoffs over ocean vs. land
                        if AEW_lon[existing, (slc_num-1)] <= -17: #If over the ocean
                            back_cutoff = back_cutoff_ocean
                            back_cutoff_long = back_cutoff_long_ocean
                        else:
                            back_cutoff = back_cutoff_land
                            back_cutoff_long = back_cutoff_long_land
                        fwd = temp_fwd_list3[min_i]
                        if force_forward == True: 
                            if fwd == False:
                                if AEW_lon[existing, (slc_num-1)] != AEW_lon[existing, (slc_num-stuck_thresh)] or cut_stuck == False:
                                    if temp_dist_list3[min_i]<= back_cutoff:
                                        AEW_lon[existing, slc_num] = AEW_lon[existing, (slc_num-1)]
                                        AEW_lat[existing, slc_num] = temp_lat_arr3[min_i]
                            else:   
                                if temp_lon_arr3[min_i] == AEW_lon[existing, (slc_num-stuck_thresh)]: #NEW, TEST FOR STUCK THRESH
                                    continue
                                AEW_lon[existing, slc_num] = temp_lon_arr3[min_i]
                                AEW_lat[existing, slc_num] = temp_lat_arr3[min_i]
                        else:
                            AEW_lon[existing, slc_num] = temp_lon_arr3[min_i]
                            AEW_lat[existing, slc_num] = temp_lat_arr3[min_i]

                    elif any(temp_dist_list6) and np.isnan(AEW_lon[existing,slc_num-1]): #If not, are there any on the 2 timestep list?
                        min_i = np.argmin(temp_dist_list6)
                        if AEW_lon[existing, (slc_num-2)] <= -17: #If over the ocean
                            back_cutoff = back_cutoff_ocean
                            back_cutoff_long = back_cutoff_long_ocean
                        else:
                            back_cutoff = back_cutoff_land
                            back_cutoff_long = back_cutoff_long_land
                            #TESTING CUT OUT BELOW IF TOO RESTRICTIVE
                        if np.abs(temp_lat_arr6[min_i]-AEW_lat[existing, (slc_num-2)])>= land_lat_limit*2:
                            continue
                        fwd = temp_fwd_list6[min_i]
                        if force_forward == True:
                            if fwd == False:
                                if AEW_lon[existing, (slc_num-2)] != AEW_lon[existing, (slc_num-stuck_thresh)] or cut_stuck == False:
                                    if temp_dist_list6[min_i]<= back_cutoff_long and long_edit==True:
                                        AEW_lon[existing, slc_num] = AEW_lon[existing, (slc_num-2)]
                                        AEW_lat[existing, slc_num] = temp_lat_arr6[min_i]
                            else:                    
                                AEW_lon[existing, slc_num] = temp_lon_arr6[min_i]
                                AEW_lat[existing, slc_num] = temp_lat_arr6[min_i]
                        else:
                            AEW_lon[existing, slc_num] = temp_lon_arr6[min_i]
                            AEW_lat[existing, slc_num] = temp_lat_arr6[min_i]     
        if merge_wave == True:
             for existing in range(np.shape(AEW_lon)[0]): #First iterate over all the waves
                        AEW_lon_slc = AEW_lon[existing,:]
                        AEW_lat_slc = AEW_lat[existing,:]

                        if len(AEW_lon_slc[~np.isnan(AEW_lon_slc)]) <= extrapolate_time:
                            continue

                        for new_existing in range(existing+1,np.shape(AEW_lon)[0]):
                            new_AEW_lon_slc = AEW_lon[new_existing, :]
                            new_AEW_lat_slc = AEW_lat[new_existing, :]


                            lon1 = AEW_lon_slc[slc_num]
                            lat1 = AEW_lat_slc[slc_num]
                            lon2 = new_AEW_lon_slc[slc_num]
                            lat2 = new_AEW_lat_slc[slc_num]

                            two_wave_distance = haversine(lon1, lat1, lon2, lat2)
                            new_wave_data = new_AEW_lon_slc[~np.isnan(new_AEW_lon_slc)]
                            try: 
                                new_wave_lon = np.abs(new_wave_data[-1] - new_wave_data[0])
                            except: 
                                new_wave_lon = 0

                            if len(new_AEW_lon_slc[~np.isnan(new_AEW_lon_slc)]) <= extrapolate_time and two_wave_distance<merge_distance and new_wave_lon<10:
                                AEW_lon[new_existing, :] = np.ones(np.shape(AEW_lon[0,:]))*np.NaN
                                AEW_lat[new_existing, :] = np.ones(np.shape(AEW_lon[0,:]))*np.NaN
                                continue
                            #print(two_wave_distance)
                            if ~np.isnan(two_wave_distance) and two_wave_distance < merge_distance: #If wave is found to be mergin
                                data_lon = AEW_lon_slc[~np.isnan(AEW_lon_slc)][-(extra_it+1)::].reshape(-1,1)
                                #The wave is merging! We keep the wave that is moving forward the fastest

                                lon_speed1 = np.mean(np.diff(AEW_lon_slc[~np.isnan(AEW_lon_slc)][-extra_it_60:]))
                                lon_speed2 = np.mean(np.diff(new_AEW_lon_slc[~np.isnan(AEW_lon_slc)][-extra_it_60:]))

                                #Merge list -- first column is kept, second column is not
                                if np.abs(lon_speed1) >= np.abs(lon_speed2):
                                    temp_merge = [existing, new_existing]
                                    AEW_lon[new_existing, slc_num] = AEW_lon[existing, slc_num]
                                    AEW_lat[new_existing, slc_num] = AEW_lat[existing, slc_num]
                                else:
                                    temp_merge = [new_existing, existing]
                                    AEW_lon[existing, slc_num] = AEW_lon[new_existing, slc_num]
                                    AEW_lat[existing, slc_num] = AEW_lat[new_existing, slc_num]           

                                #Finally, concatinate 
                                if not final_merge_list.any():
                                    final_merge_list = np.array(temp_merge).reshape(1,2)
                                else:
                                    final_merge_list = np.vstack((final_merge_list, temp_merge))

        if (speed_limit == True) and 'AEW_lon' in locals():
             for existing in range(np.shape(AEW_lon)[0]): #First iterate over all the waves
                        AEW_lon_slc = AEW_lon[existing,:]
                        AEW_lat_slc = AEW_lat[existing,:]

                        if np.isnan(AEW_lon_slc[slc_num]):
                            continue

                        curv_lon,n = find_nearest(lon, AEW_lon_slc[slc_num])
                        curv_lat,n = find_nearest(lat, AEW_lat_slc[slc_num]) 

                        curv_val = curv_vort_data[curv_lat, curv_lon]

                        if (np.abs(AEW_lon_slc[slc_num]-AEW_lon_slc[slc_num-stuck_thresh])<= stuck_lon or (AEW_lon_slc[slc_num] - AEW_lon_slc[slc_num-stuck_thresh])> stuck_lon) and AEW_lon_slc[slc_num]<= extra_lon and AEW_lat_slc[slc_num]<extra_lat_start and curv_val <= speed_curv_thresh:
                            AEW_lon[existing, slc_num] = np.NaN
                            AEW_lat[existing, slc_num] = np.NaN
                            AEW_lon[existing, (slc_num-1)] = np.NaN
                            AEW_lat[existing, (slc_num-1)] = np.NaN                                

    # --------------------------------------------------------------------------
    #END OF MAIN ITREATION LOOP

    ##### DETELE EMPTY DATASETS, DUPLICATES, AND RUN ANY LINEAR EXTENSIONS
    # --------------------------------------------------------------------------
    del_list_first = []
    for row in range(np.shape(AEW_lon)[0]):
        data = AEW_lon[row,:]
        if np.isnan(data).all():
            del_list_first.append(row)
    AEW_lon = np.delete(AEW_lon, del_list_first, axis = 0)
    AEW_lat = np.delete(AEW_lat, del_list_first, axis = 0)

    new_AEW_lon, new_AEW_lat = cleanup_AEW(AEW_lon, AEW_lat, temp_res, days_remove)
    AEW_lat = new_AEW_lat
    AEW_lon = new_AEW_lon

    #Run smoothing savgol filter
    reg_list = [num2date(i, time_units, only_use_cftime_datetimes = False) for i in time]
    data_slc = 2

    if duplicate_removal == True: 
        rm_dup_list, merge_list = AEW_duplicate(AEW_lon,1, 1, 0.7)

        AEW_lon= np.delete(AEW_lon, rm_dup_list, axis = 0)
        AEW_lat = np.delete(AEW_lat, rm_dup_list, axis = 0)

        #Adjust merge list values to be accurate
        for i in range(len(merge_list)):
            reduce_by = len(np.array(rm_dup_list)<merge_list[i])
            merge_list[i] = merge_list[i] - reduce_by
        #print(merge_list)

    AEW_lon_filter = savgol_filter(AEW_lon, smooth_len, 2, mode = 'nearest')
    AEW_lat_filter = savgol_filter(AEW_lat, smooth_len, 2, mode = 'nearest')
    #plt.plot(AEW_lon_filter[data_slc,:], AEW_lat_filter[data_slc,:], 'b')

    #Run Linear Extender to replace missing times after filtering
    if run_extender == True:
        for row in range(np.shape(AEW_lon)[0]):
            data_lon_filter = AEW_lon_filter[row,:]
            data_lat_filter = AEW_lat_filter[row,:]
            reg11,reg22, data_out_lon = extend_AEW(data_lon_filter, int(np.floor(smooth_len/3)))
            reg12,reg22, data_out_lat = extend_AEW(data_lat_filter, int(np.floor(smooth_len/3)))
            AEW_lon_filter[row,:] = data_out_lon
            AEW_lat_filter[row,:] = data_out_lat

    reg_list = [num2date(i, time_units, only_use_cftime_datetimes = False) for i in time]

    # -------------------------------------------------------------------------
    # END DATA CLEANUP
    ##### SAVE OUTPUT NETCDF4 DATAFILE
    # -------------------------------------------------------------------------
    if save_data == True:
        outfile = 'TRACKING/prelim_aew_tracks_out.nc'

        nlat = np.size(lat)
        nlon = np.size(lon)
        nsys = np.shape(AEW_lon)[0]
        system_in = np.arange(nsys)+1

        # open a netCDF file to write
        ncout = Dataset(outfile, 'w', format='NETCDF4')

        # define axis size
        ncout.createDimension('time', None)  # unlimited
        ncout.createDimension('system', nsys)
        ncout.createDimension('latitude', nlat)
        ncout.createDimension('longitude', nlon)


        # create latitude axis
        latitude = ncout.createVariable('latitude', np.dtype('float64').char, ('latitude'))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'latitude'
        latitude.units = 'degrees_north'
        latitude.axis = 'Y'

        # create longitude axis
        longitude = ncout.createVariable('longitude', np.dtype('float64').char, ('longitude'))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'longitude'
        longitude.units = 'degrees_east'
        longitude.axis = 'X'

        # create time axis
        time_data = ncout.createVariable('time', np.dtype('float64').char, ('time'))
        time_data.long_name = 'time'
        time_data.units = time_units
        time_data.calendar = 'gregorian'
        time_data.axis = 'T'

        # create system
        system_data = ncout.createVariable('system', np.dtype('float64').char, ('system'))
        system_data.standard_name = 'system'
        system_data.long_name = 'system'
        system_data.units = 'number'
        system_data.axis = 'Y'


        # create variable array
        AEW_longitude = ncout.createVariable('AEW_lon', np.dtype('float64').char, ('system', 'time'))
        AEW_longitude.long_name = 'Longitude of AEW Tracks'
        AEW_longitude.units = 'degrees'

        AEW_latitude = ncout.createVariable('AEW_lat', np.dtype('float64').char, ('system', 'time'))
        AEW_latitude.long_name = 'Latitude of AEW Tracks'
        AEW_latitude.units = 'degrees'

        AEW_longitude_smooth = ncout.createVariable('AEW_lon_smooth', np.dtype('float64').char, ('system', 'time'))
        AEW_longitude_smooth.long_name = 'Smoothed Longitude of AEW Tracks'
        AEW_longitude_smooth.units = 'degrees'

        AEW_latitude_smooth = ncout.createVariable('AEW_lat_smooth', np.dtype('float64').char, ('system', 'time'))
        AEW_latitude_smooth.long_name = 'Smoothed Latitude of AEW Tracks'
        AEW_latitude_smooth.units = 'degrees'

        curv_data_mean = ncout.createVariable('curv_data_mean', np.dtype('float64').char, ('time','longitude'))
        curv_data_mean.long_name = 'Averaged Non-Divergent Curvature Vorticity (5-20N)'
        curv_data_mean.units = 's**-1' 

        # copy axis from original dataset
        time_data[:] = time[:]
        system_data[:] = system_in[:]
        longitude[:] = lon[:]
        latitude[:] = lat[:]
        AEW_longitude[:,:] = AEW_lon[:,:]
        AEW_latitude[:,:] = AEW_lat[:,:]
        AEW_longitude_smooth[:,:] = AEW_lon_filter[:,:]
        AEW_latitude_smooth[:,:] = AEW_lat_filter[:,:]
        curv_data_mean[:,:] = mean_array[:,:]
        ncout.close()

    # -------------------------------------------------------------------------
    # END DATA SAVE

    
