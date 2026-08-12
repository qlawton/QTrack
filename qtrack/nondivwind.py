#!/usr/bin/env python3


def compute_nondiv_wind(input_file, output_file="wind_700_helmholtz.nc", lon_bounds=None, lat_bounds=None):
    """
    Compute non-divergent component of wind for the 700hPa level. Note that this requires a full global grid, as it utilizes spherical harmonics (windspharm: doi.org/10.5334/jors.129). Not required for tracking, and can be skipped.

    Parameters
    ----------
    lon_bounds, lat_bounds : (min, max), optional
        Crop the output to these bounds. Both default to None, i.e. keep the full
        global grid the spherical harmonics were computed on. This used to be a
        hard-coded crop to 120W-60E and 20S-60N, which discarded the entire
        Pacific and Indian Ocean and so made global tracking impossible.
    """

    import matplotlib as mpl
    import numpy as np
    import xarray as xr
    from netCDF4 import Dataset
    from numpy import dtype
    from windspharm.standard import VectorWind
    from windspharm.tools import order_latdim, prep_data, recover_data

    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx, array[idx]

    mpl.rcParams["mathtext.default"] = "regular"

    file_in = input_file
    outfile = output_file

    SAVE_OUTPUT = True
    # ----- READ IN THE DATA -----

    def non_divergent_wind(file_in):
        nc_wind_file = file_in
        wind_data = xr.open_dataset(nc_wind_file).fillna(0)
        nc_load = Dataset(nc_wind_file, "r")

        # Read zonal and meridional wind components from file using the netCDF4
        # module. Here I will use the 700-hPa level since it will be useful for
        # AEW wave tracking
        uwnd = wind_data["u"].values
        lons = wind_data["longitude"].values
        lats = wind_data["latitude"].values
        vwnd = wind_data["v"].values
        tm = nc_load.variables["time"][:]
        tm_units = nc_load.variables["time"].units

        # The standard interface requires that latitude and longitude be the leading
        # dimensions of the input wind components, and that wind components must be
        # either 2D or 3D arrays. The data read in is 3D and has latitude and
        # longitude as the last dimensions. The bundled tools can make the process of
        # re-shaping the data a lot easier to manage.
        uwnd, uwnd_info = prep_data(uwnd, "tyx")
        vwnd, vwnd_info = prep_data(vwnd, "tyx")

        # It is also required that the latitude dimension is north-to-south. Again the
        # bundled tools make this easy.
        lats, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)

        # Create a VectorWind instance to handle the computation of streamfunction and
        # velocity potential.
        w = VectorWind(uwnd, vwnd)
        del uwnd, vwnd
        # Compute the streamfunction and velocity potential. Also use the bundled
        # tools to re-shape the outputs to the 4D shape of the wind components as they
        # were read off files.
        upsi, vpsi = w.nondivergentcomponent()
        upsi = recover_data(upsi, uwnd_info)
        vpsi = recover_data(vpsi, vwnd_info)
        # upsi = recover_data(upsi, uwnd_info)
        # vpsi = recover_data(vpsi, vwnd_info)
        # del(uwnd_info, vwnd_info, w, uwnd, vwnd)

        # OPTIONALLY CUT DOWN DATA TO A REGION OF INTEREST
        # Note that cropping longitude leaves a non-periodic grid, which stops the
        # tracker following waves past the crop edges. Leave both bounds at None
        # (the default) for global tracking.
        if lon_bounds is not None:
            lon_st, n = find_nearest(lons, lon_bounds[0])
            lon_end, n = find_nearest(lons, lon_bounds[-1])
            lons = lons[lon_st:lon_end]
            upsi = upsi[:, :, lon_st:lon_end]
            vpsi = vpsi[:, :, lon_st:lon_end]
        if lat_bounds is not None:
            # lats run north to south here, so the max bound comes first
            lat_st, n = find_nearest(lats, max(lat_bounds))
            lat_end, n = find_nearest(lats, min(lat_bounds))
            lats = lats[lat_st:lat_end]
            upsi = upsi[:, lat_st:lat_end, :]
            vpsi = vpsi[:, lat_st:lat_end, :]

        # %% WRITE NETCDF4
        # ------------------
        # write netCDF file
        # ------------------

        if SAVE_OUTPUT:
            nlat = np.size(lats)
            nlon = np.size(lons)
            # open a netCDF file to write
            ncout = Dataset(outfile, "w", format="NETCDF4")

            # define axis size
            ncout.createDimension("time", None)  # unlimited
            ncout.createDimension("latitude", nlat)
            ncout.createDimension("longitude", nlon)

            # create time axis
            time = ncout.createVariable("time", dtype("float64").char, ("time",))
            time.long_name = "time"
            time.units = tm_units
            time.calendar = "gregorian"
            time.axis = "T"

            # create latitude axis
            latitude = ncout.createVariable("latitude", dtype("float64").char, ("latitude"))
            latitude.standard_name = "latitude"
            latitude.long_name = "latitude"
            latitude.units = "degrees_north"
            latitude.axis = "Y"

            # create longitude axis
            longitude = ncout.createVariable("longitude", dtype("float64").char, ("longitude"))
            longitude.standard_name = "longitude"
            longitude.long_name = "longitude"
            longitude.units = "degrees_east"
            longitude.axis = "X"

            # create variable array
            upsi_out = ncout.createVariable("upsi", dtype("float64").char, ("time", "latitude", "longitude"))
            upsi_out.long_name = "Zonal Non-Divergent Wind Component at 700hPa"
            upsi_out.units = "m s**-1"

            # create variable array
            vpsi_out = ncout.createVariable("vpsi", dtype("float64").char, ("time", "latitude", "longitude"))
            vpsi_out.long_name = "Meridional Non-Divergent Wind Component at 700hPa"
            vpsi_out.units = "m s**-1"

            # copy axis from original dataset
            time[:] = tm[:]
            longitude[:] = lons[:]
            latitude[:] = lats[:]
            upsi_out[:, :, :] = upsi[:, :, :]
            vpsi_out[:, :, :] = vpsi[:, :, :]

            # print(ncout)
            ncout.close()
            # del upsi, vpsi, w, uwnd_info, vwnd_info

    non_divergent_wind(file_in)
