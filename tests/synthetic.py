"""Synthetic global wind fields for tests that must not depend on downloads.

Builds a single cyclonic vortex embedded in a uniform easterly flow and
translates it westward across the international dateline, which is exactly the
case the tracker used to lose.
"""

import numpy as np
import xarray as xr

KM_PER_DEG_LAT = 110.57
KM_PER_DEG_LON_EQ = 111.32


def _lon_delta(lon_b, lon_a):
    """Signed shortest separation, duplicated here so tests stay independent of qtrack.geo."""
    return (np.asarray(lon_b) - np.asarray(lon_a) + 180.0) % 360.0 - 180.0


def vortex_winds(lon2d, lat2d, lon_c, lat_c, vmax=15.0, rmax_km=500.0, u_background=-8.0):
    """Cyclonic (counterclockwise) vortex on a uniform easterly background.

    Tangential wind follows ``vmax * (r/rmax) * exp(0.5 * (1 - (r/rmax)**2))``,
    which peaks at ``rmax_km`` and decays smoothly, so the field is differentiable
    everywhere and the curvature vorticity has a clean single maximum.
    """
    dx_km = _lon_delta(lon2d, lon_c) * KM_PER_DEG_LON_EQ * np.cos(np.deg2rad(lat_c))
    dy_km = (lat2d - lat_c) * KM_PER_DEG_LAT
    r = np.hypot(dx_km, dy_km)

    scaled = r / rmax_km
    v_tan = vmax * scaled * np.exp(0.5 * (1.0 - scaled**2))

    with np.errstate(invalid="ignore", divide="ignore"):
        u = np.where(r > 0, -v_tan * dy_km / r, 0.0)
        v = np.where(r > 0, v_tan * dx_km / r, 0.0)

    return u + u_background, v


def westward_vortex_dataset(
    path,
    n_times=40,
    start_lon=-150.0,
    center_lat=12.0,
    speed_ms=7.0,
    temporal_res_h=6,
    lat_bounds=(30.0, -5.0),
    res=1.0,
    start="2010-08-01",
    **vortex_kwargs,
):
    """Write a global netCDF of 700 hPa winds with one vortex crossing the dateline.

    Longitude is a full global 1-degree axis; latitude runs north to south, which
    is what ``prep_data`` produces. The vortex starts at ``start_lon`` and moves
    due west at ``speed_ms``, so with the defaults it crosses 180 roughly
    two-thirds of the way through the run -- well after the three days a track
    needs before extrapolation is allowed to kick in.

    Returns the array of true centre longitudes (unwrapped, so monotonically
    decreasing) for comparison against the recovered track.
    """
    lon = np.arange(-180.0, 180.0, res)
    lat = np.arange(lat_bounds[0], lat_bounds[1] - 0.5 * res, -res)
    time = np.arange(n_times) * np.timedelta64(temporal_res_h, "h").astype("timedelta64[ns]") + np.datetime64(start, "ns")

    deg_per_step = speed_ms * temporal_res_h * 3600.0 / (KM_PER_DEG_LON_EQ * np.cos(np.deg2rad(center_lat)) * 1000.0)
    true_lon = start_lon - deg_per_step * np.arange(n_times)

    lon2d, lat2d = np.meshgrid(lon, lat)
    u = np.empty((n_times, lat.size, lon.size))
    v = np.empty_like(u)
    for i, lon_c in enumerate(true_lon):
        u[i], v[i] = vortex_winds(lon2d, lat2d, lon_c, center_lat, **vortex_kwargs)

    ds = xr.Dataset(
        {
            "u": (("time", "latitude", "longitude"), u),
            "v": (("time", "latitude", "longitude"), v),
        },
        coords={"time": time, "latitude": lat, "longitude": lon},
    )
    ds["u"].attrs["units"] = "m s**-1"
    ds["v"].attrs["units"] = "m s**-1"
    encoding = {"time": {"units": f"hours since {start} 00:00:00", "calendar": "gregorian", "dtype": "float64"}}
    ds.to_netcdf(path, encoding=encoding)
    return true_lon
