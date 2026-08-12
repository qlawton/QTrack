"""
Basin-dependent tracker settings, expressed as longitude windows.

The original tracker encoded "where am I, and how should I behave here?" as three
bare longitude comparisons scattered through ``run_tracking`` and
``run_postprocessing``:

* ``lon <= -17``  -- over the ocean rather than over Africa, so use the wider
  backward-search cutoffs and the full track-reconnection window;
* ``lon <= -20``  -- far enough offshore to switch from the Elless-and-Torn
  banding method to linear extrapolation, and to apply the forward speed limit;
* ``lon <= -60``  -- into the Caribbean, so tighten the extrapolation distance and
  stop falling back to the banding method.

Those tests silently inverted west of the dateline: a wave at 150E fails
``lon <= -20`` and so never extrapolates, while a wave at 150W passes
``lon <= -60`` and picks up Caribbean settings in the middle of the Pacific.

This module replaces them with a table of :class:`Region` records, each covering a
longitude window and carrying the settings that used to be gated on the scalar.
Windows are modular, so they may cross the dateline, and the table is ordered:
:func:`region_for` returns the first match, falling back to the entry with no
window at all.
"""

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

from .geo import lon_in_range

__all__ = [
    "Region",
    "region_for",
    "resolve_regions",
    "legacy_regions",
    "global_regions",
    "over_land_anywhere",
]

# Longitude of the West African coast; the historical land/ocean divide.
AFRICA_COAST_LON = -17.0
# East edge of the African land window (the historical initiation east bound).
AFRICA_EAST_LON = 40.0


@dataclass
class Region:
    """Tracker settings that apply within one longitude window.

    ``lon_range`` is a ``(west, east)`` pair traversed *eastward* from ``west``,
    so ``(170, -170)`` is the 20-degree window straddling the dateline. It is
    half-open -- exclusive at ``west``, inclusive at ``east`` -- which reproduces
    the ``lon <= threshold`` / ``lon > threshold`` semantics of the code it
    replaces. ``lon_range=None`` marks the fallback region, matching anywhere.
    """

    name: str
    lon_range: Optional[Tuple[float, float]] = None

    # -- Elless-and-Torn banding / continuation settings ------------------
    over_land: bool = False
    back_cutoff: float = 300.0  # km, "backward" (eastward) match limit, 1 step back
    back_cutoff_long: float = 300.0  # km, same but 2 steps back
    lat_jump_limit: float = 3.0  # deg latitude a wave may shift per timestep

    # -- Extrapolation settings -------------------------------------------
    extrapolate: bool = True
    extra_dist: float = 700.0  # km an extrapolated centre may sit from the last point
    extra_it: int = 5  # timesteps used for the linear fit
    centroid_rad_extra: float = 600.0  # km, centroid radius for extrapolated guesses
    lock_extrapolation: bool = False  # once extrapolated here, never fall back to banding

    # -- Other behaviour ---------------------------------------------------
    apply_speed_limit: bool = True
    allow_initiation: bool = True
    reconnect_tm_steps: int = 4  # postprocessing reconnection search window

    def contains(self, lon) -> bool:
        """True when ``lon`` falls in this region's window (always True for the fallback)."""
        if self.lon_range is None:
            return True
        return bool(lon_in_range(lon, self.lon_range[0], self.lon_range[1], inclusive="east"))


def region_for(lon, regions: Sequence[Region]) -> Region:
    """First region in ``regions`` whose window contains ``lon``.

    A non-finite ``lon`` (a track gap) resolves to the fallback region so callers
    never have to special-case NaN.
    """
    if lon != lon:  # NaN
        return _fallback(regions)
    for region in regions:
        if region.lon_range is not None and region.contains(lon):
            return region
    return _fallback(regions)


def _fallback(regions: Sequence[Region]) -> Region:
    for region in regions:
        if region.lon_range is None:
            return region
    raise ValueError("region table has no fallback entry (a Region with lon_range=None)")


def over_land_anywhere(lon_series, regions: Sequence[Region]) -> bool:
    """True when any point of a track falls in a region flagged ``over_land``.

    Replaces ``np.nanmax(AEW_lon[row, :]) > -17``, which is meaningless once
    longitudes wrap. Within the historical -180...40 domain the two agree.
    """
    for lon in lon_series:
        if lon != lon:
            continue
        if region_for(lon, regions).over_land:
            return True
    return False


def _land_region(name, lon_range, back_cutoff_land, back_cutoff_long_land, land_lat_limit, connect_step):
    """Africa-style settings: banding method only, tight backward search."""
    return Region(
        name=name,
        lon_range=lon_range,
        over_land=True,
        back_cutoff=back_cutoff_land,
        back_cutoff_long=back_cutoff_long_land,
        lat_jump_limit=land_lat_limit,
        extrapolate=False,
        apply_speed_limit=False,
        reconnect_tm_steps=connect_step,
    )


def _ocean_region(
    name,
    lon_range,
    back_cutoff_ocean,
    back_cutoff_long_ocean,
    land_lat_limit,
    extrapolate,
    extra_dist,
    extra_it,
    centroid_rad_extra,
    lock_extrapolation=False,
    reconnect_tm_steps=4,
):
    return Region(
        name=name,
        lon_range=lon_range,
        over_land=False,
        back_cutoff=back_cutoff_ocean,
        back_cutoff_long=back_cutoff_long_ocean,
        lat_jump_limit=land_lat_limit,
        extrapolate=extrapolate,
        extra_dist=extra_dist,
        extra_it=extra_it,
        centroid_rad_extra=centroid_rad_extra,
        lock_extrapolation=lock_extrapolation,
        apply_speed_limit=extrapolate,
        reconnect_tm_steps=reconnect_tm_steps,
    )


def legacy_regions(
    africa_coast_lon=AFRICA_COAST_LON,
    africa_east_lon=AFRICA_EAST_LON,
    extrap_longitude_start=-20.0,
    carib_longitude_start=-60.0,
    extrap_dist=700.0,
    extrap_dist_carib=500.0,
    extra_it_20=5,
    extra_it_60=5,
    centroid_radius=600.0,
    centroid_radius_carib=None,
    back_cutoff_land=100.0,
    back_cutoff_ocean=300.0,
    back_cutoff_long_land=100.0,
    back_cutoff_long_ocean=300.0,
    land_lat_limit=3.0,
    connect_step=2,
):
    """Region table reproducing the pre-``global-adjusted`` Atlantic behaviour.

    Faithful within the historical tracking domain (-180 to 40). East of 40E the
    old code applied the "land" branch of ``lon <= -17``; here that longitude band
    is unreachable, because it is covered by the Caribbean window that wraps
    across the dateline. ``tests/test_regions.py`` pins the equivalence over
    -180...40.
    """
    if centroid_radius_carib is None:
        centroid_radius_carib = centroid_radius
    return [
        _land_region("africa_land", (africa_coast_lon, africa_east_lon), back_cutoff_land, back_cutoff_long_land, land_lat_limit, connect_step),
        _ocean_region(
            "east_atlantic",
            (extrap_longitude_start, africa_coast_lon),
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=False,
            extra_dist=extrap_dist,
            extra_it=extra_it_20,
            centroid_rad_extra=centroid_radius,
        ),
        _ocean_region(
            "atlantic",
            (carib_longitude_start, extrap_longitude_start),
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=True,
            extra_dist=extrap_dist,
            extra_it=extra_it_20,
            centroid_rad_extra=centroid_radius,
        ),
        _ocean_region(
            "caribbean",
            (africa_east_lon, carib_longitude_start),  # wraps eastward across the dateline
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=True,
            extra_dist=extrap_dist_carib,
            extra_it=extra_it_60,
            centroid_rad_extra=centroid_radius_carib,
            lock_extrapolation=True,
        ),
        # Unreachable given the windows above, but a table must have a fallback.
        _ocean_region(
            "open_ocean",
            None,
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=True,
            extra_dist=extrap_dist,
            extra_it=extra_it_20,
            centroid_rad_extra=centroid_radius,
        ),
    ]


def global_regions(
    africa_coast_lon=AFRICA_COAST_LON,
    africa_east_lon=AFRICA_EAST_LON,
    extrap_longitude_start=-20.0,
    carib_longitude_start=-60.0,
    carib_west_lon=-100.0,
    extrap_dist=700.0,
    extrap_dist_carib=500.0,
    extra_it_20=5,
    extra_it_60=5,
    centroid_radius=600.0,
    centroid_radius_carib=None,
    back_cutoff_land=100.0,
    back_cutoff_ocean=300.0,
    back_cutoff_long_land=100.0,
    back_cutoff_long_ocean=300.0,
    land_lat_limit=3.0,
    connect_step=2,
):
    """Default table: the Atlantic sector unchanged, general ocean settings elsewhere.

    Africa through the Caribbean keeps exactly the windows and settings the
    tracker has always used, so results in that sector are unaffected. The
    Caribbean window is now bounded on the west at ``carib_west_lon`` (default
    100W) instead of running to the array edge; everything west of that -- the
    East Pacific, West Pacific and Indian Ocean -- falls through to a general
    ocean region carrying the Atlantic settings.

    That is deliberate: the Caribbean regime is the *restrictive* one (a shorter
    extrapolation leash, and ``lock_extrapolation`` permanently disabling the
    banding fallback). Applying it around the whole globe would stop waves being
    re-acquired after they weaken and re-intensify, which is exactly what
    continuous global tracking needs to do.
    """
    if centroid_radius_carib is None:
        centroid_radius_carib = centroid_radius
    return [
        _land_region("africa_land", (africa_coast_lon, africa_east_lon), back_cutoff_land, back_cutoff_long_land, land_lat_limit, connect_step),
        _ocean_region(
            "east_atlantic",
            (extrap_longitude_start, africa_coast_lon),
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=False,
            extra_dist=extrap_dist,
            extra_it=extra_it_20,
            centroid_rad_extra=centroid_radius,
        ),
        _ocean_region(
            "atlantic",
            (carib_longitude_start, extrap_longitude_start),
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=True,
            extra_dist=extrap_dist,
            extra_it=extra_it_20,
            centroid_rad_extra=centroid_radius,
        ),
        _ocean_region(
            "caribbean",
            (carib_west_lon, carib_longitude_start),
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=True,
            extra_dist=extrap_dist_carib,
            extra_it=extra_it_60,
            centroid_rad_extra=centroid_radius_carib,
            lock_extrapolation=True,
        ),
        _ocean_region(
            "global_ocean",
            None,
            back_cutoff_ocean,
            back_cutoff_long_ocean,
            land_lat_limit,
            extrapolate=True,
            extra_dist=extrap_dist,
            extra_it=extra_it_20,
            centroid_rad_extra=centroid_radius,
        ),
    ]


_PRESETS = {
    "global": global_regions,
    "atlantic": legacy_regions,
    "legacy": legacy_regions,
}


def resolve_regions(regions, **kwargs):
    """Turn the ``regions`` argument of ``run_tracking`` into a list of :class:`Region`.

    Accepts ``None`` (the global default), a preset name (``"global"``,
    ``"atlantic"``/``"legacy"``), or an explicit sequence of :class:`Region`
    records. ``kwargs`` are the tracker's scalar settings, used to build a preset.
    """
    if regions is None:
        return global_regions(**kwargs)
    if isinstance(regions, str):
        try:
            builder = _PRESETS[regions.lower()]
        except KeyError:
            raise ValueError(f"Unknown region preset {regions!r}. Choose from {sorted(_PRESETS)} or pass a list of Region objects.") from None
        return builder(**kwargs)
    regions = list(regions)
    if not regions:
        raise ValueError("regions must not be empty")
    for region in regions:
        if not isinstance(region, Region):
            raise TypeError(f"regions must contain Region objects, got {type(region).__name__}")
    _fallback(regions)  # raises if the caller forgot a fallback entry
    return regions
