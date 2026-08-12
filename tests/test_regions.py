"""The region table must be a faithful refactor of the old scalar longitude tests.

``legacy_regions()`` encodes exactly what ``run_tracking`` used to do with
``lon <= -17``, ``lon <= extra_lon`` and ``lon <= extra_transition``. These tests
sweep the historical tracking domain (-180 to 40) and assert the table gives the
same answer at every gridpoint, so the switch to a region lookup can be reviewed
independently of the switch to global defaults.
"""

import numpy as np
import pytest

from qtrack import regions as R

# The historical scalar settings, from run_tracking's settings block.
AFRICA_COAST = -17.0
EXTRA_LON = -20.0
EXTRA_TRANSITION = -60.0
BACK_CUTOFF_LAND = 100.0
BACK_CUTOFF_OCEAN = 300.0
EXTRAP_DIST = 700.0
EXTRAP_DIST_CARIB = 500.0
CONNECT_STEP = 2

# Every gridpoint the old tracker could reach, at 1-degree centres.
DOMAIN = np.arange(-179.5, 40.0, 1.0)


@pytest.fixture
def table():
    return R.legacy_regions()


def test_every_domain_longitude_resolves(table):
    for lon in DOMAIN:
        assert R.region_for(lon, table) is not None


@pytest.mark.parametrize("lon", DOMAIN)
def test_land_ocean_back_cutoff_matches_legacy(table, lon):
    # Old: if AEW_lon <= -17: ocean cutoffs, else land cutoffs
    expected = BACK_CUTOFF_OCEAN if lon <= AFRICA_COAST else BACK_CUTOFF_LAND
    assert R.region_for(lon, table).back_cutoff == expected


@pytest.mark.parametrize("lon", DOMAIN)
def test_over_land_flag_matches_legacy(table, lon):
    assert R.region_for(lon, table).over_land == (lon > AFRICA_COAST)


@pytest.mark.parametrize("lon", DOMAIN)
def test_extrapolation_gate_matches_legacy(table, lon):
    # Old: extrapolation runs only where AEW_lon <= extra_lon
    assert R.region_for(lon, table).extrapolate == (lon <= EXTRA_LON)


@pytest.mark.parametrize("lon", DOMAIN)
def test_speed_limit_gate_matches_legacy(table, lon):
    # Old: speed limit applies only where AEW_lon_slc <= extra_lon
    assert R.region_for(lon, table).apply_speed_limit == (lon <= EXTRA_LON)


@pytest.mark.parametrize("lon", DOMAIN)
def test_caribbean_extra_dist_matches_legacy(table, lon):
    # Old: if AEW_lon <= extra_transition: extra_distance_60 else extra_distance
    expected = EXTRAP_DIST_CARIB if lon <= EXTRA_TRANSITION else EXTRAP_DIST
    assert R.region_for(lon, table).extra_dist == expected


@pytest.mark.parametrize("lon", DOMAIN)
def test_reconnect_window_matches_legacy(table, lon):
    # Old: elif lon_first >= -17 and tm_step > connect_step: skip
    expected = CONNECT_STEP if lon >= AFRICA_COAST else 4
    # The two tests disagree only at exactly -17, which no 1-degree centre hits.
    if lon != AFRICA_COAST:
        assert R.region_for(lon, table).reconnect_tm_steps == expected


def test_scalar_kwargs_move_the_boundaries():
    table = R.legacy_regions(extrap_longitude_start=-30.0, carib_longitude_start=-80.0)
    assert R.region_for(-25.0, table).extrapolate is False
    assert R.region_for(-35.0, table).extrapolate is True
    assert R.region_for(-70.0, table).extra_dist == 700.0
    assert R.region_for(-90.0, table).extra_dist == 500.0


class TestSeamCrossingWindows:
    def test_legacy_caribbean_window_wraps_the_dateline(self):
        table = R.legacy_regions()
        assert R.region_for(-179.0, table).name == "caribbean"
        assert R.region_for(-61.0, table).name == "caribbean"
        assert R.region_for(-59.0, table).name == "atlantic"

    def test_user_window_may_straddle_the_dateline(self):
        west_pacific = R.Region(name="west_pacific", lon_range=(150.0, -150.0), extra_dist=900.0)
        table = [west_pacific] + R.global_regions()
        assert R.region_for(170.0, table).name == "west_pacific"
        assert R.region_for(-170.0, table).name == "west_pacific"
        assert R.region_for(120.0, table).name == "global_ocean"


class TestGlobalTable:
    def test_atlantic_sector_is_unchanged_from_legacy(self):
        legacy = R.legacy_regions()
        glob = R.global_regions()
        # Africa through the Caribbean must behave identically.
        for lon in np.arange(-99.5, 40.0, 1.0):
            a = R.region_for(lon, legacy)
            b = R.region_for(lon, glob)
            assert a.name == b.name, lon
            assert a.extrapolate == b.extrapolate
            assert a.extra_dist == b.extra_dist
            assert a.back_cutoff == b.back_cutoff
            assert a.lock_extrapolation == b.lock_extrapolation

    def test_pacific_gets_the_general_ocean_settings(self):
        glob = R.global_regions()
        for lon in (-150.0, -179.0, 179.0, 90.0, 150.0):
            region = R.region_for(lon, glob)
            assert region.name == "global_ocean"
            assert region.extrapolate is True
            assert region.extra_dist == EXTRAP_DIST
            # The Caribbean lock would stop waves being re-acquired after
            # weakening, which global tracking needs to do.
            assert region.lock_extrapolation is False

    def test_africa_still_special(self):
        glob = R.global_regions()
        africa = R.region_for(0.0, glob)
        assert africa.over_land is True
        assert africa.extrapolate is False
        assert africa.back_cutoff == BACK_CUTOFF_LAND
        assert africa.reconnect_tm_steps == CONNECT_STEP


class TestRegionFor:
    def test_nan_falls_back(self):
        table = R.global_regions()
        assert R.region_for(np.nan, table).name == "global_ocean"

    def test_table_without_fallback_is_rejected(self):
        with pytest.raises(ValueError):
            R.region_for(0.0, [R.Region(name="x", lon_range=(10.0, 20.0))])


class TestOverLandAnywhere:
    def test_matches_legacy_nanmax_test(self):
        table = R.legacy_regions()
        over_africa = np.array([-30.0, -20.0, -10.0, np.nan])
        not_over_africa = np.array([-80.0, -70.0, np.nan, -60.0])
        assert R.over_land_anywhere(over_africa, table) is True
        assert R.over_land_anywhere(not_over_africa, table) is False
        # ... which is what np.nanmax(lon) > -17 used to say
        assert bool(np.nanmax(over_africa) > AFRICA_COAST) is True
        assert bool(np.nanmax(not_over_africa) > AFRICA_COAST) is False

    def test_all_nan_track(self):
        assert R.over_land_anywhere(np.full(3, np.nan), R.global_regions()) is False


class TestResolveRegions:
    def test_default_is_global(self):
        assert R.resolve_regions(None)[-1].name == "global_ocean"

    def test_preset_names(self):
        assert R.resolve_regions("global")[-1].name == "global_ocean"
        assert R.resolve_regions("atlantic")[-1].name == "open_ocean"
        assert R.resolve_regions("legacy")[-1].name == "open_ocean"

    def test_unknown_preset(self):
        with pytest.raises(ValueError, match="Unknown region preset"):
            R.resolve_regions("pacific")

    def test_explicit_list_passes_through(self):
        table = R.global_regions()
        assert R.resolve_regions(table) == table

    def test_list_needs_a_fallback(self):
        with pytest.raises(ValueError):
            R.resolve_regions([R.Region(name="x", lon_range=(10.0, 20.0))])

    def test_rejects_non_regions(self):
        with pytest.raises(TypeError):
            R.resolve_regions([{"name": "x"}])

    def test_kwargs_reach_the_builder(self):
        table = R.resolve_regions("global", extrap_dist=1234.0)
        assert R.region_for(-40.0, table).extra_dist == 1234.0
