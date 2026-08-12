import numpy as np
import pytest

from qtrack import geo


class TestWrap:
    def test_wrap180_matches_prep_data_convention(self):
        # prep_data uses (lon + 180) % 360 - 180, i.e. [-180, 180)
        assert geo.wrap180(-180.0) == -180.0
        assert geo.wrap180(180.0) == -180.0
        assert geo.wrap180(181.0) == -179.0
        assert geo.wrap180(-181.0) == 179.0
        assert geo.wrap180(0.0) == 0.0

    def test_wrap180_array(self):
        out = geo.wrap180(np.array([-190.0, 10.0, 200.0]))
        np.testing.assert_allclose(out, [170.0, 10.0, -160.0])

    def test_wrap360(self):
        assert geo.wrap360(-1.0) == 359.0
        assert geo.wrap360(360.0) == 0.0


class TestLonDelta:
    def test_across_the_dateline(self):
        assert geo.lon_delta(-179.0, 179.0) == pytest.approx(2.0)
        assert geo.lon_delta(179.0, -179.0) == pytest.approx(-2.0)

    def test_ordinary_separation(self):
        assert geo.lon_delta(-15.0, -20.0) == pytest.approx(5.0)
        assert geo.lon_delta(-25.0, -20.0) == pytest.approx(-5.0)

    def test_nan_propagates(self):
        assert np.isnan(geo.lon_delta(np.nan, 10.0))

    def test_array_against_scalar(self):
        out = geo.lon_delta(np.array([179.0, -179.0]), -180.0)
        np.testing.assert_allclose(out, [-1.0, 1.0])

    def test_returns_scalar_for_scalar_input(self):
        assert isinstance(geo.lon_delta(1.0, 2.0), float)


class TestIsWestward:
    def test_westward_across_the_seam(self):
        # Westward is decreasing longitude, so a wave at 179.5W (-179.5)
        # continuing west enters the West Pacific at 179.5E (+179.5).
        assert geo.is_westward(179.5, -179.5) is True
        # The old raw test got this backwards: 179.5 <= -179.5 is False.
        assert (179.5 <= -179.5) is False

    def test_eastward_across_the_seam(self):
        assert geo.is_westward(-179.5, 179.5) is False

    def test_matches_old_behaviour_away_from_the_seam(self):
        for new, old in [(-30.0, -20.0), (-10.0, -20.0), (-20.0, -20.0)]:
            assert geo.is_westward(new, old) == (new <= old)

    def test_nan_is_not_westward(self):
        assert geo.is_westward(np.nan, -20.0) is False


class TestLonInRange:
    def test_window_straddling_the_dateline(self):
        assert geo.lon_in_range(175.0, 170.0, -170.0) is True
        assert geo.lon_in_range(-175.0, 170.0, -170.0) is True
        assert geo.lon_in_range(0.0, 170.0, -170.0) is False

    def test_window_not_straddling_the_dateline(self):
        assert geo.lon_in_range(0.0, -170.0, 170.0) is True
        assert geo.lon_in_range(175.0, -170.0, 170.0) is False

    def test_endpoints(self):
        assert geo.lon_in_range(-17.0, -17.0, 40.0, inclusive="both") is True
        assert geo.lon_in_range(-17.0, -17.0, 40.0, inclusive="east") is False
        assert geo.lon_in_range(40.0, -17.0, 40.0, inclusive="east") is True
        assert geo.lon_in_range(40.0, -17.0, 40.0, inclusive="west") is False

    def test_equal_bounds_is_the_whole_globe(self):
        for lon in (-180.0, -60.0, 0.0, 179.0):
            assert geo.lon_in_range(lon, -180.0, 180.0) is True

    def test_nan_is_never_inside(self):
        assert geo.lon_in_range(np.nan, -180.0, 180.0) is False

    def test_reproduces_legacy_left_right_bounds(self):
        # Old code: lon > -180 and lon < 40 (strict on both ends).
        lons = np.arange(-179.5, 179.5, 1.0)
        expected = (lons > -180.0) & (lons < 40.0)
        got = geo.lon_in_range(lons, -180.0, 40.0, inclusive="neither")
        np.testing.assert_array_equal(got, expected)


class TestLonInAnyRange:
    def test_none_means_unrestricted(self):
        assert geo.lon_in_any_range(123.0, None) is True

    def test_single_tuple_accepted(self):
        assert geo.lon_in_any_range(0.0, (-35.0, 40.0)) is True

    def test_multiple_windows(self):
        ranges = [(-35.0, 40.0), (160.0, -160.0)]
        assert geo.lon_in_any_range(170.0, ranges) is True
        assert geo.lon_in_any_range(-170.0, ranges) is True
        assert geo.lon_in_any_range(-100.0, ranges) is False


class TestUnwrap:
    def test_round_trip_through_wrap180(self):
        raw = np.array([178.0, 179.0, -180.0, -179.0, -178.0])
        np.testing.assert_allclose(geo.wrap180(geo.unwrap_lon(raw)), raw)

    def test_westward_crossing_is_monotonic(self):
        raw = np.array([-178.0, -179.0, 180.0, 179.0, 178.0])
        out = geo.unwrap_lon(geo.wrap180(raw))
        assert np.all(np.diff(out) < 0)
        np.testing.assert_allclose(np.diff(out), -1.0)

    def test_nan_gaps_are_preserved_and_bridged(self):
        raw = np.array([179.0, np.nan, -179.0, -178.0])
        out = geo.unwrap_lon(raw)
        assert np.isnan(out[1])
        # The two-degree jump across the gap must not become -358.
        assert out[2] - out[0] == pytest.approx(2.0)

    def test_all_nan(self):
        out = geo.unwrap_lon(np.full(4, np.nan))
        assert np.isnan(out).all()

    def test_single_value(self):
        out = geo.unwrap_lon(np.array([np.nan, 12.0, np.nan]))
        assert out[1] == 12.0

    def test_rejects_2d(self):
        with pytest.raises(ValueError):
            geo.unwrap_lon(np.zeros((2, 3)))


class TestOnUnwrapped:
    def test_interpolation_across_the_seam(self):
        # A gap in the middle of a seam crossing must interpolate the short way.
        series = np.array([179.0, np.nan, -179.0])
        filled = geo.on_unwrapped(lambda x: np.interp(np.arange(3), [0, 2], [x[0], x[2]]), series)
        assert filled[1] == pytest.approx(-180.0)


class TestGridDetection:
    def test_global_one_degree(self):
        lon = np.arange(-180.0, 180.0, 1.0)
        assert geo.is_global_lon(lon) is True
        assert geo.grid_resolution(lon) == pytest.approx(1.0)

    def test_global_half_degree(self):
        lon = np.arange(0.0, 360.0, 0.5)
        assert geo.is_global_lon(geo.wrap180(lon)) is True

    def test_regional_is_not_global(self):
        lon = np.arange(-120.0, 60.0, 1.0)
        assert geo.is_global_lon(lon) is False

    def test_degenerate_axis(self):
        assert geo.is_global_lon(np.array([1.0])) is False


class TestFindNearestLon:
    def test_wraps_instead_of_clamping(self):
        lon = np.arange(-180.0, 180.0, 1.0)
        idx, val = geo.find_nearest_lon(lon, -181.0)
        assert val == 179.0
        assert idx == len(lon) - 1

    def test_ordinary_lookup(self):
        lon = np.arange(-180.0, 180.0, 1.0)
        idx, val = geo.find_nearest_lon(lon, -17.4)
        assert val == -17.0
        assert lon[idx] == -17.0


class TestPadding:
    def test_pad_and_unpad_round_trip(self):
        arr = np.arange(12.0).reshape(3, 4)
        padded = geo.pad_lon(arr, 2)
        assert padded.shape == (3, 8)
        np.testing.assert_allclose(geo.unpad_lon(padded, 2), arr)

    def test_pad_wraps_the_longitude_axis(self):
        arr = np.array([[0.0, 1.0, 2.0, 3.0]])
        padded = geo.pad_lon(arr, 1)
        np.testing.assert_allclose(padded, [[3.0, 0.0, 1.0, 2.0, 3.0, 0.0]])

    def test_zero_pad_is_a_no_op(self):
        arr = np.arange(4.0)
        np.testing.assert_allclose(geo.pad_lon(arr, 0), arr)

    def test_pad_lon_coord_stays_monotonic(self):
        lon = np.arange(-180.0, 180.0, 1.0)
        padded = geo.pad_lon_coord(lon, 5)
        assert np.all(np.diff(padded) > 0)
        np.testing.assert_allclose(np.diff(padded), 1.0)
        assert padded[0] == -185.0
        assert padded[-1] == 184.0

    def test_pad_lon_coord_keeps_gradient_uniform(self):
        # The whole point: dx must not spike at the seam.
        lon = np.arange(-180.0, 180.0, 1.0)
        naive = np.concatenate((lon[-5:], lon, lon[:5]))
        assert np.abs(np.gradient(naive)).max() > 100.0
        padded = geo.pad_lon_coord(lon, 5)
        np.testing.assert_allclose(np.gradient(padded), 1.0)


class TestWrapIndexDelta:
    def test_across_the_seam(self):
        assert geo.wrap_index_delta(1 - 359, 360) == 2
        assert geo.wrap_index_delta(359 - 1, 360) == -2

    def test_ordinary(self):
        assert geo.wrap_index_delta(10 - 8, 360) == 2
