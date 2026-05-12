from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from fullpos import (
    average_interpolate,
    bilinear_interpolate,
    horizontal_interpolate,
    nearest_interpolate,
    quadratic12_interpolate,
    vertical_capabilities,
    vertical_interpolate,
)
from fullpos.errors import FullposNotImplementedError
from fullpos.grids import gaussian_latitudes, octahedral_pl, regular_longitudes


def test_bilinear_interpolate_matches_source_grid_points() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    obj = xr.DataArray(
        np.arange(8 * 16, dtype=np.float64).reshape(8, 16),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    out = bilinear_interpolate(obj, target_lats=target_lats, target_lons=target_lons)

    np.testing.assert_allclose(out.values, obj.values[2:-2], atol=1.0e-12)


def test_quadratic12_interpolate_matches_source_grid_points() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    obj = xr.DataArray(
        np.arange(8 * 16, dtype=np.float64).reshape(8, 16),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    out = quadratic12_interpolate(obj, target_lats=target_lats, target_lons=target_lons)

    np.testing.assert_allclose(out.values, obj.values[2:-2], atol=1.0e-11)


def test_nearest_interpolate_matches_source_grid_points() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    obj = xr.DataArray(
        np.arange(8 * 16, dtype=np.float64).reshape(8, 16),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    out = nearest_interpolate(obj, target_lats=target_lats, target_lons=target_lons)

    np.testing.assert_allclose(out.values, obj.values[2:-2], atol=1.0e-12)


def test_average_interpolate_matches_fullpos_2x2_halo_mean_at_grid_points() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    values = np.arange(8 * 16, dtype=np.float64).reshape(8, 16)
    obj = xr.DataArray(values, dims=("latitude", "longitude"), coords={"latitude": lats, "longitude": lons})
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    out = average_interpolate(obj, target_lats=target_lats, target_lons=target_lons, average_radius=1)

    expected = np.empty_like(values[2:-2])
    for row in range(2, 6):
        for col in range(values.shape[1]):
            expected[row - 2, col] = np.mean(
                [
                    values[row, col],
                    values[row, (col + 1) % values.shape[1]],
                    values[row + 1, col],
                    values[row + 1, (col + 1) % values.shape[1]],
                ]
            )
    np.testing.assert_allclose(out.values, expected, atol=1.0e-12)


@pytest.mark.parametrize("method", ["bilinear", "quadratic12"])
def test_horizontal_interpolate_supports_packed_reduced_o_grid_points(method) -> None:
    pl = octahedral_pl(4)
    lats = gaussian_latitudes(pl.size)
    rows = [row * 1000.0 + np.arange(row_nlon, dtype=np.float64) for row, row_nlon in enumerate(pl)]
    values = np.concatenate(rows)
    target_lats = []
    target_lons = []
    expected = []
    for row in range(2, pl.size - 2):
        row_lons = regular_longitudes(int(pl[row]))
        target_lats.append(np.full(row_lons.shape, lats[row]))
        target_lons.append(row_lons)
        expected.append(rows[row])

    out = horizontal_interpolate(
        values,
        source_grid="O4",
        target_lats=np.concatenate(target_lats),
        target_lons=np.concatenate(target_lons),
        method=method,
    )

    np.testing.assert_allclose(out, np.concatenate(expected), atol=1.0e-10)


def test_horizontal_interpolate_packed_reduced_data_array_uses_grib_pl() -> None:
    pl = octahedral_pl(4)
    lats = gaussian_latitudes(pl.size)
    rows = [row * 1000.0 + np.arange(row_nlon, dtype=np.float64) for row, row_nlon in enumerate(pl)]
    values = np.stack([np.concatenate(rows), np.concatenate(rows) + 100.0])
    obj = xr.DataArray(
        values,
        dims=("time", "values"),
        coords={"time": [0, 1]},
        attrs={"GRIB_pl": pl},
        name="t",
    )
    target_lats = []
    target_lons = []
    expected = []
    for row in range(2, pl.size - 2):
        row_lons = regular_longitudes(int(pl[row]))
        target_lats.append(np.full(row_lons.shape, lats[row]))
        target_lons.append(row_lons)
        expected.append(rows[row])

    out = horizontal_interpolate(
        obj,
        target_lats=np.concatenate(target_lats),
        target_lons=np.concatenate(target_lons),
        method="bilinear",
        chunks={"time": 1},
    )

    assert out.dims == ("time", "points")
    np.testing.assert_allclose(out.isel(time=0).values, np.concatenate(expected), atol=1.0e-10)
    np.testing.assert_allclose(out.isel(time=1).values, np.concatenate(expected) + 100.0, atol=1.0e-10)


def test_horizontal_interpolate_data_array_supports_chunks() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    base = np.arange(8 * 16, dtype=np.float64).reshape(8, 16)
    obj = xr.DataArray(
        np.stack([base, base + 100.0]),
        dims=("time", "latitude", "longitude"),
        coords={"time": [0, 1], "latitude": lats, "longitude": lons},
        name="t",
        attrs={"units": "K"},
    )
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    out = horizontal_interpolate(
        obj,
        target_lats=target_lats,
        target_lons=target_lons,
        method="bilinear",
        chunks={"time": 1},
    )

    assert out.dims == ("time", "target_y", "target_x")
    assert out.name == "t"
    assert out.attrs["units"] == "K"
    np.testing.assert_allclose(out.isel(time=0).values, base[2:-2], atol=1.0e-12)
    np.testing.assert_allclose(out.isel(time=1).values, base[2:-2] + 100.0, atol=1.0e-12)


def test_horizontal_interpolate_dataset_skips_non_horizontal_variables() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    base = np.arange(8 * 16, dtype=np.float64).reshape(8, 16)
    ds = xr.Dataset(
        data_vars={
            "t": (("latitude", "longitude"), base),
            "surface": ("time", np.array([1.0, 2.0])),
        },
        coords={"latitude": lats, "longitude": lons, "time": [0, 1]},
    )
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    out = horizontal_interpolate(
        ds,
        target_lats=target_lats,
        target_lons=target_lons,
        method="bilinear",
    )

    assert out["t"].dims == ("target_y", "target_x")
    assert out["surface"].dims == ("time",)
    np.testing.assert_allclose(out["t"].values, base[2:-2], atol=1.0e-12)


def test_masked_horizontal_interpolation_is_still_unimplemented() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    obj = xr.DataArray(
        np.ones((8, 16), dtype=np.float64),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )

    with pytest.raises(FullposNotImplementedError):
        horizontal_interpolate(
            obj,
            target_lats=np.array([lats[3]]),
            target_lons=np.array([lons[0]]),
            source_mask=np.ones_like(obj.values, dtype=bool),
        )


def test_horizontal_interpolate_validates_bad_xarray_chunks() -> None:
    obj = xr.DataArray(
        np.ones((8, 16)),
        dims=("latitude", "longitude"),
        coords={"latitude": gaussian_latitudes(8), "longitude": regular_longitudes(16)},
    )

    with pytest.raises(ValueError, match="not present"):
        horizontal_interpolate(
            obj,
            target_lats=np.array([0.0]),
            target_lons=np.array([0.0]),
            chunks={"level": 1},
        )


def test_horizontal_interpolate_validates_dataset_variables() -> None:
    ds = xr.Dataset(
        data_vars={"surface": ("time", np.array([1.0, 2.0]))},
        coords={"time": [0, 1]},
    )

    with pytest.raises(ValueError, match="surface"):
        horizontal_interpolate(
            ds,
            target_lats=np.array([0.0]),
            target_lons=np.array([0.0]),
            variables=["surface"],
        )


def test_horizontal_interpolate_rejects_backend_argument() -> None:
    obj = xr.DataArray(
        np.ones((8, 16)),
        dims=("latitude", "longitude"),
        coords={"latitude": gaussian_latitudes(8), "longitude": regular_longitudes(16)},
    )

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        horizontal_interpolate(
            obj,
            target_lats=np.array([0.0]),
            target_lons=np.array([0.0]),
            backend="python",
        )


def test_vertical_interpolation_reports_native_pressure_target() -> None:
    capabilities = vertical_capabilities()

    assert capabilities["pressure"] == "native"
    with pytest.raises(TypeError, match="requires values input"):
        vertical_interpolate(target="pressure")
