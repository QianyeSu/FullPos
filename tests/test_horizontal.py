from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from fullpos import (
    average_interpolate,
    bilinear_interpolate,
    horizontal_interpolate,
    land_sea_mask_to_grid,
    masked_surface_interpolate,
    nearest_interpolate,
    quadratic12_interpolate,
    vertical_capabilities,
    vertical_interpolate,
)
from fullpos.grids import gaussian_latitudes, octahedral_pl, parse_grid, regular_longitudes
from fullpos.native import add_native_runtime_dir


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


def test_quadratic12_shape_preserving_clamps_native_overshoot() -> None:
    lats = gaussian_latitudes(16)
    lons = regular_longitudes(32)
    values = np.zeros((16, 32), dtype=np.float64)
    values[7, 10] = 1.0
    obj = xr.DataArray(
        values,
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )
    target_lats = np.linspace(lats[5], lats[10], 21)
    target_lons = np.linspace(lons[8], lons[13], 21)
    target_lats, target_lons = np.meshgrid(target_lats, target_lons, indexing="ij")

    unconstrained = horizontal_interpolate(
        obj,
        target_lats=target_lats,
        target_lons=target_lons,
        method="quadratic12",
    )
    monotonic = horizontal_interpolate(
        obj,
        target_lats=target_lats,
        target_lons=target_lons,
        method="quadratic12",
        shape_preserving=True,
    )

    assert float(unconstrained.min()) < 0.0
    assert float(monotonic.min()) >= 0.0
    assert float(monotonic.max()) <= 1.0


def test_quadratic12_monotonic_method_alias_matches_shape_preserving_flag() -> None:
    lats = gaussian_latitudes(16)
    lons = regular_longitudes(32)
    values = np.zeros((16, 32), dtype=np.float64)
    values[7, 10] = 1.0
    obj = xr.DataArray(
        values,
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )
    target_lats = np.linspace(lats[5], lats[10], 21)
    target_lons = np.linspace(lons[8], lons[13], 21)
    target_lats, target_lons = np.meshgrid(target_lats, target_lons, indexing="ij")

    via_flag = horizontal_interpolate(
        obj,
        target_lats=target_lats,
        target_lons=target_lons,
        method="quadratic12",
        shape_preserving=True,
    )
    via_alias = horizontal_interpolate(
        obj,
        target_lats=target_lats,
        target_lons=target_lons,
        method="quadratic12_monotonic",
    )

    np.testing.assert_allclose(via_alias.values, via_flag.values, atol=1.0e-14)


def test_quadratic12_monotonic_supports_regular_gaussian_target_grid() -> None:
    lats = gaussian_latitudes(16)
    lons = regular_longitudes(32)
    values = np.zeros((16, 32), dtype=np.float64)
    values[7, 10] = 1.0
    obj = xr.DataArray(
        values,
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
        attrs={"GRIB_N": 8, "GRIB_gridType": "regular_gg"},
        name="q",
    )

    out = horizontal_interpolate(
        obj,
        source_grid="F8",
        target_grid="F16",
        method="quadratic12_monotonic",
    )

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (32, 64)
    assert out.attrs["GRIB_N"] == 16
    assert out.attrs["GRIB_gridType"] == "regular_gg"
    assert out.attrs["GRIB_numberOfPoints"] == 32 * 64
    assert "method=quadratic12_monotonic" in out.attrs["history"]
    assert float(out.min()) >= 0.0
    assert float(out.max()) <= 1.0


def test_native_field_major_batch_matches_legacy_batch_layout() -> None:
    """The field-major fast path must be only a layout optimization."""
    from fullpos.interpolation.kernels import (
        _prepare_horizontal_row_batch,
        _regular_stencil_safe_mask,
    )

    add_native_runtime_dir()
    import fullpos._ectrans as _ectrans

    source_grid = parse_grid("O32")
    target_grid = parse_grid("F48")
    values = np.arange(5 * source_grid.size, dtype=np.float64).reshape(5, source_grid.size)
    values = values / 1000.0
    fields, nloen, source_lats = _prepare_horizontal_row_batch(
        values,
        source_lats=None,
        source_lons=None,
        source_pl=source_grid.pl,
        source_grid=None,
    )
    target_lats, target_lons = np.meshgrid(
        gaussian_latitudes(target_grid.nlat),
        regular_longitudes(target_grid.work_nlon),
        indexing="ij",
    )
    flat_target_lats = target_lats.reshape(-1)
    flat_target_lons = target_lons.reshape(-1)
    safe = _regular_stencil_safe_mask(flat_target_lats, source_lats)

    for method, monotonic in (
        ("bilinear", False),
        ("quadratic12", False),
        ("quadratic12", True),
    ):
        legacy = _ectrans.horizontal_regular_kernel_batch(
            np.asfortranarray(fields.T, dtype=np.float64),
            np.asfortranarray(nloen, dtype=np.int32),
            np.asfortranarray(np.deg2rad(source_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_target_lats[safe]), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_target_lons[safe]), dtype=np.float64),
            method,
            int(monotonic),
        )
        field_major = _ectrans.horizontal_regular_kernel_batch_field_major(
            np.ascontiguousarray(fields, dtype=np.float64),
            np.asfortranarray(nloen, dtype=np.int32),
            np.asfortranarray(np.deg2rad(source_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_target_lats[safe]), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_target_lons[safe]), dtype=np.float64),
            method,
            int(monotonic),
        )

        np.testing.assert_array_equal(field_major, np.asarray(legacy).T)


def test_quadratic12_monotonic_supports_packed_o_grid_to_regular_gaussian_target() -> None:
    pl = octahedral_pl(8)
    rows = []
    for row, row_nlon in enumerate(pl):
        values = np.zeros(int(row_nlon), dtype=np.float64)
        if row == 7:
            values[10] = 1.0
        rows.append(values)
    obj = xr.DataArray(
        np.concatenate(rows),
        dims=("values",),
        attrs={"GRIB_pl": pl, "GRIB_N": 8, "GRIB_gridType": "reduced_gg"},
        name="q",
    )

    out = horizontal_interpolate(
        obj,
        source_grid="O8",
        target_grid="F16",
        method="quadratic12_monotonic",
    )

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (32, 64)
    assert out.attrs["GRIB_N"] == 16
    assert out.attrs["GRIB_gridType"] == "regular_gg"
    assert "GRIB_pl" not in out.attrs
    assert float(out.min()) >= 0.0
    assert float(out.max()) <= 1.0


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


def test_horizontal_interpolate_n_alias_uses_grib_pl_when_present() -> None:
    pl = octahedral_pl(4)
    lats = gaussian_latitudes(pl.size)
    rows = [row * 1000.0 + np.arange(row_nlon, dtype=np.float64) for row, row_nlon in enumerate(pl)]
    obj = xr.DataArray(
        np.concatenate(rows),
        dims=("values",),
        attrs={"GRIB_pl": pl, "GRIB_N": 4, "GRIB_gridType": "reduced_gg"},
        name="q",
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
        source_grid="N4",
        target_lats=np.concatenate(target_lats),
        target_lons=np.concatenate(target_lons),
        method="quadratic12_monotonic",
    )

    np.testing.assert_allclose(out.values, np.concatenate(expected), atol=1.0e-10)


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


def test_horizontal_interpolate_chunks_use_batched_native_path(monkeypatch) -> None:
    pl = octahedral_pl(8)
    rows = [row * 1000.0 + np.arange(row_nlon, dtype=np.float64) for row, row_nlon in enumerate(pl)]
    base = np.concatenate(rows)
    obj = xr.DataArray(
        np.stack([base, base + 10.0, base * 0.5]),
        dims=("hybrid", "values"),
        coords={"hybrid": [1, 2, 3]},
        attrs={"GRIB_pl": pl, "GRIB_N": 8, "GRIB_gridType": "reduced_gg"},
        name="q",
    )

    from fullpos.interpolation import horizontal as horizontal_module

    calls = {"batch": 0, "single": 0}
    original_batch = horizontal_module.horizontal_regular_kernel_batch
    original_single = horizontal_module.horizontal_regular_kernel

    def counted_batch(*args, **kwargs):
        calls["batch"] += 1
        return original_batch(*args, **kwargs)

    def counted_single(*args, **kwargs):
        calls["single"] += 1
        return original_single(*args, **kwargs)

    monkeypatch.setattr(
        horizontal_module,
        "horizontal_regular_kernel_batch",
        counted_batch,
    )
    monkeypatch.setattr(
        horizontal_module,
        "horizontal_regular_kernel",
        counted_single,
    )

    batched = horizontal_interpolate(
        obj,
        source_grid="O8",
        target_grid="F16",
        method="quadratic12_monotonic",
        chunks={"hybrid": 3},
    )
    single_fields = [
        horizontal_interpolate(
            obj.isel(hybrid=i),
            source_grid="O8",
            target_grid="F16",
            method="quadratic12_monotonic",
        ).values
        for i in range(obj.sizes["hybrid"])
    ]

    assert batched.dims == ("hybrid", "latitude", "longitude")
    np.testing.assert_allclose(batched.values, np.stack(single_fields), atol=1.0e-10)
    assert calls == {"batch": 1, "single": obj.sizes["hybrid"]}


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


def test_source_mask_is_supported_for_native_average_and_outputs_nan() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    values = np.arange(8 * 16, dtype=np.float64).reshape(8, 16)
    obj = xr.DataArray(
        values,
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )
    source_mask = np.ones_like(values, dtype=bool)
    source_mask[3, 0] = False
    source_mask[3, 1] = False
    source_mask[4, 0] = False
    source_mask[4, 1] = False

    out = horizontal_interpolate(
        obj,
        target_lats=np.array([lats[3]]),
        target_lons=np.array([lons[0]]),
        method="average",
        source_mask=source_mask,
    )

    assert np.isnan(out.values[0])


def test_target_mask_forces_masked_output_to_nan() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    obj = xr.DataArray(
        np.ones((8, 16), dtype=np.float64),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )

    out = horizontal_interpolate(
        obj,
        target_lats=np.array([lats[3], lats[3]]),
        target_lons=np.array([lons[0], lons[1]]),
        method="nearest",
        source_mask=np.ones_like(obj.values, dtype=bool),
        target_mask=np.array([True, False]),
    )

    assert np.isfinite(out.values[0])
    assert np.isnan(out.values[1])


def test_source_mask_is_rejected_for_fpint4_fpint12_paths() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    obj = xr.DataArray(
        np.ones((8, 16), dtype=np.float64),
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
    )

    with pytest.raises(ValueError, match="source_mask"):
        horizontal_interpolate(
            obj,
            target_lats=np.array([lats[3]]),
            target_lons=np.array([lons[0]]),
            source_mask=np.ones_like(obj.values, dtype=bool),
        )


def test_land_sea_mask_to_grid_samples_regular_ll_mask_to_packed_o_grid() -> None:
    lsm = xr.DataArray(
        np.array([[1.0, 0.0], [1.0, 0.0]]),
        dims=("latitude", "longitude"),
        coords={"latitude": [45.0, -45.0], "longitude": [0.0, 180.0]},
    )

    sea = land_sea_mask_to_grid(lsm, target_grid="O2", kind="sea")
    land = land_sea_mask_to_grid(lsm, target_grid="O2", kind="land")

    assert sea.dims == ("values",)
    assert sea.size == int(octahedral_pl(2).sum())
    assert bool(sea.any())
    assert bool(land.any())
    np.testing.assert_array_equal(sea.values, ~land.values)


def test_masked_surface_interpolate_wraps_land_sea_masks() -> None:
    lats = gaussian_latitudes(8)
    lons = regular_longitudes(16)
    values = np.full((8, 16), 280.0, dtype=np.float64)
    values[3:5, 0:2] = np.nan
    field = xr.DataArray(
        values,
        dims=("latitude", "longitude"),
        coords={"latitude": lats, "longitude": lons},
        name="sst",
    )
    lsm = xr.DataArray(
        np.array([[1.0, 0.0], [1.0, 0.0]]),
        dims=("latitude", "longitude"),
        coords={"latitude": [45.0, -45.0], "longitude": [0.0, 180.0]},
    )

    out = masked_surface_interpolate(
        field,
        land_sea_mask=lsm,
        source_grid="F4",
        target_grid="F4",
        method="average",
    )

    assert out.name == "sst"
    assert out.attrs["fullpos_surface_mask_kind"] == "sea"
    assert np.isnan(out.values).any()
    assert np.isfinite(out.values).any()


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


def test_horizontal_interpolate_rejects_shape_preserving_for_non_fpint12_methods() -> None:
    obj = xr.DataArray(
        np.ones((8, 16)),
        dims=("latitude", "longitude"),
        coords={"latitude": gaussian_latitudes(8), "longitude": regular_longitudes(16)},
    )

    with pytest.raises(ValueError, match="shape_preserving=True"):
        horizontal_interpolate(
            obj,
            target_lats=np.array([0.0]),
            target_lons=np.array([0.0]),
            method="bilinear",
            shape_preserving=True,
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
