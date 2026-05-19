import pytest
import xarray as xr
import numpy as np

from fullpos import (
    Regridder,
    regrid,
    regrid_values,
    spectral_interpolate,
    spectral_interpolate_values,
    spectral_regrid,
    spectral_regrid_values,
)
from fullpos.grids import gaussian_latitudes, octahedral_pl, regular_longitudes


def test_regrid_requires_target_grid() -> None:
    with pytest.raises(ValueError, match="target_grid is required"):
        regrid(object())


def test_regrid_rejects_unsupported_method() -> None:
    with pytest.raises(ValueError, match="method must be"):
        regrid(object(), target_grid="O480", method="nearest")


def test_regrid_rejects_non_xarray_input() -> None:
    with pytest.raises(TypeError, match="xarray DataArray or Dataset"):
        regrid(object(), target_grid="O480")


def test_regrid_accepts_small_xarray_input_until_backend_dependency() -> None:
    obj = xr.DataArray(np.ones((8, 16), dtype=np.float32), dims=("latitude", "longitude"))
    out = regrid(obj, source_grid="N4", target_grid="N4")
    assert out.shape == obj.shape
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regrid_supports_small_regular_to_octahedral() -> None:
    obj = xr.DataArray(np.ones((8, 16), dtype=np.float32), dims=("latitude", "longitude"))
    out = regrid(obj, source_grid="N4", target_grid="O4")
    assert out.dims == ("values",)
    assert out.size == 208
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regrid_supports_small_octahedral_to_regular() -> None:
    obj = xr.DataArray(np.ones(208, dtype=np.float32), dims=("values",))
    out = regrid(obj, source_grid="O4", target_grid="N4")
    assert out.dims == ("latitude", "longitude")
    assert out.shape == (8, 16)
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regrid_supports_regular_latlon_target() -> None:
    obj = xr.DataArray(
        np.ones((8, 16), dtype=np.float32),
        dims=("latitude", "longitude"),
        coords={
            "latitude": gaussian_latitudes(8),
            "longitude": regular_longitudes(16),
        },
        attrs={"units": "K"},
    )

    out = regrid(obj, source_grid="F4", target_grid="LL1.0")

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (180, 360)
    assert out.attrs["GRIB_gridType"] == "regular_ll"
    assert out.attrs["GRIB_numberOfPoints"] == 180 * 360
    assert out.attrs["GRIB_iDirectionIncrementInDegrees"] == 1.0
    assert out.attrs["GRIB_jDirectionIncrementInDegrees"] == 1.0
    assert "target_grid=LL1" in out.attrs["history"]
    np.testing.assert_allclose(out["latitude"].values[[0, -1]], [89.5, -89.5])
    np.testing.assert_allclose(out["longitude"].values[[0, -1]], [0.0, 359.0])
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regrid_dataset_regular_latlon_target_selects_variables() -> None:
    ds = xr.Dataset(
        data_vars={
            "t": (("latitude", "longitude"), np.ones((8, 16), dtype=np.float32)),
            "surface": ("time", np.arange(2, dtype=np.float32)),
        },
        coords={
            "latitude": gaussian_latitudes(8),
            "longitude": regular_longitudes(16),
            "time": [0, 1],
        },
    )

    out = regrid(ds, source_grid="F4", target_grid="LL2.5", variables=["t"])

    assert out["t"].dims == ("latitude", "longitude")
    assert out["t"].shape == (72, 144)
    assert "surface" not in out
    assert out.attrs["GRIB_gridType"] == "regular_ll"
    assert out["t"].attrs["GRIB_gridType"] == "regular_ll"
    np.testing.assert_allclose(out["t"].values, 1.0, atol=1.0e-5)


def test_regrid_regular_latlon_target_accepts_named_chunks() -> None:
    obj = xr.DataArray(
        np.ones((2, 3, 8, 16), dtype=np.float32),
        dims=("time", "hybrid", "latitude", "longitude"),
        coords={
            "time": [0, 1],
            "hybrid": [1, 2, 3],
            "latitude": gaussian_latitudes(8),
            "longitude": regular_longitudes(16),
        },
    )

    out = regrid(
        obj,
        source_grid="F4",
        target_grid="LL2.5",
        chunks={"time": 1, "hybrid": 2},
    )

    assert out.dims == ("time", "hybrid", "latitude", "longitude")
    assert out.shape == (2, 3, 72, 144)
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regrid_infers_octahedral_from_grib_attrs() -> None:
    obj = xr.DataArray(
        np.ones(208, dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": octahedral_pl(4),
        },
    )
    out = regrid(obj, target_grid="O4", chunk_size=2)
    assert out.dims == ("values",)
    assert out.size == 208
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regrid_updates_output_grid_attrs() -> None:
    obj = xr.DataArray(
        np.ones(208, dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": octahedral_pl(4),
            "GRIB_numberOfPoints": 208,
        },
    )
    out = regrid(obj, target_grid="N4")
    assert out.attrs["GRIB_gridType"] == "regular_gg"
    assert out.attrs["GRIB_N"] == 4
    assert out.attrs["GRIB_numberOfPoints"] == 128
    assert "GRIB_pl" not in out.attrs
    assert "fullpos regrid" in out.attrs["history"]
    assert "source_grid=O4" in out.attrs["history"]
    assert "target_grid=N4" in out.attrs["history"]


def test_regrid_prepends_existing_history() -> None:
    obj = xr.DataArray(
        np.ones(208, dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": octahedral_pl(4),
            "history": "previous operation",
        },
    )

    out = regrid(obj, target_grid="N4")

    assert out.attrs["history"].splitlines()[1] == "previous operation"


def test_regrid_rejects_bad_reduced_grid_metadata() -> None:
    obj = xr.DataArray(
        np.ones(128, dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": np.full(8, 16),
        },
    )
    with pytest.raises(ValueError, match="not an ECMWF octahedral O-grid"):
        regrid(obj, target_grid="O4")


def test_regrid_values_supports_leading_dims_and_chunks() -> None:
    values = np.ones((3, 208), dtype=np.float32)
    out = regrid_values(values, source_grid="O4", target_grid="N4", chunk_size=2)
    assert out.shape == (3, 8, 16)
    np.testing.assert_allclose(out, 1.0, atol=1.0e-5)


def test_regrid_values_rejects_invalid_chunk_size() -> None:
    values = np.ones((2, 208), dtype=np.float32)
    with pytest.raises(ValueError, match="chunk_size"):
        regrid_values(values, source_grid="O4", target_grid="O4", chunk_size=0)


def test_regrid_values_rejects_missing_values_by_default() -> None:
    values = np.ones((2, 208), dtype=np.float32)
    values[0, 0] = np.nan

    with pytest.raises(ValueError, match="requires finite input"):
        regrid_values(values, source_grid="O4", target_grid="O4", chunk_size=1)


def test_regrid_values_can_preserve_nan_propagating_backend_behavior() -> None:
    values = np.ones((1, 208), dtype=np.float32)
    values[0, 0] = np.nan

    out = regrid_values(
        values,
        source_grid="O4",
        target_grid="O4",
        chunk_size=1,
        missing_policy="ignore",
    )

    assert np.isnan(out).any()


def test_spectral_interpolate_supports_gaussian_targets() -> None:
    obj = xr.DataArray(np.ones((8, 16), dtype=np.float32), dims=("latitude", "longitude"))

    out = spectral_interpolate(obj, source_grid="N4", target_grid="O4")

    assert out.dims == ("values",)
    assert out.size == 208
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_spectral_interpolate_rejects_regular_latlon_targets() -> None:
    obj = xr.DataArray(np.ones((8, 16), dtype=np.float32), dims=("latitude", "longitude"))

    with pytest.raises(ValueError, match="only supports Gaussian"):
        spectral_interpolate(obj, source_grid="F4", target_grid="LL1.0")


def test_spectral_interpolate_values_supports_gaussian_targets() -> None:
    values = np.ones((3, 208), dtype=np.float32)

    out = spectral_interpolate_values(values, source_grid="O4", target_grid="N4", chunk_size=2)

    assert out.shape == (3, 8, 16)
    np.testing.assert_allclose(out, 1.0, atol=1.0e-5)


def test_spectral_regrid_alias_supports_gaussian_targets() -> None:
    obj = xr.DataArray(np.ones((8, 16), dtype=np.float32), dims=("latitude", "longitude"))

    out = spectral_regrid(obj, source_grid="N4", target_grid="O4")

    assert out.dims == ("values",)
    assert out.size == 208
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_spectral_regrid_values_alias_supports_gaussian_targets() -> None:
    values = np.ones((3, 208), dtype=np.float32)

    out = spectral_regrid_values(values, source_grid="O4", target_grid="N4", chunk_size=2)

    assert out.shape == (3, 8, 16)
    np.testing.assert_allclose(out, 1.0, atol=1.0e-5)


def test_masked_regrid_values_preserves_constant_with_missing_points() -> None:
    values = np.ones(208, dtype=np.float32)
    values[::7] = np.nan

    out = regrid_values(values, source_grid="O4", target_grid="N4", method="masked")

    assert out.shape == (8, 16)
    assert np.isfinite(out).any()
    np.testing.assert_allclose(out[np.isfinite(out)], 1.0)


def test_masked_regrid_data_array_marks_output_method() -> None:
    obj = xr.DataArray(
        np.ones(208, dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": octahedral_pl(4),
        },
    )
    obj.values[::7] = np.nan

    out = regrid(obj, target_grid="N4", method="masked")

    assert out.dims == ("latitude", "longitude")
    assert out.attrs["fullpos_regrid_method"] == "masked"
    assert np.isfinite(out.values).any()


def test_regrid_values_restores_non_trailing_horizontal_axis() -> None:
    values = np.ones((208, 3), dtype=np.float32)
    out = regrid_values(values, source_grid="O4", target_grid="N4", axis=0, chunk_size=2)
    assert out.shape == (8, 16, 3)
    np.testing.assert_allclose(out, 1.0, atol=1.0e-5)


def test_regridder_reuses_grid_pair_for_numpy_values() -> None:
    regridder = Regridder("O4", "N4", chunk_size=2)
    values = np.ones((3, 208), dtype=np.float32)
    out = regridder.regrid_values(values)
    assert out.shape == (3, 8, 16)
    np.testing.assert_allclose(out, 1.0, atol=1.0e-5)


def test_regridder_reuses_grid_pair_for_data_array() -> None:
    regridder = Regridder("O4", "N4", chunk_size=2)
    obj = xr.DataArray(np.ones((3, 208), dtype=np.float32), dims=("hybrid", "values"))
    out = regridder.regrid_data_array(obj)
    assert out.dims == ("hybrid", "latitude", "longitude")
    assert out.shape == (3, 8, 16)
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regridder_exposes_grid_metadata() -> None:
    regridder = Regridder("O4", "N4", chunk_size=2)

    assert regridder.source_shape == (208,)
    assert regridder.target_shape == (8, 16)
    assert regridder.source_dims == ("values",)
    assert regridder.target_dims == ("latitude", "longitude")
    assert set(regridder.target_coords) == {"latitude", "longitude"}
    assert regridder.info["source_grid"] == "O4"
    assert regridder.info["target_grid"] == "N4"
    assert regridder.info["method"] == "linear"


def test_regridder_from_dataarray_infers_source_grid() -> None:
    obj = xr.DataArray(
        np.ones(208, dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": octahedral_pl(4),
        },
    )

    regridder = Regridder.from_dataarray(obj, target_grid="N4")
    out = regridder.regrid_data_array(obj)

    assert regridder.source_grid.name == "O4"
    assert out.shape == (8, 16)
    np.testing.assert_allclose(out.values, 1.0, atol=1.0e-5)


def test_regridder_from_dataarray_infers_classic_reduced_source_grid() -> None:
    pl = np.array([10, 12, 14, 16, 16, 14, 12, 10], dtype=np.int64)
    obj = xr.DataArray(
        np.ones(int(pl.sum()), dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": pl,
        },
    )

    regridder = Regridder.from_dataarray(obj, target_grid="F4")

    assert regridder.source_grid.name == "N4"
    assert regridder.source_grid.kind == "classic_reduced"
    assert regridder.source_grid.size == int(pl.sum())


def test_regridder_from_dataarray_uses_classic_reduced_grid_for_explicit_n_alias() -> None:
    pl = np.array([10, 12, 14, 16, 16, 14, 12, 10], dtype=np.int64)
    obj = xr.DataArray(
        np.ones(int(pl.sum()), dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": pl,
        },
    )

    regridder = Regridder.from_dataarray(obj, target_grid="F4", source_grid="N4")

    assert regridder.source_grid.name == "N4"
    assert regridder.source_grid.kind == "classic_reduced"
    assert regridder.target_grid.name == "F4"


def test_regridder_from_dataarray_uses_classic_reduced_grid_for_explicit_n_target() -> None:
    pl = np.array([10, 12, 14, 16, 16, 14, 12, 10], dtype=np.int64)
    obj = xr.DataArray(
        np.ones(int(pl.sum()), dtype=np.float32),
        dims=("values",),
        attrs={
            "GRIB_gridType": "reduced_gg",
            "GRIB_N": 4,
            "GRIB_pl": pl,
        },
    )

    regridder = Regridder.from_dataarray(obj, target_grid="N4", source_grid="N4")

    assert regridder.source_grid.kind == "classic_reduced"
    assert regridder.target_grid.name == "N4"
    assert regridder.target_grid.kind == "classic_reduced"
    assert regridder.target_grid.size == int(pl.sum())


def test_regridder_from_dataarray_uses_preserved_source_pl_for_explicit_n_target() -> None:
    pl = np.array([10, 12, 14, 16, 16, 14, 12, 10], dtype=np.int64)
    obj = xr.DataArray(
        np.ones((8, 16), dtype=np.float32),
        dims=("latitude", "longitude"),
        attrs={
            "GRIB_gridType": "regular_gg",
            "GRIB_N": 4,
            "fullpos_source_grid": "N4",
            "fullpos_source_grid_kind": "classic_reduced",
            "fullpos_source_GRIB_N": 4,
            "fullpos_source_GRIB_gridType": "reduced_gg",
            "fullpos_source_GRIB_numberOfPoints": int(pl.sum()),
            "fullpos_source_GRIB_pl": pl,
        },
    )

    regridder = Regridder.from_dataarray(obj, target_grid="N4", source_grid="F4")

    assert regridder.source_grid.name == "F4"
    assert regridder.target_grid.name == "N4"
    assert regridder.target_grid.kind == "classic_reduced"
    assert regridder.target_grid.size == int(pl.sum())


def test_regridder_from_dataset_infers_source_grid() -> None:
    ds = xr.Dataset(
        data_vars={
            "t": (
                ("values",),
                np.ones(208, dtype=np.float32),
                {
                    "GRIB_gridType": "reduced_gg",
                    "GRIB_N": 4,
                    "GRIB_pl": octahedral_pl(4),
                },
            )
        }
    )

    regridder = Regridder.from_dataset(ds, target_grid="N4")

    assert regridder.source_grid.name == "O4"
    assert regridder.target_shape == (8, 16)


def test_regridder_supports_masked_method_for_numpy_values() -> None:
    regridder = Regridder("O4", "N4", method="masked", chunk_size=2)
    values = np.ones(208, dtype=np.float32)
    values[::7] = np.nan

    out = regridder.regrid_values(values)

    assert out.shape == (8, 16)
    assert np.isfinite(out).any()
    np.testing.assert_allclose(out[np.isfinite(out)], 1.0)


def test_regridder_supports_masked_method_for_data_array() -> None:
    regridder = Regridder("O4", "N4", method="masked", chunk_size=2)
    obj = xr.DataArray(np.ones(208, dtype=np.float32), dims=("values",))
    obj.values[::7] = np.nan

    out = regridder.regrid_data_array(obj)

    assert out.dims == ("latitude", "longitude")
    assert out.attrs["fullpos_regrid_method"] == "masked"
    assert np.isfinite(out.values).any()


def test_regrid_dataset_regrids_multiple_horizontal_variables() -> None:
    ds = xr.Dataset(
        data_vars={
            "t": (("time", "hybrid", "values"), np.ones((2, 3, 208), dtype=np.float32)),
            "q": (("time", "hybrid", "values"), np.full((2, 3, 208), 2.0, dtype=np.float32)),
            "surface": ("time", np.arange(2, dtype=np.float32)),
        },
        coords={
            "time": [0, 1],
            "hybrid": [1, 2, 3],
            "values": np.arange(208),
        },
    )

    out = regrid(ds, source_grid="O4", target_grid="N4", chunk_size=2)

    assert out["t"].dims == ("time", "hybrid", "latitude", "longitude")
    assert out["q"].shape == (2, 3, 8, 16)
    assert out["surface"].dims == ("time",)
    assert "values" not in out.coords
    assert "latitude" in out.coords
    assert "longitude" in out.coords
    assert "fullpos regrid" in out.attrs["history"]
    assert "variables=t,q,surface" in out.attrs["history"]
    np.testing.assert_allclose(out["t"].values, 1.0, atol=1.0e-5)
    np.testing.assert_allclose(out["q"].values, 2.0, atol=1.0e-5)


def test_regrid_dataset_explicit_non_horizontal_variable_errors() -> None:
    ds = xr.Dataset(
        data_vars={"surface": ("time", np.arange(2, dtype=np.float32))},
        coords={"time": [0, 1]},
    )

    with pytest.raises(ValueError, match="surface"):
        regrid(ds, target_grid="N4", variables=["surface"])
