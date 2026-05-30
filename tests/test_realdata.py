from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from fullpos import (
    horizontal_interpolate,
    masked_surface_interpolate,
    regrid,
    spectral_filter,
    vertical_interpolate,
)
from fullpos._vertical.pressure import (
    _infer_reference_pressure_from_ak,
    _midlevel_coefficients_from_half_levels,
    prepare_pressure_request,
)
from fullpos._vertical.validation import pressure_metric_summary


REAL_O320 = Path(
    r"L:\ERA5_test\era5_reanalysis_model_level_20250102_packing_CCSDS_O320.grib2"
)
REAL_SURFACE_O96 = Path(
    r"L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19950710_hourly_O96.grib"
)
REAL_SURFACE_O96_MATCHING = Path(
    r"L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19781201_hourly_O96.grib"
)
REAL_OROGRAPHY_O96 = Path(
    r"L:\ERA5_Complete\ERA5_Reanalysis_surface_geopotential_O96.grib"
)
REAL_LSM = Path(
    r"L:\ERA5_Complete\ERA5_Land_Sea_Mask.grib"
)
REAL_MODEL_O96 = Path(
    r"L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2"
)
REAL_MODEL_F96 = Path(
    r"L:\ERA5_test\era5_reanalysis_model_level_20250101_packing_CCSDS_F96.grib2"
)


def _optional_real_era5_pressure_paths() -> tuple[Path, Path]:
    pytest.importorskip("cfgrib")
    model = Path(os.environ.get("FULLPOS_ERA5_MODEL_FILE") or REAL_MODEL_O96)
    surface = Path(os.environ.get("FULLPOS_ERA5_SURFACE_FILE") or REAL_SURFACE_O96_MATCHING)
    if not model.exists() or not surface.exists():
        pytest.skip(
            "local ERA5 pressure smoke inputs not found; set FULLPOS_ERA5_MODEL_FILE "
            "and FULLPOS_ERA5_SURFACE_FILE or place the expected O96 samples locally"
        )
    return model, surface


@pytest.mark.skipif(not REAL_O320.exists(), reason="local ERA5 O320 GRIB sample not found")
def test_real_o320_temperature_to_o480_has_finite_packed_output() -> None:
    ds = xr.open_dataset(
        REAL_O320,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    field = ds["t"].isel(time=0, hybrid=0)

    out = regrid(field, source_grid="O320", target_grid="O480")

    assert out.dims == ("values",)
    assert out.size == 938880
    assert np.isfinite(out.values).all()


@pytest.mark.skipif(not REAL_O320.exists(), reason="local ERA5 O320 GRIB sample not found")
def test_real_o320_temperature_to_ll1_has_regular_latlon_output() -> None:
    ds = xr.open_dataset(
        REAL_O320,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    field = ds["t"].isel(time=0, hybrid=0)

    out = regrid(field, source_grid="O320", target_grid="LL1.0", chunk_size=1)

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (180, 360)
    assert out.attrs["GRIB_gridType"] == "regular_ll"
    np.testing.assert_allclose(out["latitude"].values[[0, -1]], [89.5, -89.5])
    np.testing.assert_allclose(out["longitude"].values[[0, -1]], [0.0, 359.0])
    assert np.isfinite(out.values).all()


@pytest.mark.skipif(not REAL_SURFACE_O96.exists(), reason="local ERA5 O96 surface GRIB sample not found")
def test_real_o96_sst_bitmap_is_rejected_for_spectral_regridding() -> None:
    ds = xr.open_dataset(
        REAL_SURFACE_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sst", "typeOfLevel": "surface"},
            "read_keys": [
                "gridType",
                "N",
                "pl",
                "numberOfPoints",
                "bitmapPresent",
                "numberOfMissing",
            ],
        },
    )
    field = ds["sst"].isel(time=0) if "time" in ds["sst"].dims else ds["sst"]

    assert field.attrs.get("GRIB_bitmapPresent") == 1
    assert not np.isfinite(field.values).all()
    with pytest.raises(ValueError, match="masked surface interpolation"):
        regrid(field, target_grid="N160")


@pytest.mark.skipif(not REAL_SURFACE_O96.exists(), reason="local ERA5 O96 surface GRIB sample not found")
def test_real_o96_sst_masked_regridding_has_finite_sea_points() -> None:
    ds = xr.open_dataset(
        REAL_SURFACE_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sst", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "bitmapPresent"],
        },
    )
    field = ds["sst"].isel(time=0) if "time" in ds["sst"].dims else ds["sst"]

    out = regrid(field, target_grid="N160", method="masked", chunk_size=8)

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (320, 640)
    assert np.isfinite(out.values).sum() > np.isfinite(field.values).sum()
    assert np.isnan(out.values).any()


@pytest.mark.skipif(
    not (REAL_SURFACE_O96_MATCHING.exists() and REAL_LSM.exists()),
    reason="local ERA5 O96 SST and land-sea-mask GRIB samples not found",
)
def test_real_o96_sst_land_sea_masked_average_to_f160_preserves_sea_only_output() -> None:
    surface = xr.open_dataset(
        REAL_SURFACE_O96_MATCHING,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sst", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "bitmapPresent"],
        },
    )
    lsm = xr.open_dataset(
        REAL_LSM,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "lsm"},
            "read_keys": ["gridType", "N", "numberOfPoints"],
        },
    )
    field = surface["sst"].isel(time=0) if "time" in surface["sst"].dims else surface["sst"]
    land_sea_mask = lsm["lsm"].isel(time=0) if "time" in lsm["lsm"].dims else lsm["lsm"]

    out = masked_surface_interpolate(
        field,
        land_sea_mask=land_sea_mask,
        source_grid="O96",
        target_grid="F160",
        kind="sea",
        method="average",
    )

    assert out.dims == ("latitude", "longitude")
    assert out.shape == (320, 640)
    assert out.attrs["fullpos_surface_mask_kind"] == "sea"
    assert out.attrs["fullpos_surface_mask_source_grid"] == "O96"
    assert out.attrs["fullpos_surface_mask_target_grid"] == "F160"
    assert np.isfinite(out.values).sum() > 0
    assert np.isnan(out.values).sum() > 0
    assert not np.any(np.abs(out.values[np.isfinite(out.values)]) >= 1.0e10)
    assert 250.0 < float(np.nanmin(out.values)) < 320.0
    assert 250.0 < float(np.nanmean(out.values)) < 320.0
    assert 250.0 < float(np.nanmax(out.values)) < 320.0


@pytest.mark.skipif(not REAL_MODEL_O96.exists(), reason="local ERA5 O96 model-level GRIB sample not found")
def test_real_o96_temperature_spectral_filter_has_finite_output() -> None:
    ds = xr.open_dataset(
        REAL_MODEL_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    field = ds["t"].isel(time=0, hybrid=0)

    out = spectral_filter(field, grid="O96", ntrunc=47, chunk_size=8)
    diff = out.values - field.values

    assert out.dims == ("values",)
    assert out.size == 40320
    assert np.isfinite(out.values).all()
    assert 0.0 < float(np.nanmean(np.abs(diff))) < 10.0


@pytest.mark.skipif(not REAL_MODEL_F96.exists(), reason="local ERA5 F96 model-level GRIB sample not found")
def test_real_f96_temperature_fullpos_horizontal_interpolation_matches_grid_points() -> None:
    ds = xr.open_dataset(
        REAL_MODEL_F96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "numberOfPoints", "packingType"],
        },
    )
    field = ds["t"].isel(time=0, hybrid=0).load()
    lats = field.latitude.values
    lons = field.longitude.values
    target_lats, target_lons = np.meshgrid(lats[2:-2], lons, indexing="ij")

    expected = field.isel(latitude=slice(2, -2)).values
    for method in ("nearest", "bilinear", "quadratic12"):
        out = horizontal_interpolate(
            field,
            target_lats=target_lats,
            target_lons=target_lons,
            method=method,
        )
        np.testing.assert_allclose(out.values, expected, atol=1.0e-10)

    averaged = horizontal_interpolate(
        field,
        target_lats=target_lats,
        target_lons=target_lons,
        method="average",
        average_radius=1,
    )
    values = field.values
    columns = np.arange(values.shape[1])
    expected_average = np.empty_like(expected)
    for out_row, src_row in enumerate(range(2, values.shape[0] - 2)):
        expected_average[out_row, :] = (
            values[src_row, columns]
            + values[src_row, (columns + 1) % values.shape[1]]
            + values[src_row + 1, columns]
            + values[src_row + 1, (columns + 1) % values.shape[1]]
        ) / 4.0
    np.testing.assert_allclose(averaged.values, expected_average, atol=1.0e-10)


@pytest.mark.skipif(
    not (REAL_MODEL_O96.exists() and REAL_SURFACE_O96_MATCHING.exists()),
    reason="local ERA5 O96 model/surface GRIB samples not found",
)
def test_real_o96_pressure_interpolation_uses_native_fullpos_backend() -> None:
    model = xr.open_dataset(
        REAL_MODEL_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "pv", "packingType"],
        },
    )["t"].load()
    surface = xr.open_dataset(
        REAL_SURFACE_O96_MATCHING,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sp", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints"],
        },
    )["sp"].load()

    out = vertical_interpolate(
        model,
        target="pressure",
        levels=[100000.0, 85000.0, 50000.0],
        chunks={"time": 1, "values": 10000},
        surface_pressure=surface,
    )

    assert out.dims == ("time", "pressure", "values")
    assert out.shape == (4, 3, 40320)
    assert out.attrs["vertical_backend"] == "FULLPOS"
    assert out.attrs["vertical_target"] == "pressure"
    assert out.attrs["pressure_units"] == "Pa"
    np.testing.assert_allclose(out["pressure"].values, [100000.0, 85000.0, 50000.0])
    assert np.isfinite(out.values).all()


def test_real_era5_pressure_vertical_dataset_smoke_uses_apache_for_tuvq() -> None:
    model_path, surface_path = _optional_real_era5_pressure_paths()

    model = xr.open_dataset(
        model_path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "pv", "packingType"],
        },
    )[["t", "u", "v", "q"]].isel(time=slice(0, 1)).load()
    surface = xr.open_dataset(
        surface_path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sp", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints"],
        },
    )["sp"].sel(time=model["time"]).load()

    levels = [20000.0, 30000.0, 50000.0]
    out = vertical_interpolate(
        model,
        target="pressure",
        levels=levels,
        variables=["t", "u", "v", "q"],
        chunks={"time": 1, "values": 10000},
        surface_pressure=surface,
    )
    chunked = vertical_interpolate(
        model,
        target="pressure",
        levels=levels,
        variables=["t", "u", "v", "q"],
        chunks={"time": 1, "values": 2000},
        surface_pressure=surface,
    )

    assert set(out.data_vars) == {"t", "u", "v", "q"}
    np.testing.assert_allclose(out["pressure"].values, levels)
    for name in ("t", "u", "v", "q"):
        assert out[name].dims == ("time", "pressure", "values")
        assert out[name].attrs["vertical_backend"] == "FULLPOS"
        assert out[name].attrs["vertical_target"] == "pressure"
        assert out[name].attrs["vertical_native_path"] == "APACHE"
        assert out[name].attrs["pressure_units"] == "Pa"
        finite_fraction = float(np.isfinite(out[name].values).mean())
        assert finite_fraction > 0.99

    for name in ("t", "u", "v", "q"):
        metrics = pressure_metric_summary(chunked[name], out[name], level_dim="pressure")
        assert metrics["overall"]["count"] == out[name].size
        assert metrics["overall"]["rmse"] <= 1.0e-8
        assert metrics["overall"]["max_abs"] <= 1.0e-8


@pytest.mark.skipif(
    not (REAL_MODEL_O96.exists() and REAL_SURFACE_O96_MATCHING.exists()),
    reason="local ERA5 O96 model/surface GRIB samples not found",
)
def test_real_o96_pressure_dataset_uses_apache_lescale() -> None:
    model = xr.open_dataset(
        REAL_MODEL_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "pv", "packingType"],
        },
    )[["t", "u", "v", "q"]].isel(time=slice(0, 1)).load()
    surface = xr.open_dataset(
        REAL_SURFACE_O96_MATCHING,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sp", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints"],
        },
    )["sp"].sel(time=model["time"]).load()
    target_surface = surface * 0.95

    out = vertical_interpolate(
        model,
        target="pressure",
        levels=[30000.0, 50000.0, 85000.0],
        variables=["t", "u", "v", "q"],
        chunks={"time": 1, "values": 10000},
        surface_pressure=surface,
        target_surface_pressure=target_surface,
        lescale=True,
    )
    chunked = vertical_interpolate(
        model,
        target="pressure",
        levels=[30000.0, 50000.0, 85000.0],
        variables=["t", "u", "v", "q"],
        chunks={"time": 1, "values": 2000},
        surface_pressure=surface,
        target_surface_pressure=target_surface,
        lescale=True,
    )

    assert set(out.data_vars) == {"t", "u", "v", "q"}
    for name in ("t", "u", "v", "q"):
        assert out[name].dims == ("time", "pressure", "values")
        assert out[name].attrs["vertical_backend"] == "FULLPOS"
        assert out[name].attrs["vertical_native_path"] == "APACHE"
        assert out[name].attrs["vertical_lescale"] == "enabled"
        assert np.isfinite(out[name].values).all()
        metrics = pressure_metric_summary(chunked[name], out[name], level_dim="pressure")
        assert metrics["overall"]["count"] == out[name].size
        assert metrics["overall"]["rmse"] <= 1.0e-8
        assert metrics["overall"]["max_abs"] <= 1.0e-8


@pytest.mark.skipif(
    not (REAL_MODEL_O96.exists() and REAL_SURFACE_O96_MATCHING.exists() and REAL_OROGRAPHY_O96.exists()),
    reason="local ERA5 O96 model/surface/orography GRIB samples not found",
)
def test_real_o96_flight_level_matches_equivalent_height_above_sea() -> None:
    model = xr.open_dataset(
        REAL_MODEL_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "pv", "packingType"],
        },
    )["t"].isel(time=slice(0, 1)).load()
    surface = xr.open_dataset(
        REAL_SURFACE_O96_MATCHING,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sp", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints"],
        },
    )["sp"].sel(time=model["time"]).load()
    orography = xr.open_dataset(
        REAL_OROGRAPHY_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "z", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints"],
        },
    )["z"].load()
    if "time" in orography.dims:
        orography = orography.isel(time=0, drop=True)

    flight = vertical_interpolate(
        model,
        target="flight_level",
        levels=[100.0, 200.0, 300.0, 350.0],
        chunks={"time": 1, "values": 10000},
        surface_pressure=surface,
        orography_geopotential=orography,
    )
    height = vertical_interpolate(
        model,
        target="height_above_sea",
        levels=[3048.0, 6096.0, 9144.0, 10668.0],
        chunks={"time": 1, "values": 10000},
        surface_pressure=surface,
        orography_geopotential=orography,
    )

    assert flight.dims == ("time", "flight_level", "values")
    assert flight.shape == (1, 4, 40320)
    assert flight.attrs["vertical_backend"] == "FULLPOS"
    assert flight.attrs["vertical_target"] == "flight_level"
    assert flight.attrs["flight_level_units"] == "hundreds of feet"
    assert np.isfinite(flight.values).all()
    np.testing.assert_allclose(flight["flight_level"].values, [100.0, 200.0, 300.0, 350.0])
    np.testing.assert_allclose(flight.values, height.values)


@pytest.mark.skipif(
    not (REAL_MODEL_O96.exists() and REAL_SURFACE_O96_MATCHING.exists()),
    reason="local ERA5 O96 model/surface GRIB samples not found",
)
def test_real_o96_skyborn_pressure_reference_runs_on_packed_era5_input() -> None:
    pytest.importorskip("skyborn.interp.interpolation")
    from skyborn.interp.interpolation import interp_hybrid_to_pressure

    model = xr.open_dataset(
        REAL_MODEL_O96,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "pv", "packingType"],
        },
    )["t"].load()
    surface = xr.open_dataset(
        REAL_SURFACE_O96_MATCHING,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "sp", "typeOfLevel": "surface"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints"],
        },
    )["sp"].load()

    request = prepare_pressure_request(
        model,
        levels=[100000.0, 85000.0, 50000.0],
        surface_pressure=surface,
    )
    hyam_values, hybm_values = _midlevel_coefficients_from_half_levels(
        request.ak,
        request.bk,
    )
    p0 = _infer_reference_pressure_from_ak(request.ak)
    hyam = xr.DataArray(
        hyam_values,
        dims=(request.hybrid_dim,),
        coords={request.hybrid_dim: model[request.hybrid_dim].values},
        name="hyam",
    )
    hybm = xr.DataArray(
        hybm_values,
        dims=(request.hybrid_dim,),
        coords={request.hybrid_dim: model[request.hybrid_dim].values},
        name="hybm",
    )

    out = interp_hybrid_to_pressure(
        data=model,
        ps=request.surface_pressure,
        hyam=hyam,
        hybm=hybm,
        p0=p0,
        new_levels=request.levels,
        lev_dim=request.hybrid_dim,
        method="linear",
    )

    assert out.dims == ("time", "plev", "values")
    assert out.shape == (4, 3, 40320)
    assert p0 == 1.0
    np.testing.assert_allclose(out["plev"].values, request.levels)
    assert out.attrs["GRIB_gridType"] == "reduced_gg"
    finite_counts = [int(np.isfinite(out.isel(plev=index).values).sum()) for index in range(out.sizes["plev"])]
    assert 0 < finite_counts[0] < out.isel(plev=0).size
    assert finite_counts[0] < finite_counts[1] < finite_counts[2]
    assert finite_counts[2] == out.isel(plev=2).size
