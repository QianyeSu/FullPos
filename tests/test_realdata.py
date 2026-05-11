from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from fullpos import horizontal_interpolate, regrid, spectral_filter, vertical_interpolate
from fullpos._vertical.pressure import (
    _infer_reference_pressure_from_ak,
    _midlevel_coefficients_from_half_levels,
    prepare_pressure_request,
)


REAL_O320 = Path(
    r"L:\ERA5_test\era5_reanalysis_model_level_20250102_packing_CCSDS_O320.grib2"
)
REAL_SURFACE_O96 = Path(
    r"L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19950710_hourly_O96.grib"
)
REAL_SURFACE_O96_MATCHING = Path(
    r"L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19781201_hourly_O96.grib"
)
REAL_MODEL_O96 = Path(
    r"L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2"
)
REAL_MODEL_F96 = Path(
    r"L:\ERA5_test\era5_reanalysis_model_level_20250101_packing_CCSDS_F96.grib2"
)


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
    np.testing.assert_allclose(out["pressure"].values, [100000.0, 85000.0, 50000.0])
    assert np.isfinite(out.values).all()


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
