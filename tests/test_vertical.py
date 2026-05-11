from __future__ import annotations

import numpy as np
import pytest
import xarray as xr

from fullpos import diagnose_potential_vorticity, vertical_capabilities, vertical_interpolate
from fullpos.grids import gaussian_latitudes, regular_longitudes
from fullpos._vertical.pressure import (
    _infer_reference_pressure_from_ak,
    _midlevel_coefficients_from_half_levels,
    _select_hybrid_coefficients_for_levels,
    prepare_pressure_request,
)
from fullpos._vertical.validation import metric_block, pressure_metric_summary


def test_vertical_capabilities_reports_native_pressure_target() -> None:
    capabilities = vertical_capabilities()

    assert capabilities["model_level"] == "native_pp_chain"
    assert capabilities["potential_temperature"] == "native_pp_chain"
    assert capabilities["potential_vorticity"] == "native_ppltp"
    assert capabilities["pressure"] == "native"
    assert capabilities["temperature"] == "native_pp_chain"


def test_pressure_vertical_interpolate_rejects_missing_levels() -> None:
    data = xr.DataArray(
        np.ones((2, 3, 4)),
        dims=("time", "hybrid", "values"),
    )

    with pytest.raises(TypeError):
        vertical_interpolate(data, target="pressure")


def test_pressure_vertical_interpolate_validates_levels() -> None:
    data = xr.DataArray(
        np.ones((2, 3, 4)),
        dims=("time", "hybrid", "values"),
    )

    with pytest.raises(ValueError, match="positive"):
        vertical_interpolate(data, target="pressure", levels=[1000.0, -500.0])


def test_prepare_pressure_request_extracts_grib_pv_and_aligns_surface_pressure() -> None:
    data = xr.DataArray(
        np.ones((2, 3, 4)),
        dims=("time", "hybrid", "values"),
        coords={
            "time": np.array(["2000-01-01T00:00:00", "2000-01-01T06:00:00"], dtype="datetime64[ns]"),
        },
        attrs={"GRIB_pv": np.array([0.0, 1.0, 2.0, 3.0, 0.0, 0.3, 0.6, 1.0])},
        name="t",
    )
    sp = xr.DataArray(
        np.arange(4 * 4, dtype=np.float64).reshape(4, 4),
        dims=("time", "values"),
        coords={
            "time": np.array(
                [
                    "2000-01-01T00:00:00",
                    "2000-01-01T03:00:00",
                    "2000-01-01T06:00:00",
                    "2000-01-01T09:00:00",
                ],
                dtype="datetime64[ns]",
            ),
        },
        name="sp",
    )

    request = prepare_pressure_request(
        data,
        levels=[100000.0, 85000.0, 50000.0],
        surface_pressure=sp,
        chunks={"time": 1, "hybrid": 3},
    )

    np.testing.assert_allclose(request.ak, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(request.bk, [0.0, 0.3, 0.6, 1.0])
    assert request.hybrid_dim == "hybrid"
    assert request.surface_pressure.shape == (2, 4)
    np.testing.assert_array_equal(
        request.surface_pressure["time"].values,
        data["time"].values,
    )


def test_prepare_pressure_request_accepts_external_hybrid_coefficients_and_lnsp() -> None:
    data = xr.DataArray(
        np.ones((2, 3, 4)),
        dims=("time", "hybrid", "values"),
        coords={
            "time": np.array(["2000-01-01T00:00:00", "2000-01-01T06:00:00"], dtype="datetime64[ns]"),
        },
        name="q",
    )
    lnsp = xr.DataArray(
        np.log(np.full((2, 4), 90000.0)),
        dims=("time", "values"),
        coords={"time": data["time"].values},
        name="lnsp",
        attrs={"GRIB_shortName": "lnsp"},
    )
    coeffs = xr.Dataset(
        {
            "hyai": ("nhyi", np.array([0.0, 1.0, 2.0, 3.0])),
            "hybi": ("nhyi", np.array([0.0, 0.3, 0.6, 1.0])),
        }
    )

    request = prepare_pressure_request(
        data,
        levels=[85000.0],
        surface_pressure=lnsp,
        hybrid_coefficients=coeffs,
    )

    np.testing.assert_allclose(request.surface_pressure.values, 90000.0)
    np.testing.assert_allclose(request.ak, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(request.bk, [0.0, 0.3, 0.6, 1.0])


def test_midlevel_coefficients_collapse_adjacent_half_levels() -> None:
    hyam, hybm = _midlevel_coefficients_from_half_levels(
        np.array([0.0, 1.0, 3.0, 7.0]),
        np.array([0.0, 0.2, 0.6, 1.0]),
    )

    np.testing.assert_allclose(hyam, [0.5, 2.0, 5.0])
    np.testing.assert_allclose(hybm, [0.1, 0.4, 0.8])


def test_hybrid_coefficients_are_sliced_for_contiguous_level_subset() -> None:
    data = xr.DataArray(
        np.ones((3, 2)),
        dims=("hybrid", "values"),
        coords={"hybrid": [2, 3, 4]},
    )
    ak, bk = _select_hybrid_coefficients_for_levels(
        data,
        hybrid_dim="hybrid",
        ak=np.arange(6, dtype=np.float64),
        bk=np.arange(10, 16, dtype=np.float64),
        expected_half_levels=4,
    )

    np.testing.assert_allclose(ak, [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(bk, [11.0, 12.0, 13.0, 14.0])


def test_infer_reference_pressure_detects_pressure_unit_hybrid_a() -> None:
    assert _infer_reference_pressure_from_ak(np.array([0.0, 200.0, 5000.0])) == 1.0
    assert _infer_reference_pressure_from_ak(np.array([0.0, 0.1, 0.5])) == 100000.0
    assert _infer_reference_pressure_from_ak(np.array([0.0, 0.1, 0.5]), units="Pa") == 1.0


def test_vertical_validation_metrics_ignore_pairwise_nonfinite_values() -> None:
    metrics = metric_block(
        np.array([1.0, 2.0, np.nan, 4.0]),
        np.array([2.0, 2.0, 3.0, np.inf]),
    )

    assert metrics["count"] == 2
    assert metrics["rmse"] == pytest.approx(np.sqrt(0.5))
    assert metrics["mae"] == pytest.approx(0.5)
    assert metrics["max_abs"] == pytest.approx(1.0)
    assert metrics["bias"] == pytest.approx(0.5)


def test_pressure_metric_summary_reports_overall_and_per_level_errors() -> None:
    reference = xr.DataArray(
        np.array(
            [
                [[1.0, 2.0], [3.0, 4.0]],
                [[5.0, 6.0], [7.0, 8.0]],
            ]
        ),
        dims=("time", "plev", "values"),
        coords={"plev": [30000.0, 50000.0]},
    )
    candidate = reference.copy()
    candidate.loc[{"plev": 30000.0}] = candidate.sel(plev=30000.0) + 1.0

    summary = pressure_metric_summary(reference, candidate)

    assert summary["overall"]["count"] == 8
    assert summary["overall"]["rmse"] == pytest.approx(np.sqrt(0.5))
    assert summary["overall"]["max_abs"] == pytest.approx(1.0)
    assert summary["per_level"]["30000"]["rmse"] == pytest.approx(1.0)
    assert summary["per_level"]["50000"]["rmse"] == pytest.approx(0.0)


def test_pressure_vertical_interpolate_validates_dataset_variables_and_chunks() -> None:
    ds = xr.Dataset(
        {
            "t": xr.DataArray(np.ones((2, 3, 4)), dims=("time", "hybrid", "values")),
            "z": xr.DataArray(np.ones((2, 4)), dims=("time", "values")),
        }
    )

    with pytest.raises(ValueError, match="supported hybrid/model-level dimension"):
        vertical_interpolate(
            ds,
            target="pressure",
            levels=[850.0],
            variables=["z"],
        )

    with pytest.raises(ValueError, match="dimensions not present"):
        vertical_interpolate(
            ds,
            target="pressure",
            levels=[850.0],
            variables=["t"],
            chunks={"member": 1},
        )


def test_pressure_vertical_interpolate_returns_native_fullpos_output() -> None:
    data = xr.DataArray(
        np.array(
            [
                [[210.0, 220.0, 230.0, 240.0], [260.0, 270.0, 280.0, 290.0], [300.0, 310.0, 320.0, 330.0]],
                [[211.0, 221.0, 231.0, 241.0], [261.0, 271.0, 281.0, 291.0], [301.0, 311.0, 321.0, 331.0]],
            ]
        ),
        dims=("time", "hybrid", "values"),
        coords={"time": [0, 1]},
        attrs={"GRIB_pv": np.array([0.0, 20000.0, 60000.0, 100000.0, 0.0, 0.0, 0.0, 0.0])},
        name="q",
    )
    sp = xr.DataArray(np.full((2, 4), 90000.0), dims=("time", "values"), coords={"time": [0, 1]})

    out = vertical_interpolate(
        data,
        target="pressure",
        levels=[30000.0, 80000.0],
        chunks={"time": 1},
        surface_pressure=sp,
    )

    assert out.dims == ("time", "pressure", "values")
    assert out.shape == (2, 2, 4)
    assert out.attrs["vertical_backend"] == "FULLPOS"
    np.testing.assert_allclose(out["pressure"].values, [30000.0, 80000.0])
    assert np.isfinite(out.values).all()


def test_pressure_vertical_interpolate_dataset_uses_fullpos_wind_pair() -> None:
    coords = {"time": [0], "hybrid": [1, 2, 3], "values": [0, 1]}
    attrs = {"GRIB_pv": np.array([0.0, 20000.0, 60000.0, 100000.0, 0.0, 0.0, 0.0, 0.0])}
    ds = xr.Dataset(
        {
            "u": xr.DataArray(np.ones((1, 3, 2)), dims=("time", "hybrid", "values"), coords=coords, attrs=attrs),
            "v": xr.DataArray(np.full((1, 3, 2), 2.0), dims=("time", "hybrid", "values"), coords=coords, attrs=attrs),
            "q": xr.DataArray(np.full((1, 3, 2), 0.001), dims=("time", "hybrid", "values"), coords=coords, attrs=attrs),
        }
    )
    sp = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})

    out = vertical_interpolate(
        ds,
        target="pressure",
        levels=[30000.0, 80000.0],
        variables=["u", "v", "q"],
        chunks={"time": 1},
        surface_pressure=sp,
    )

    assert set(out.data_vars) == {"u", "v", "q"}
    assert out["u"].dims == ("time", "pressure", "values")
    np.testing.assert_allclose(out["u"].values, 1.0)
    np.testing.assert_allclose(out["v"].values, 2.0)
    np.testing.assert_allclose(out["q"].values, 0.001)


def test_native_column_pressure_matches_constant_pressure_levels() -> None:
    from fullpos import _vertical_native

    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    ps = np.array([90000.0, 90000.0], dtype=np.float64)
    levels = np.array([30000.0, 80000.0], dtype=np.float64)
    targets = np.asfortranarray(np.tile(levels, (2, 1)))
    values = np.asfortranarray(
        np.array([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]], dtype=np.float64)
    )

    by_level = _vertical_native.pressure_ppq(values, ak, bk, ps, levels)
    by_column = _vertical_native.column_pressure_ppq(values, ak, bk, ps, targets)

    np.testing.assert_allclose(by_column, by_level)


def test_native_hybrid_pressure_matches_column_pressure_targets() -> None:
    from fullpos import _vertical_native

    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    ps = np.array([90000.0, 95000.0], dtype=np.float64)
    target_ak = np.array([20000.0, 40000.0, 90000.0], dtype=np.float64)
    target_bk = np.zeros(3, dtype=np.float64)
    targets = np.asfortranarray(np.tile([30000.0, 65000.0], (2, 1)))
    values = np.asfortranarray(
        np.array([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]], dtype=np.float64)
    )

    by_column = _vertical_native.column_pressure_ppq(values, ak, bk, ps, targets)
    by_hybrid = _vertical_native.hybrid_pressure_ppq(values, ak, bk, ps, target_ak, target_bk)

    np.testing.assert_allclose(by_hybrid, by_column)


def test_native_theta_pressures_matches_gptet_formula_for_linear_profile() -> None:
    from fullpos import _vertical_native

    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    ps = np.array([90000.0, 95000.0], dtype=np.float64)
    pressure_full = np.array([10000.05, 40000.0, 80000.0], dtype=np.float64)
    kappa = 287.0596736665907 / 1004.7095955714683
    theta = np.array([[340.0, 300.0, 260.0], [350.0, 310.0, 270.0]], dtype=np.float64)
    temperature = np.asfortranarray(theta * (pressure_full[None, :] / 100000.0) ** kappa)
    levels = np.array([320.0, 330.0], dtype=np.float64)

    out = _vertical_native.theta_pressures(temperature, ak, bk, ps, levels)

    np.testing.assert_allclose(out[0], [25000.025, 17500.0375], rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(out[1], [32500.0125, 25000.025], rtol=1e-12, atol=1e-9)


def test_potential_temperature_vertical_interpolate_uses_native_fullpos_output() -> None:
    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    ps = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})
    pressure_full = np.array([10000.05, 40000.0, 80000.0], dtype=np.float64)
    kappa = 287.0596736665907 / 1004.7095955714683
    theta = np.array([[[260.0, 270.0], [300.0, 310.0], [340.0, 350.0]]], dtype=np.float64)
    temperature_values = theta * (pressure_full[None, :, None] / 100000.0) ** kappa
    attrs = {"GRIB_pv": np.concatenate([ak, bk])}
    temp = xr.DataArray(
        temperature_values,
        dims=("time", "hybrid", "values"),
        coords={"time": [0], "hybrid": [1, 2, 3], "values": [0, 1]},
        attrs=attrs,
        name="t",
    )
    q = xr.DataArray(
        np.ones((1, 3, 2), dtype=np.float64),
        dims=temp.dims,
        coords=temp.coords,
        attrs=attrs,
        name="q",
    )

    out = vertical_interpolate(
        q,
        target="potential_temperature",
        levels=[280.0, 320.0],
        surface_pressure=ps,
        temperature=temp,
    )

    assert out.dims == ("time", "potential_temperature", "values")
    assert out.shape == (1, 2, 2)
    assert out.attrs["vertical_backend"] == "FULLPOS"
    assert out.attrs["vertical_target"] == "potential_temperature"
    np.testing.assert_allclose(out["potential_temperature"].values, [280.0, 320.0])
    np.testing.assert_allclose(out.values, 1.0)


def test_potential_temperature_dataset_uses_temperature_variable_and_wind_pair() -> None:
    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    attrs = {"GRIB_pv": np.concatenate([ak, bk])}
    coords = {"time": [0], "hybrid": [1, 2, 3], "values": [0, 1]}
    temperature = xr.DataArray(
        np.full((1, 3, 2), 260.0),
        dims=("time", "hybrid", "values"),
        coords=coords,
        attrs=attrs,
        name="t",
    )
    ds = xr.Dataset(
        {
            "t": temperature,
            "u": xr.DataArray(np.ones((1, 3, 2)), dims=temperature.dims, coords=coords, attrs=attrs),
            "v": xr.DataArray(np.full((1, 3, 2), 2.0), dims=temperature.dims, coords=coords, attrs=attrs),
        }
    )
    sp = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})

    out = vertical_interpolate(
        ds,
        target="potential_temperature",
        levels=[280.0],
        variables=["u", "v"],
        surface_pressure=sp,
    )

    assert set(out.data_vars) == {"u", "v"}
    assert out["u"].dims == ("time", "potential_temperature", "values")
    np.testing.assert_allclose(out["u"].values, 1.0)
    np.testing.assert_allclose(out["v"].values, 2.0)


def test_native_temperature_pressures_returns_finite_targets() -> None:
    from fullpos import _vertical_native

    ak = np.array([0.0, 5000.0, 20000.0, 50000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(5, dtype=np.float64)
    ps = np.array([90000.0, 95000.0], dtype=np.float64)
    temperature = np.asfortranarray(
        np.array(
            [
                [215.0, 235.0, 255.0, 275.0],
                [220.0, 240.0, 260.0, 280.0],
            ],
            dtype=np.float64,
        )
    )

    out = _vertical_native.temperature_pressures(
        temperature,
        ak,
        bk,
        ps,
        np.array([250.0, 270.0], dtype=np.float64),
    )

    assert out.shape == (2, 2)
    assert np.isfinite(out).all()
    assert np.all(out > 0.0)


def test_temperature_vertical_interpolate_dataset_uses_temperature_variable_and_wind_pair() -> None:
    ak = np.array([0.0, 5000.0, 20000.0, 50000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(5, dtype=np.float64)
    attrs = {"GRIB_pv": np.concatenate([ak, bk])}
    coords = {"time": [0], "hybrid": [1, 2, 3, 4], "values": [0, 1]}
    temperature = xr.DataArray(
        np.array([[[215.0, 220.0], [235.0, 240.0], [255.0, 260.0], [275.0, 280.0]]]),
        dims=("time", "hybrid", "values"),
        coords=coords,
        attrs=attrs,
        name="t",
    )
    ds = xr.Dataset(
        {
            "t": temperature,
            "u": xr.DataArray(np.ones((1, 4, 2)), dims=temperature.dims, coords=coords, attrs=attrs),
            "v": xr.DataArray(np.full((1, 4, 2), 2.0), dims=temperature.dims, coords=coords, attrs=attrs),
            "q": xr.DataArray(np.full((1, 4, 2), 0.001), dims=temperature.dims, coords=coords, attrs=attrs),
        }
    )
    sp = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})

    out = vertical_interpolate(
        ds,
        target="temperature",
        levels=[250.0, 270.0],
        variables=["u", "v", "q"],
        surface_pressure=sp,
    )

    assert set(out.data_vars) == {"u", "v", "q"}
    assert out["u"].dims == ("time", "temperature", "values")
    np.testing.assert_allclose(out["temperature"].values, [250.0, 270.0])
    np.testing.assert_allclose(out["u"].values, 1.0)
    np.testing.assert_allclose(out["v"].values, 2.0)
    np.testing.assert_allclose(out["q"].values, 0.001)


def test_native_potential_vorticity_pressures_returns_finite_targets() -> None:
    from fullpos import _vertical_native

    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    ps = np.array([90000.0, 95000.0], dtype=np.float64)
    coriolis = np.array([1.0e-4, 1.2e-4], dtype=np.float64)
    pv = np.asfortranarray(
        np.array(
            [
                [1.0, 4.0, 8.0],
                [1.5, 4.5, 8.5],
            ],
            dtype=np.float64,
        )
    )

    out = _vertical_native.potential_vorticity_pressures(
        pv,
        ak,
        bk,
        ps,
        coriolis,
        np.array([2.0, 5.0], dtype=np.float64),
    )

    assert out.shape == (2, 2)
    assert np.isfinite(out).all()
    assert np.all(out > 0.0)


def test_native_gppvo_diagnoses_pv_and_theta() -> None:
    from fullpos import _vertical_native

    ak = np.array([10000.0, 30000.0, 70000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    ps = np.array([90000.0, 95000.0], dtype=np.float64)
    coriolis = np.array([1.0e-4, 1.2e-4], dtype=np.float64)
    u = np.asfortranarray(np.array([[10.0, 15.0, 20.0], [12.0, 17.0, 22.0]], dtype=np.float64))
    v = np.asfortranarray(np.array([[5.0, 7.0, 9.0], [6.0, 8.0, 10.0]], dtype=np.float64))
    t = np.asfortranarray(np.array([[260.0, 270.0, 280.0], [262.0, 272.0, 282.0]], dtype=np.float64))
    vort = np.asfortranarray(np.full((2, 3), 1.0e-5, dtype=np.float64))
    tm = np.asfortranarray(np.full((2, 3), 1.0e-6, dtype=np.float64))
    tl = np.asfortranarray(np.full((2, 3), -1.0e-6, dtype=np.float64))
    spm = np.array([0.01, 0.02], dtype=np.float64)
    spl = np.array([-0.01, -0.02], dtype=np.float64)
    kappa = np.asfortranarray(np.full((2, 3), 287.0596736665907 / 1004.7095955714683, dtype=np.float64))

    pv, theta = _vertical_native.diagnose_potential_vorticity(
        u,
        v,
        t,
        vort,
        tm,
        tl,
        spm,
        spl,
        kappa,
        ak,
        bk,
        ps,
        coriolis,
    )

    assert pv.shape == (2, 3)
    assert theta.shape == (2, 3)
    assert np.isfinite(pv).all()
    assert np.isfinite(theta).all()
    assert not np.allclose(pv, 0.0)
    assert np.all(theta > t)


def test_native_ectrans_vector_scalar_diagnostics_returns_finite_derivatives() -> None:
    from fullpos import _ectrans
    from fullpos.native import add_native_runtime_dir

    add_native_runtime_dir()
    pl = np.full(8, 16, dtype=np.int32)
    points = int(pl.sum())
    u = np.ascontiguousarray(np.ones((2, points), dtype=np.float64))
    v = np.ascontiguousarray(np.zeros((2, points), dtype=np.float64))
    scalars = np.ascontiguousarray(np.vstack([np.ones(points), np.ones(points)]), dtype=np.float64)

    vort, div, ns, ew = _ectrans.vector_scalar_diagnostics(u, v, scalars, pl, 3)

    assert vort.shape == (2, points)
    assert div.shape == (2, points)
    assert ns.shape == (2, points)
    assert ew.shape == (2, points)
    assert np.isfinite(vort).all()
    assert np.isfinite(div).all()
    assert np.isfinite(ns).all()
    assert np.isfinite(ew).all()
    np.testing.assert_allclose(ns, 0.0, atol=1.0e-12)
    np.testing.assert_allclose(ew, 0.0, atol=1.0e-12)


def test_native_gprcp_kappa_uses_fullpos_constants() -> None:
    from fullpos import _vertical_native

    q = np.zeros((3, 2), dtype=np.float64, order="F")
    kappa = _vertical_native.gprcp_kappa(q)

    assert kappa.shape == (3, 2)
    np.testing.assert_allclose(kappa, 0.2857142857142857, atol=1.0e-14)


def test_diagnose_potential_vorticity_auto_prepares_native_inputs() -> None:
    n = 4
    ak = np.array([10000.0, 30000.0, 70000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    attrs = {
        "GRIB_pv": np.concatenate([ak, bk]),
        "GRIB_gridType": "regular_gg",
        "GRIB_N": n,
        "GRIB_numberOfPoints": 8 * n * n,
    }
    coords = {
        "time": [0],
        "hybrid": [1, 2, 3],
        "latitude": gaussian_latitudes(2 * n),
        "longitude": regular_longitudes(4 * n),
    }
    shape = (1, 3, 2 * n, 4 * n)
    lat = np.deg2rad(coords["latitude"])[None, None, :, None]
    u = xr.DataArray(10.0 + np.zeros(shape), dims=("time", "hybrid", "latitude", "longitude"), coords=coords, attrs=attrs, name="u")
    v = xr.DataArray(2.0 * np.cos(lat) + np.zeros(shape), dims=u.dims, coords=coords, attrs=attrs, name="v")
    t = xr.DataArray(260.0 + 5.0 * np.sin(lat) + np.zeros(shape), dims=u.dims, coords=coords, attrs=attrs, name="t")
    q = xr.zeros_like(t)
    sp = xr.DataArray(
        np.full((1, 2 * n, 4 * n), 90000.0),
        dims=("time", "latitude", "longitude"),
        coords={"time": [0], "latitude": coords["latitude"], "longitude": coords["longitude"]},
    )

    out = diagnose_potential_vorticity(
        u=u,
        v=v,
        temperature=t,
        surface_pressure=sp,
        specific_humidity=q,
        source_grid="N4",
        chunks={"time": 1},
    )

    assert out.potential_vorticity.shape == shape
    assert out.potential_temperature.shape == shape
    assert out.potential_vorticity.attrs["vertical_native_path"] == "GPPVO"
    assert out.potential_temperature.attrs["vertical_backend"] == "FULLPOS"
    assert np.isfinite(out.potential_vorticity.values).all()
    assert np.isfinite(out.potential_temperature.values).all()


def test_diagnose_potential_vorticity_xarray_wrapper_uses_fullpos_gppvo() -> None:
    ak = np.array([0.0, 30000.0, 70000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    attrs = {"GRIB_pv": np.concatenate([ak, bk])}
    coords = {"time": [0], "hybrid": [1, 2, 3], "values": [0, 1]}
    dims = ("time", "hybrid", "values")
    u = xr.DataArray(np.array([[[10.0, 12.0], [15.0, 17.0], [20.0, 22.0]]]), dims=dims, coords=coords, attrs=attrs)
    v = xr.DataArray(np.array([[[5.0, 6.0], [7.0, 8.0], [9.0, 10.0]]]), dims=dims, coords=coords, attrs=attrs)
    t = xr.DataArray(np.array([[[260.0, 262.0], [270.0, 272.0], [280.0, 282.0]]]), dims=dims, coords=coords, attrs=attrs, name="t")
    vort = xr.full_like(t, 1.0e-5)
    tm = xr.full_like(t, 1.0e-6)
    tl = xr.full_like(t, -1.0e-6)
    kappa = xr.full_like(t, 287.0596736665907 / 1004.7095955714683)
    sp = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})
    coriolis = xr.DataArray(np.array([[1.0e-4, 1.2e-4]]), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})
    spm = xr.DataArray(np.array([[0.01, 0.02]]), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})
    spl = xr.DataArray(np.array([[-0.01, -0.02]]), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})

    out = diagnose_potential_vorticity(
        u=u,
        v=v,
        temperature=t,
        relative_vorticity=vort,
        temperature_meridional_gradient=tm,
        temperature_zonal_gradient=tl,
        surface_pressure=sp,
        coriolis=coriolis,
        surface_pressure_meridional_gradient=spm,
        surface_pressure_zonal_gradient=spl,
        kappa=kappa,
    )

    assert out.potential_vorticity.dims == dims
    assert out.potential_temperature.dims == dims
    assert out.potential_vorticity.attrs["vertical_backend"] == "FULLPOS"
    assert out.potential_temperature.attrs["vertical_native_path"] == "GPPVO"
    assert np.isfinite(out.potential_vorticity.values).all()
    assert np.isfinite(out.potential_temperature.values).all()


def test_potential_vorticity_vertical_interpolate_uses_native_fullpos_output() -> None:
    ak = np.array([0.0, 20000.0, 60000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    coords = {"time": [0], "hybrid": [1, 2, 3], "values": [0, 1]}
    sp = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})
    pv = xr.DataArray(
        np.array([[[1.0, 1.5], [4.0, 4.5], [8.0, 8.5]]], dtype=np.float64),
        dims=("time", "hybrid", "values"),
        coords=coords,
        name="pv",
    )
    q = xr.DataArray(
        np.full((1, 3, 2), 7.0, dtype=np.float64),
        dims=("time", "hybrid", "values"),
        coords=coords,
        attrs={"GRIB_pv": np.concatenate([ak, bk])},
        name="q",
    )
    ds = xr.Dataset({"q": q, "pv": pv})
    coriolis = xr.DataArray(np.array([[1.0e-4, 1.1e-4]], dtype=np.float64), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})

    out = vertical_interpolate(
        ds,
        target="potential_vorticity",
        levels=[2.0, 5.0],
        variables=["q"],
        surface_pressure=sp,
        potential_vorticity=ds["pv"],
        coriolis=coriolis,
    )

    assert set(out.data_vars) == {"q"}
    assert out["q"].dims == ("time", "potential_vorticity", "values")
    assert out["q"].shape == (1, 2, 2)
    assert out["q"].attrs["vertical_backend"] == "FULLPOS"
    assert out["q"].attrs["vertical_target"] == "potential_vorticity"
    np.testing.assert_allclose(out["q"]["potential_vorticity"].values, [2.0, 5.0])
    np.testing.assert_allclose(out["q"].values, 7.0)


def test_potential_vorticity_vertical_interpolate_auto_diagnoses_native_pv() -> None:
    n = 4
    ak = np.array([10000.0, 30000.0, 70000.0, 100000.0], dtype=np.float64)
    bk = np.zeros(4, dtype=np.float64)
    attrs = {
        "GRIB_pv": np.concatenate([ak, bk]),
        "GRIB_gridType": "regular_gg",
        "GRIB_N": n,
        "GRIB_numberOfPoints": 8 * n * n,
    }
    coords = {
        "time": [0],
        "hybrid": [1, 2, 3],
        "latitude": gaussian_latitudes(2 * n),
        "longitude": regular_longitudes(4 * n),
    }
    shape = (1, 3, 2 * n, 4 * n)
    lat = np.deg2rad(coords["latitude"])[None, None, :, None]
    u = xr.DataArray(10.0 + np.zeros(shape), dims=("time", "hybrid", "latitude", "longitude"), coords=coords, attrs=attrs, name="u")
    v = xr.DataArray(2.0 * np.cos(lat) + np.zeros(shape), dims=u.dims, coords=coords, attrs=attrs, name="v")
    t = xr.DataArray(260.0 + 5.0 * np.sin(lat) + np.zeros(shape), dims=u.dims, coords=coords, attrs=attrs, name="t")
    q = xr.zeros_like(t).rename("q")
    q.attrs = dict(attrs)
    sp = xr.DataArray(
        np.full((1, 2 * n, 4 * n), 90000.0),
        dims=("time", "latitude", "longitude"),
        coords={"time": [0], "latitude": coords["latitude"], "longitude": coords["longitude"]},
    )
    coriolis = xr.full_like(sp, 1.0e-4)
    ds = xr.Dataset({"u": u, "v": v, "t": t, "q": q})

    out = vertical_interpolate(
        ds,
        target="potential_vorticity",
        levels=[2.0e-6, 3.0e-6],
        variables=["q"],
        surface_pressure=sp,
        coriolis=coriolis,
        source_grid="N4",
        chunks={"time": 1},
    )

    assert set(out.data_vars) == {"q"}
    assert out["q"].dims == ("time", "potential_vorticity", "latitude", "longitude")
    assert out["q"].shape == (1, 2, 2 * n, 4 * n)
    assert out["q"].attrs["vertical_backend"] == "FULLPOS"
    assert out["q"].attrs["vertical_target"] == "potential_vorticity"
    np.testing.assert_allclose(out["q"]["potential_vorticity"].values, [2.0e-6, 3.0e-6])
    assert np.isfinite(out["q"].values).all()


def test_model_level_vertical_interpolate_uses_native_fullpos_output() -> None:
    data = xr.DataArray(
        np.array(
            [
                [[210.0, 220.0, 230.0, 240.0], [260.0, 270.0, 280.0, 290.0], [300.0, 310.0, 320.0, 330.0]],
                [[211.0, 221.0, 231.0, 241.0], [261.0, 271.0, 281.0, 291.0], [301.0, 311.0, 321.0, 331.0]],
            ]
        ),
        dims=("time", "hybrid", "values"),
        coords={"time": [0, 1]},
        attrs={"GRIB_pv": np.array([0.0, 20000.0, 60000.0, 100000.0, 0.0, 0.0, 0.0, 0.0])},
        name="q",
    )
    sp = xr.DataArray(np.full((2, 4), 90000.0), dims=("time", "values"), coords={"time": [0, 1]})

    out = vertical_interpolate(
        data,
        target="model_level",
        chunks={"time": 1},
        surface_pressure=sp,
    )

    assert out.dims == ("time", "model_level", "values")
    assert out.shape == data.shape
    assert out.attrs["vertical_backend"] == "FULLPOS"
    assert out.attrs["vertical_target"] == "model_level"
    np.testing.assert_array_equal(out["model_level"].values, [1, 2, 3])
    assert np.isfinite(out.values).all()


def test_model_level_vertical_interpolate_dataset_uses_fullpos_wind_pair() -> None:
    coords = {"time": [0], "hybrid": [1, 2, 3], "values": [0, 1]}
    attrs = {"GRIB_pv": np.array([0.0, 20000.0, 60000.0, 100000.0, 0.0, 0.0, 0.0, 0.0])}
    ds = xr.Dataset(
        {
            "u": xr.DataArray(np.ones((1, 3, 2)), dims=("time", "hybrid", "values"), coords=coords, attrs=attrs),
            "v": xr.DataArray(np.full((1, 3, 2), 2.0), dims=("time", "hybrid", "values"), coords=coords, attrs=attrs),
            "q": xr.DataArray(np.full((1, 3, 2), 0.001), dims=("time", "hybrid", "values"), coords=coords, attrs=attrs),
        }
    )
    sp = xr.DataArray(np.full((1, 2), 90000.0), dims=("time", "values"), coords={"time": [0], "values": [0, 1]})

    out = vertical_interpolate(
        ds,
        target="model_level",
        variables=["u", "v", "q"],
        chunks={"time": 1},
        surface_pressure=sp,
    )

    assert set(out.data_vars) == {"u", "v", "q"}
    assert out["u"].dims == ("time", "model_level", "values")
    np.testing.assert_allclose(out["u"].values, 1.0)
    np.testing.assert_allclose(out["v"].values, 2.0)
    np.testing.assert_allclose(out["q"].values, 0.001)
