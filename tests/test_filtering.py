from __future__ import annotations

import numpy as np
import xarray as xr

from fullpos import spectral_filter, spectral_fit, spectral_synthesis


def test_spectral_filter_reduces_high_frequency_variance() -> None:
    lon = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    values = np.tile(np.sin(6.0 * lon), (8, 1)).astype(np.float32)
    obj = xr.DataArray(values, dims=("latitude", "longitude"))

    out = spectral_filter(obj, grid="N4", ntrunc=1, chunk_size=1)

    assert out.shape == obj.shape
    assert float(np.std(out.values)) < float(np.std(obj.values))
    assert "fullpos regrid" in out.attrs["history"]
    assert "ntrunc=1" in out.attrs["history"]


def test_spectral_filter_supports_numpy_values() -> None:
    lon = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    field = np.tile(np.sin(6.0 * lon), (8, 1))
    values = np.stack([field, field]).astype(np.float32)

    out = spectral_filter(values, grid="N4", ntrunc=1, axis=(-2, -1), chunk_size=1)

    assert out.shape == values.shape
    assert float(np.std(out)) < float(np.std(values))


def test_spectral_filter_dataset_skips_non_horizontal_variables() -> None:
    lon = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    values = np.tile(np.sin(6.0 * lon), (8, 1)).astype(np.float32)
    ds = xr.Dataset(
        data_vars={
            "field": (("latitude", "longitude"), values),
            "scalar": ((), np.float32(1.0)),
        }
    )

    out = spectral_filter(ds, grid="N4", ntrunc=1, chunk_size=1)

    assert out["field"].shape == values.shape
    assert out["scalar"].shape == ()
    assert "spectral_filter" in out.attrs["history"]


def test_spectral_fit_synthesis_roundtrip_regular_grid() -> None:
    field = np.full((8, 16), 3.0)

    coeffs = spectral_fit(field, grid="N4", ntrunc=3, chunk_size=1)
    out = spectral_synthesis(coeffs, grid="N4", ntrunc=3, chunk_size=1)

    assert coeffs.shape == (20,)
    assert out.shape == field.shape
    np.testing.assert_allclose(out, field, rtol=2e-5, atol=2e-5)


def test_spectral_fit_synthesis_supports_leading_dimensions() -> None:
    field = np.full((8, 16), 2.0)
    values = np.stack([field, 2.0 * field]).reshape(1, 2, 8, 16)

    coeffs = spectral_fit(values, grid="N4", ntrunc=3, axis=(-2, -1), chunk_size=1)
    out = spectral_synthesis(coeffs, grid="N4", ntrunc=3, chunk_size=1)

    assert coeffs.shape == (1, 2, 20)
    assert out.shape == values.shape
    np.testing.assert_allclose(out, values, rtol=2e-5, atol=2e-5)


def test_spectral_fit_synthesis_matches_spectral_filter() -> None:
    lon = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    field = np.tile(np.sin(6.0 * lon), (8, 1))
    values = np.stack([field, field]).astype(np.float32)

    coeffs = spectral_fit(values, grid="N4", ntrunc=1, axis=(-2, -1), chunk_size=1)
    out = spectral_synthesis(coeffs, grid="N4", ntrunc=1, chunk_size=1)
    filtered = spectral_filter(values, grid="N4", ntrunc=1, axis=(-2, -1), chunk_size=1)

    np.testing.assert_allclose(out, filtered, rtol=2e-5, atol=2e-5)
