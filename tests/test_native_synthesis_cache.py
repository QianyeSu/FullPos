from __future__ import annotations

import numpy as np

from fullpos.grids import parse_grid
from fullpos.native import add_native_runtime_dir
from fullpos.spectral import spectral_synthesis_batch, spectral_wind_synthesis


def test_spectral_synthesis_batch_is_reproducible_across_calls() -> None:
    coeffs = np.arange(20.0, dtype=np.float64).reshape(1, -1)

    out1 = spectral_synthesis_batch(coeffs, grid="N4", ntrunc=3)
    out2 = spectral_synthesis_batch(coeffs, grid="N4", ntrunc=3)

    np.testing.assert_allclose(out1, out2, rtol=0.0, atol=0.0)


def test_native_vordiv_synthesis_is_reproducible_across_calls() -> None:
    grid = parse_grid("N4")
    pl = np.full(grid.nlat, grid.work_nlon, dtype=np.int32)
    coeffs = np.arange(20.0, dtype=np.float64).reshape(1, -1)

    add_native_runtime_dir()
    from fullpos import _ectrans

    out1 = _ectrans.vordiv_synthesis(coeffs, coeffs * 0.5, pl, 3)
    out2 = _ectrans.vordiv_synthesis(coeffs, coeffs * 0.5, pl, 3)

    assert out1.shape == (2, grid.size)
    np.testing.assert_allclose(out1, out2, rtol=0.0, atol=0.0)


def test_spectral_wind_synthesis_returns_u_v_fields() -> None:
    coeffs = np.arange(20.0, dtype=np.float64).reshape(1, -1)

    u1, v1 = spectral_wind_synthesis(coeffs, coeffs * 0.5, grid="N4", ntrunc=3)
    u2, v2 = spectral_wind_synthesis(coeffs, coeffs * 0.5, grid="N4", ntrunc=3)

    assert u1.shape == (1, 8, 16)
    assert v1.shape == (1, 8, 16)
    np.testing.assert_allclose(u1, u2, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(v1, v2, rtol=0.0, atol=0.0)


def test_spectral_synthesis_batch_is_reproducible_for_multi_field_inputs() -> None:
    coeffs = np.arange(40.0, dtype=np.float64).reshape(2, -1)

    out1 = spectral_synthesis_batch(coeffs, grid="N4", ntrunc=3)
    out2 = spectral_synthesis_batch(coeffs, grid="N4", ntrunc=3)

    np.testing.assert_allclose(out1, out2, rtol=0.0, atol=0.0)

