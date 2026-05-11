from __future__ import annotations

import numpy as np

from fullpos import generic_spectral_filter
from fullpos.filters import (
    SpectralFilter,
    expand_profile_to_coefficients,
    fpfilter_profile,
    gaussian_filter_profile,
    load_filter_matrix,
    load_spectral_filter,
    low_pass_filter_profile,
    save_filter_matrix,
    save_spectral_filter,
)


def test_gaussian_filter_profile_matches_fpfilter_formula() -> None:
    profile = gaussian_filter_profile(4, cutoff=3, exponent=2.0)
    n = np.arange(5, dtype=np.float64)
    expected = np.exp(-0.5 * 2.0 * n**2 / 3.0**2)

    np.testing.assert_allclose(profile, expected)


def test_low_pass_filter_profile_matches_fpfilter_formula() -> None:
    profile = low_pass_filter_profile(5, cutoff=3, selectivity=2.0, exponent=0.5)
    n = np.arange(6, dtype=np.float64)
    expected = 0.5 * (1.0 - np.tanh(np.exp(-0.5 * 2.0) * (n - 3.0)))

    np.testing.assert_allclose(profile, expected)


def test_fpfilter_profile_accepts_fullpos_aliases() -> None:
    np.testing.assert_allclose(
        fpfilter_profile(4, kind="classic"),
        gaussian_filter_profile(4),
    )
    np.testing.assert_allclose(
        fpfilter_profile(4, kind="THX", cutoff=2),
        low_pass_filter_profile(4, cutoff=2),
    )


def test_expand_profile_to_coefficients_uses_total_wave_number() -> None:
    weights = expand_profile_to_coefficients(np.array([10.0, 20.0, 30.0]))

    np.testing.assert_allclose(weights, [10.0, 10.0, 20.0, 20.0, 30.0, 30.0, 20.0, 20.0, 30.0, 30.0, 30.0, 30.0])


def test_generic_spectral_filter_low_pass_reduces_high_frequency_variance() -> None:
    lon = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    values = np.tile(np.sin(6.0 * lon), (8, 1)).astype(np.float32)

    out = generic_spectral_filter(
        values,
        grid="N4",
        filter_kind="low_pass",
        ntrunc=3,
        cutoff=1,
        axis=(-2, -1),
        chunk_size=1,
    )

    assert out.shape == values.shape
    assert float(np.std(out)) < float(np.std(values))


def test_spectral_filter_object_applies_to_coefficients() -> None:
    filt = SpectralFilter.low_pass(2, cutoff=1)
    coeffs = np.ones((2, filt.nspec2), dtype=np.float64)

    out = filt.apply_coefficients(coeffs)

    assert out.shape == coeffs.shape
    np.testing.assert_allclose(out[0], filt.coefficient_weights)
    assert filt.to_metadata()["layout"] == "ectrans_global_real_imag_by_total_wave_number"


def test_spectral_filter_object_applies_to_grid_values() -> None:
    lon = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    values = np.tile(np.sin(6.0 * lon), (8, 1)).astype(np.float32)
    filt = SpectralFilter.low_pass(3, cutoff=1)

    out = filt.apply(values, grid="N4", axis=(-2, -1), chunk_size=1)

    assert out.shape == values.shape
    assert float(np.std(out)) < float(np.std(values))


def test_generic_spectral_filter_accepts_filter_object() -> None:
    values = np.ones((8, 16), dtype=np.float32)
    filt = SpectralFilter.gaussian(3)

    out = generic_spectral_filter(values, grid="N4", filter_kind=filt, axis=(-2, -1), chunk_size=1)

    assert out.shape == values.shape
    np.testing.assert_allclose(out, values, atol=1.0e-5)


def test_spectral_filter_roundtrip_io(tmp_path) -> None:
    path = tmp_path / "filter.npz"
    filt = SpectralFilter.low_pass(3, cutoff=2, selectivity=4.0, exponent=0.5)

    save_spectral_filter(path, filt)
    loaded = load_spectral_filter(path)

    assert loaded.kind == "low_pass"
    assert loaded.cutoff == 2
    np.testing.assert_allclose(loaded.profile, filt.profile)
    np.testing.assert_allclose(loaded.coefficient_weights, filt.coefficient_weights)


def test_filter_matrix_roundtrip_io(tmp_path) -> None:
    path = tmp_path / "matrix.npz"
    profile = gaussian_filter_profile(3, exponent=2.0)

    save_filter_matrix(path, profile)
    loaded = load_filter_matrix(path)

    np.testing.assert_allclose(loaded["profile"], profile)
    np.testing.assert_allclose(
        loaded["coefficient_weights"],
        expand_profile_to_coefficients(profile),
    )
    assert loaded["metadata"]["kind"] == "diagonal_profile_matrix"
