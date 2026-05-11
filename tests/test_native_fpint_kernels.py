from __future__ import annotations

import numpy as np

from fullpos.interpolation.kernels import (
    fpavg_kernel,
    fpint4_kernel,
    fpint12_kernel,
    fpnear_kernel,
)


def test_fpint4_kernel_calls_native_fullpos_weights() -> None:
    pbuf = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    kl0 = np.array([[0, 1, 3, 0]], dtype=np.int32)
    pwxx = np.array([[0.25, 0.25, 0.25, 0.25]])

    out = fpint4_kernel(pbuf, kl0=kl0, pwxx=pwxx)

    np.testing.assert_allclose(out, [[35.0]])


def test_fpint12_kernel_calls_native_fullpos_weights() -> None:
    pbuf = np.arange(1.0, 20.0)
    kl0 = np.array([[0, 1, 3, 5]], dtype=np.int32)
    pwxx = np.zeros((1, 12), dtype=np.float64)
    pwxx[0, 0] = 1.0

    out = fpint12_kernel(pbuf, kl0=kl0, pwxx=pwxx)

    np.testing.assert_allclose(out, [[2.0]])


def test_fpavg_kernel_calls_native_fullpos_halo_average() -> None:
    pbuf = np.array([10.0, 20.0, 30.0, 40.0])
    ks0 = np.array([[1, 3]], dtype=np.int32)
    pmask = np.array([[1.0]])

    out = fpavg_kernel(pbuf, ks0=ks0, pmask=pmask)

    np.testing.assert_allclose(out, [[25.0]])


def test_fpnear_kernel_calls_native_fullpos_nearest_valid_halo_point() -> None:
    pundef = 1.0e20
    pbuf = np.array([pundef, 20.0, 30.0, 40.0])
    ks0 = np.array([[1, 3]], dtype=np.int32)
    pmask = np.array([[1.0]])

    out = fpnear_kernel(pbuf, ks0=ks0, pmask=pmask, pundef=pundef)

    np.testing.assert_allclose(out, [[20.0]])


def test_fpint4_kernel_matches_openifs_formula_for_multiple_points() -> None:
    rng = np.random.default_rng(42)
    pbuf = rng.normal(size=80)
    kl0 = np.array(
        [
            [0, 5, 20, 0],
            [0, 7, 25, 0],
            [0, 9, 30, 0],
        ],
        dtype=np.int32,
    )
    pwxx = rng.random((3, 4))
    pwxx /= pwxx.sum(axis=1, keepdims=True)

    out = fpint4_kernel(pbuf, kl0=kl0, pwxx=pwxx)
    expected = np.empty((3, 1))
    for j in range(3):
        idx = [kl0[j, 1] + 1, kl0[j, 1] + 2, kl0[j, 2] + 1, kl0[j, 2] + 2]
        expected[j, 0] = np.dot(pwxx[j], pbuf[np.array(idx) - 1])

    np.testing.assert_allclose(out, expected, rtol=0.0, atol=1.0e-14)


def test_fpint12_kernel_matches_openifs_formula_for_multiple_points() -> None:
    rng = np.random.default_rng(123)
    pbuf = rng.normal(size=100)
    kl0 = np.array(
        [
            [10, 30, 50, 70],
            [12, 32, 52, 72],
        ],
        dtype=np.int32,
    )
    pwxx = rng.normal(size=(2, 12))
    out = fpint12_kernel(pbuf, kl0=kl0, pwxx=pwxx)
    expected = np.empty((2, 1))
    offsets = [
        (1, 1),
        (1, 2),
        (2, 1),
        (2, 2),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 3),
        (2, 0),
        (2, 3),
        (3, 1),
        (3, 2),
    ]
    for j in range(2):
        values = np.array([pbuf[kl0[j, row] + col - 1] for row, col in offsets])
        expected[j, 0] = np.dot(pwxx[j], values)

    np.testing.assert_allclose(out, expected, rtol=0.0, atol=1.0e-14)


def test_fpavg_and_fpnear_kernels_match_openifs_halo_logic() -> None:
    pundef = 1.0e20
    pbuf = np.array([10.0, pundef, 30.0, 40.0, 50.0, 60.0])
    ks0 = np.array([[1, 4]], dtype=np.int32)
    pmask = np.array([[1.0]])

    avg = fpavg_kernel(pbuf, ks0=ks0, pmask=pmask, pundef=pundef)
    near = fpnear_kernel(pbuf, ks0=ks0, pmask=pmask, pundef=pundef)

    np.testing.assert_allclose(avg, [[(10.0 + 40.0 + 50.0) / 3.0]])
    np.testing.assert_allclose(near, [[10.0]])
