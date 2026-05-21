from __future__ import annotations

from importlib import import_module

import numpy as np

from ..grids import GaussianGrid, gaussian_latitudes, parse_grid
from ..native import add_native_runtime_dir


def _load_ectrans():
    return import_module("fullpos._ectrans")


def fpint4_kernel(
    pbuf,
    *,
    kl0,
    pwxx,
    pwxy=None,
    kbox=None,
    kder=None,
    kdmp=None,
    ldnrst=None,
    ldmask=None,
    ldml: bool = False,
    pundef: float = 1.0e20,
) -> np.ndarray:
    """Call the native OpenIFS/FULLPOS ``FPINT4`` kernel.

    This is a low-level kernel API. Inputs must already be FULLPOS-style
    interpolation buffers, addresses, and weights.
    """
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    pbuf_arr = np.asfortranarray(pbuf, dtype=np.float64)
    kl0_arr = np.asfortranarray(kl0, dtype=np.int32)
    kfproma = kl0_arr.shape[0]
    kfields = _infer_kfields(kder, ldnrst, ldmask)
    return _ectrans.fpint4_kernel(
        pbuf_arr,
        kl0_arr,
        np.asfortranarray(pwxx, dtype=np.float64),
        _default_pwxy(kfproma, pwxy),
        _default_int_array(kbox, kfproma, 1),
        _default_int_array(kder, kfields, 1),
        _default_int_array(kdmp, max(1, kfields), 0),
        _default_int_array(ldnrst, kfields, 0),
        _default_int_array(ldmask, kfields, 0),
        int(bool(ldml)),
        float(pundef),
    )


def fpint12_kernel(
    pbuf,
    *,
    kl0,
    pwxx,
    pwxy=None,
    kbox=None,
    kder=None,
    kdmp=None,
    ldnrst=None,
    ldmono=None,
    ldmask=None,
    ldml: bool = False,
    pundef: float = 1.0e20,
) -> np.ndarray:
    """Call the native OpenIFS/FULLPOS ``FPINT12`` kernel."""
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    pbuf_arr = np.asfortranarray(pbuf, dtype=np.float64)
    kl0_arr = np.asfortranarray(kl0, dtype=np.int32)
    kfproma = kl0_arr.shape[0]
    kfields = _infer_kfields(kder, ldnrst, ldmono, ldmask)
    return _ectrans.fpint12_kernel(
        pbuf_arr,
        kl0_arr,
        np.asfortranarray(pwxx, dtype=np.float64),
        _default_pwxy(kfproma, pwxy),
        _default_int_array(kbox, kfproma, 1),
        _default_int_array(kder, kfields, 1),
        _default_int_array(kdmp, max(1, kfields), 0),
        _default_int_array(ldnrst, kfields, 0),
        _default_int_array(ldmono, kfields, 0),
        _default_int_array(ldmask, kfields, 0),
        int(bool(ldml)),
        float(pundef),
    )


def fpavg_kernel(
    pbuf,
    *,
    ks0,
    pmask,
    kmask=None,
    ldmask=None,
    pundef: float = 1.0e20,
) -> np.ndarray:
    """Call the native OpenIFS/FULLPOS ``FPAVG`` kernel."""
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    kfields = _infer_kfields(kmask, ldmask)
    return _ectrans.fpavg_kernel(
        np.asfortranarray(pbuf, dtype=np.float64),
        np.asfortranarray(ks0, dtype=np.int32),
        np.asfortranarray(pmask, dtype=np.float64),
        _default_int_array(kmask, kfields, 1),
        _default_int_array(ldmask, kfields, 0),
        float(pundef),
    )


def fpnear_kernel(
    pbuf,
    *,
    ks0,
    pmask,
    kmask=None,
    ldmask=None,
    pundef: float = 1.0e20,
) -> np.ndarray:
    """Call the native OpenIFS/FULLPOS ``FPNEAR`` kernel."""
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    kfields = _infer_kfields(kmask, ldmask)
    return _ectrans.fpnear_kernel(
        np.asfortranarray(pbuf, dtype=np.float64),
        np.asfortranarray(ks0, dtype=np.int32),
        np.asfortranarray(pmask, dtype=np.float64),
        _default_int_array(kmask, kfields, 1),
        _default_int_array(ldmask, kfields, 0),
        float(pundef),
    )


def horizontal_regular_kernel(
    values,
    *,
    source_lats=None,
    source_lons=None,
    source_pl=None,
    source_grid: str | GaussianGrid | None = None,
    target_lats,
    target_lons,
    method: str = "bilinear",
    monotonic: bool = False,
) -> np.ndarray:
    """Call native FULLPOS ``SUHOW1/SUHOW2`` plus ``FPINT4``/``FPINT12``.

    The source can be either a rectangular regular-row field with
    ``source_lats`` and ``source_lons``, or a packed reduced Gaussian field with
    ``source_pl``/``source_grid``. Coordinates are in degrees.

    ``monotonic=True`` enables the native ``FPINT12`` monotonic clamp and is
    only valid together with ``method="quadratic12"``.
    """
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    field, nloen, src_lats = _prepare_horizontal_rows(
        values,
        source_lats=source_lats,
        source_lons=source_lons,
        source_pl=source_pl,
        source_grid=source_grid,
    )
    tgt_lats = np.asarray(target_lats, dtype=np.float64)
    tgt_lons = np.asarray(target_lons, dtype=np.float64)
    if tgt_lats.shape != tgt_lons.shape:
        raise ValueError("target_lats and target_lons must have matching shapes")
    normalized_method = _normalize_regular_method(method)
    if monotonic and normalized_method != "quadratic12":
        raise ValueError("monotonic=True is only supported with method='quadratic12'")

    flat = np.asfortranarray(field, dtype=np.float64)
    flat_tgt_lats = tgt_lats.reshape(-1)
    flat_tgt_lons = tgt_lons.reshape(-1)
    safe = _regular_stencil_safe_mask(flat_tgt_lats, src_lats)
    if np.all(safe):
        out = _ectrans.horizontal_regular_kernel(
            flat,
            np.asfortranarray(nloen, dtype=np.int32),
            np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lons), dtype=np.float64),
            normalized_method,
            int(bool(monotonic)),
        )
    else:
        out = np.empty(flat_tgt_lats.size, dtype=np.float64)
        if np.any(safe):
            out[safe] = _ectrans.horizontal_regular_kernel(
                flat,
                np.asfortranarray(nloen, dtype=np.int32),
                np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
                np.asfortranarray(np.deg2rad(flat_tgt_lats[safe]), dtype=np.float64),
                np.asfortranarray(np.deg2rad(flat_tgt_lons[safe]), dtype=np.float64),
                normalized_method,
                int(bool(monotonic)),
            )
        unsafe = ~safe
        out[unsafe] = _ectrans.horizontal_halo_kernel(
            np.asfortranarray(field, dtype=np.float64),
            np.asfortranarray(nloen, dtype=np.int32),
            np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lats[unsafe]), dtype=np.float64),
            np.asfortranarray(np.deg2rad(np.nextafter(flat_tgt_lons[unsafe], np.inf)), dtype=np.float64),
            "nearest",
            1,
        )
    return np.asarray(out, dtype=np.float64).reshape(tgt_lats.shape)


def horizontal_regular_kernel_batch(
    values,
    *,
    source_lats=None,
    source_lons=None,
    source_pl=None,
    source_grid: str | GaussianGrid | None = None,
    target_lats,
    target_lons,
    method: str = "bilinear",
    monotonic: bool = False,
) -> np.ndarray:
    """Call native FULLPOS regular interpolation for multiple fields at once."""
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    fields, nloen, src_lats = _prepare_horizontal_row_batch(
        values,
        source_lats=source_lats,
        source_lons=source_lons,
        source_pl=source_pl,
        source_grid=source_grid,
    )
    tgt_lats = np.asarray(target_lats, dtype=np.float64)
    tgt_lons = np.asarray(target_lons, dtype=np.float64)
    if tgt_lats.shape != tgt_lons.shape:
        raise ValueError("target_lats and target_lons must have matching shapes")
    normalized_method = _normalize_regular_method(method)
    if monotonic and normalized_method != "quadratic12":
        raise ValueError("monotonic=True is only supported with method='quadratic12'")

    flat = np.asfortranarray(fields.T, dtype=np.float64)
    flat_tgt_lats = tgt_lats.reshape(-1)
    flat_tgt_lons = tgt_lons.reshape(-1)
    safe = _regular_stencil_safe_mask(flat_tgt_lats, src_lats)
    if np.all(safe):
        out = _ectrans.horizontal_regular_kernel_batch(
            flat,
            np.asfortranarray(nloen, dtype=np.int32),
            np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lons), dtype=np.float64),
            normalized_method,
            int(bool(monotonic)),
        )
        return np.asarray(out, dtype=np.float64).T.reshape((fields.shape[0],) + tgt_lats.shape)

    out = np.empty((fields.shape[0], flat_tgt_lats.size), dtype=np.float64)
    if np.any(safe):
        safe_out = _ectrans.horizontal_regular_kernel_batch(
            flat,
            np.asfortranarray(nloen, dtype=np.int32),
            np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lats[safe]), dtype=np.float64),
            np.asfortranarray(np.deg2rad(flat_tgt_lons[safe]), dtype=np.float64),
            normalized_method,
            int(bool(monotonic)),
        )
        out[:, safe] = np.asarray(safe_out, dtype=np.float64).T
    unsafe = ~safe
    unsafe_out = _ectrans.horizontal_halo_kernel_batch(
        flat,
        np.asfortranarray(nloen, dtype=np.int32),
        np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
        np.asfortranarray(np.deg2rad(flat_tgt_lats[unsafe]), dtype=np.float64),
        np.asfortranarray(np.deg2rad(np.nextafter(flat_tgt_lons[unsafe], np.inf)), dtype=np.float64),
        "nearest",
        1,
    )
    out[:, unsafe] = np.asarray(unsafe_out, dtype=np.float64).T
    return out.reshape((fields.shape[0],) + tgt_lats.shape)


def horizontal_halo_kernel(
    values,
    *,
    source_lats=None,
    source_lons=None,
    source_pl=None,
    source_grid: str | GaussianGrid | None = None,
    target_lats,
    target_lons,
    method: str = "nearest",
    kslwide: int = 1,
) -> np.ndarray:
    """Call native FULLPOS ``SUHOX1`` plus ``FPNEAR``/``FPAVG``.

    This is the unmasked halo interpolation path. ``kslwide`` is the FULLPOS
    halo half-width; the default ``1`` uses a 2x2 halo.
    """
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    field, nloen, src_lats = _prepare_horizontal_rows(
        values,
        source_lats=source_lats,
        source_lons=source_lons,
        source_pl=source_pl,
        source_grid=source_grid,
    )
    tgt_lats = np.asarray(target_lats, dtype=np.float64)
    tgt_lons = np.nextafter(np.asarray(target_lons, dtype=np.float64), np.inf)
    if tgt_lats.shape != tgt_lons.shape:
        raise ValueError("target_lats and target_lons must have matching shapes")
    out = _ectrans.horizontal_halo_kernel(
        np.asfortranarray(field, dtype=np.float64),
        np.asfortranarray(nloen, dtype=np.int32),
        np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
        np.asfortranarray(np.deg2rad(tgt_lats.reshape(-1)), dtype=np.float64),
        np.asfortranarray(np.deg2rad(tgt_lons.reshape(-1)), dtype=np.float64),
        _normalize_halo_method(method),
        int(kslwide),
    )
    return np.asarray(out, dtype=np.float64).reshape(tgt_lats.shape)


def horizontal_halo_kernel_batch(
    values,
    *,
    source_lats=None,
    source_lons=None,
    source_pl=None,
    source_grid: str | GaussianGrid | None = None,
    target_lats,
    target_lons,
    method: str = "nearest",
    kslwide: int = 1,
) -> np.ndarray:
    """Call native FULLPOS halo interpolation for multiple fields at once."""
    add_native_runtime_dir()
    _ectrans = _load_ectrans()

    fields, nloen, src_lats = _prepare_horizontal_row_batch(
        values,
        source_lats=source_lats,
        source_lons=source_lons,
        source_pl=source_pl,
        source_grid=source_grid,
    )
    tgt_lats = np.asarray(target_lats, dtype=np.float64)
    tgt_lons = np.nextafter(np.asarray(target_lons, dtype=np.float64), np.inf)
    if tgt_lats.shape != tgt_lons.shape:
        raise ValueError("target_lats and target_lons must have matching shapes")
    flat = np.asfortranarray(fields.T, dtype=np.float64)
    out = _ectrans.horizontal_halo_kernel_batch(
        flat,
        np.asfortranarray(nloen, dtype=np.int32),
        np.asfortranarray(np.deg2rad(src_lats), dtype=np.float64),
        np.asfortranarray(np.deg2rad(tgt_lats.reshape(-1)), dtype=np.float64),
        np.asfortranarray(np.deg2rad(tgt_lons.reshape(-1)), dtype=np.float64),
        _normalize_halo_method(method),
        int(kslwide),
    )
    return np.asarray(out, dtype=np.float64).T.reshape((fields.shape[0],) + tgt_lats.shape)


def _prepare_horizontal_rows(
    values,
    *,
    source_lats,
    source_lons,
    source_pl,
    source_grid: str | GaussianGrid | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=np.float64)
    nloen = _resolve_nloen(source_pl=source_pl, source_grid=source_grid)
    if nloen is None:
        src_lats = np.asarray(source_lats, dtype=np.float64)
        src_lons = np.asarray(source_lons, dtype=np.float64)
        if src_lats.ndim != 1 or src_lons.ndim != 1:
            raise ValueError("source_lats and source_lons must be 1D arrays")
        if arr.shape != (src_lats.size, src_lons.size):
            raise ValueError(
                "regular values must have shape "
                f"{(src_lats.size, src_lons.size)}, got {arr.shape}"
            )
        return (
            np.ascontiguousarray(arr.reshape(-1), dtype=np.float64),
            np.full(src_lats.size, src_lons.size, dtype=np.int32),
            src_lats,
        )

    if arr.ndim != 1:
        raise ValueError("packed reduced values must be a 1D array")
    if arr.size != int(nloen.sum()):
        raise ValueError(f"packed field has {arr.size} points, expected {int(nloen.sum())}")
    src_lats = (
        gaussian_latitudes(nloen.size)
        if source_lats is None
        else np.asarray(source_lats, dtype=np.float64)
    )
    if src_lats.ndim != 1 or src_lats.size != nloen.size:
        raise ValueError("source_lats must have one latitude per source row")
    return np.ascontiguousarray(arr, dtype=np.float64), nloen, src_lats


def _prepare_horizontal_row_batch(
    values,
    *,
    source_lats,
    source_lons,
    source_pl,
    source_grid: str | GaussianGrid | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    arr = np.asarray(values, dtype=np.float64)
    nloen = _resolve_nloen(source_pl=source_pl, source_grid=source_grid)
    if nloen is None:
        src_lats = np.asarray(source_lats, dtype=np.float64)
        src_lons = np.asarray(source_lons, dtype=np.float64)
        if src_lats.ndim != 1 or src_lons.ndim != 1:
            raise ValueError("source_lats and source_lons must be 1D arrays")
        if arr.ndim == 2:
            arr = arr.reshape((1,) + arr.shape)
        if arr.ndim != 3 or arr.shape[-2:] != (src_lats.size, src_lons.size):
            raise ValueError(
                "regular batch values must have shape "
                f"(nfields, {src_lats.size}, {src_lons.size}), got {arr.shape}"
            )
        fields = arr.reshape((arr.shape[0], src_lats.size * src_lons.size))
        return np.ascontiguousarray(fields, dtype=np.float64), np.full(src_lats.size, src_lons.size, dtype=np.int32), src_lats

    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.ndim != 2:
        raise ValueError("packed reduced batch values must have shape (nfields, nsrc_points)")
    if arr.shape[1] != int(nloen.sum()):
        raise ValueError(f"packed batch has {arr.shape[1]} points per field, expected {int(nloen.sum())}")
    src_lats = (
        gaussian_latitudes(nloen.size)
        if source_lats is None
        else np.asarray(source_lats, dtype=np.float64)
    )
    if src_lats.ndim != 1 or src_lats.size != nloen.size:
        raise ValueError("source_lats must have one latitude per source row")
    return np.ascontiguousarray(arr, dtype=np.float64), nloen, src_lats


def _regular_stencil_safe_mask(target_lats: np.ndarray, source_lats: np.ndarray) -> np.ndarray:
    if source_lats.size < 5:
        return np.zeros(target_lats.shape, dtype=bool)
    north = float(source_lats[2])
    south = float(source_lats[-3])
    return (target_lats <= north) & (target_lats >= south)


def _resolve_nloen(*, source_pl, source_grid: str | GaussianGrid | None) -> np.ndarray | None:
    if source_pl is not None:
        nloen = np.asarray(source_pl, dtype=np.int32).reshape(-1)
        if nloen.size == 0 or np.any(nloen <= 0):
            raise ValueError("source_pl must contain positive row lengths")
        return nloen
    if source_grid is None:
        return None
    grid = source_grid if isinstance(source_grid, GaussianGrid) else parse_grid(source_grid)
    if grid.pl is None:
        return None
    return np.asarray(grid.pl, dtype=np.int32)


def _default_pwxy(kfproma: int, pwxy) -> np.ndarray:
    if pwxy is None:
        return np.asfortranarray(np.zeros((kfproma, 12), dtype=np.float64))
    return np.asfortranarray(pwxy, dtype=np.float64)


def _default_int_array(values, size: int, fill: int) -> np.ndarray:
    if values is None:
        return np.asfortranarray(np.full(size, fill, dtype=np.int32))
    return np.asfortranarray(values, dtype=np.int32)


def _infer_kfields(*candidates) -> int:
    for candidate in candidates:
        if candidate is not None:
            return int(np.asarray(candidate).shape[0])
    return 1


def _normalize_regular_method(method: str) -> str:
    normalized = str(method).lower().replace("-", "_")
    aliases = {
        "linear": "bilinear",
        "4": "bilinear",
        "fpint4": "bilinear",
        "12": "quadratic12",
        "quadratic": "quadratic12",
        "fpint12": "quadratic12",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"bilinear", "quadratic12"}:
        raise ValueError("method must be 'bilinear'/'fpint4' or 'quadratic12'/'fpint12'")
    return normalized


def _normalize_halo_method(method: str) -> str:
    normalized = str(method).lower().replace("-", "_")
    aliases = {
        "nearest_neighbour": "nearest",
        "nearest_neighbor": "nearest",
        "fpnear": "nearest",
        "avg": "average",
        "fpavg": "average",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"nearest", "average"}:
        raise ValueError("method must be 'nearest'/'fpnear' or 'average'/'fpavg'")
    return normalized
