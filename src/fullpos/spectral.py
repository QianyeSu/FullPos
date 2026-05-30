from __future__ import annotations

import numpy as np

from .errors import FullposBackendError
from .grids import GaussianGrid, parse_grid
from .native import add_native_runtime_dir, native_runtime_info


_NATIVE_FAILURE = None
_NATIVE_LAST_OK = False


def spectral_regrid_values(
    values: np.ndarray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    ntrunc: int | None = None,
    missing_policy: str = "error",
) -> np.ndarray:
    """Regrid a single native-layout field through ECTRANS spectral coefficients."""
    arr = _to_native_global(values, source_grid)
    out = spectral_regrid_batch(
        arr.reshape(1, -1),
        source_grid=source_grid,
        target_grid=target_grid,
        ntrunc=ntrunc,
        missing_policy=missing_policy,
    )
    if target_grid.is_reduced:
        return out.reshape(-1)
    return out.reshape(target_grid.nlat, target_grid.work_nlon)


def spectral_regrid_batch(
    values: np.ndarray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    ntrunc: int | None = None,
    missing_policy: str = "error",
) -> np.ndarray:
    """Regrid a ``(nfield, npoints)`` batch using the native ECTRANS backend."""
    ntrunc = _default_ntrunc(source_grid, target_grid, ntrunc)
    _validate_batch(values, source_grid, missing_policy=missing_policy)
    add_native_runtime_dir()
    _ectrans = _import_native_backend()

    global _NATIVE_FAILURE, _NATIVE_LAST_OK
    arr = np.asarray(values)
    try:
        out = _ectrans.regrid(
            np.ascontiguousarray(arr, dtype=np.float64),
            _native_pl(source_grid),
            _native_pl(target_grid),
            int(ntrunc),
        )
    except Exception as exc:
        _NATIVE_FAILURE = exc
        _NATIVE_LAST_OK = False
        raise FullposBackendError(
            "native FULLPOS/ECTRANS backend failed during spectral regridding"
        ) from exc

    _NATIVE_FAILURE = None
    _NATIVE_LAST_OK = True
    out = np.asarray(out, dtype=np.float32)
    if target_grid.is_reduced:
        return out
    return out.reshape((arr.shape[0], target_grid.nlat, target_grid.work_nlon))


def spectral_regrid_chunks(
    values: np.ndarray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    ntrunc: int | None = None,
    chunk_size: int | None = 64,
    missing_policy: str = "error",
) -> np.ndarray:
    """Chunk a flattened field batch and concatenate native spectral regridding results."""
    arr = np.asarray(values)
    if arr.ndim != 2:
        raise ValueError("chunked values must have shape (nfield, npoints)")
    if arr.shape[0] == 0:
        raise ValueError("cannot regrid an empty field batch")

    chunks = []
    for start, stop in _chunk_slices(arr.shape[0], chunk_size):
        native = spectral_regrid_batch(
            np.ascontiguousarray(arr[start:stop]),
            source_grid=source_grid,
            target_grid=target_grid,
            ntrunc=ntrunc,
            missing_policy=missing_policy,
        )
        if native.ndim > 2:
            native = native.reshape((native.shape[0], -1))
        chunks.append(native)
    return np.concatenate(chunks, axis=0)


def spectral_fit(
    values: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
    axis: int | tuple[int, ...] = -1,
    chunk_size: int | None = 64,
    missing_policy: str = "error",
) -> np.ndarray:
    """Fit grid-point field(s) to native ECTRANS spectral coefficients.

    The returned last dimension is the ECTRANS global real/imag coefficient
    layout with length ``(ntrunc + 1) * (ntrunc + 2)``. Leading dimensions are
    preserved.
    """
    source_grid = _ensure_grid(grid)
    ntrunc = _default_fit_ntrunc(source_grid, ntrunc)
    arr = np.asarray(values)
    expected_shape = _horizontal_shape(source_grid)
    moved, leading_shape = _move_horizontal_to_end(arr, expected_shape, axis)
    flat = moved.reshape((-1, source_grid.size))
    coeffs = spectral_fit_chunks(
        flat,
        grid=source_grid,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
    )
    return coeffs.reshape(leading_shape + (coeffs.shape[-1],))


def spectral_synthesis(
    coefficients: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
    axis: int = -1,
    chunk_size: int | None = 64,
) -> np.ndarray:
    """Synthesize grid-point field(s) from native ECTRANS spectral coefficients.

    ``coefficients`` must use the same global real/imag layout returned by
    :func:`spectral_fit`.
    """
    target_grid = _ensure_grid(grid)
    arr = np.asarray(coefficients)
    if arr.ndim == 0:
        raise ValueError("coefficients must have at least one dimension")
    coeff_axis = int(axis)
    if coeff_axis < 0:
        coeff_axis += arr.ndim
    if coeff_axis < 0 or coeff_axis >= arr.ndim:
        raise ValueError("axis is out of bounds")
    moved = np.moveaxis(arr, coeff_axis, -1)
    leading_shape = moved.shape[:-1]
    coeff_count = moved.shape[-1]
    ntrunc = _default_synthesis_ntrunc(target_grid, ntrunc, coeff_count)
    flat = moved.reshape((-1, coeff_count))
    values = spectral_synthesis_chunks(
        flat,
        grid=target_grid,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
    )
    target_shape = _horizontal_shape(target_grid)
    return values.reshape(leading_shape + target_shape)


def spectral_fit_values(
    values: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
    missing_policy: str = "error",
) -> np.ndarray:
    """Fit one native-layout field and return a 1D coefficient vector."""
    source_grid = _ensure_grid(grid)
    arr = _to_native_global(values, source_grid)
    out = spectral_fit_batch(
        arr.reshape(1, -1),
        grid=source_grid,
        ntrunc=ntrunc,
        missing_policy=missing_policy,
    )
    return out.reshape(-1)


def spectral_synthesis_values(
    coefficients: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
) -> np.ndarray:
    """Synthesize one coefficient vector to native-layout grid-point values."""
    target_grid = _ensure_grid(grid)
    arr = np.asarray(coefficients)
    if arr.ndim != 1:
        raise ValueError("single-field coefficients must be a 1D array")
    out = spectral_synthesis_batch(arr.reshape(1, -1), grid=target_grid, ntrunc=ntrunc)
    if target_grid.is_reduced:
        return out.reshape(-1)
    return out.reshape(target_grid.nlat, target_grid.work_nlon)


def spectral_fit_batch(
    values: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
    missing_policy: str = "error",
) -> np.ndarray:
    """Fit a flattened ``(nfield, npoints)`` batch to spectral coefficients."""
    source_grid = _ensure_grid(grid)
    ntrunc = _default_fit_ntrunc(source_grid, ntrunc)
    _validate_batch(values, source_grid, missing_policy=missing_policy)
    add_native_runtime_dir()
    _ectrans = _import_native_backend()

    global _NATIVE_FAILURE, _NATIVE_LAST_OK
    arr = np.asarray(values)
    try:
        out = _ectrans.fit(
            np.ascontiguousarray(arr, dtype=np.float64),
            _native_pl(source_grid),
            int(ntrunc),
        )
    except Exception as exc:
        _NATIVE_FAILURE = exc
        _NATIVE_LAST_OK = False
        raise FullposBackendError(
            "native FULLPOS/ECTRANS backend failed during spectral fitting"
        ) from exc

    _NATIVE_FAILURE = None
    _NATIVE_LAST_OK = True
    return np.asarray(out, dtype=np.float64)


def spectral_synthesis_batch(
    coefficients: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
) -> np.ndarray:
    """Synthesize a flattened coefficient batch to grid-point values."""
    target_grid = _ensure_grid(grid)
    arr = np.asarray(coefficients)
    if arr.ndim != 2:
        raise ValueError("batch coefficients must have shape (nfield, nspec2)")
    if arr.shape[0] == 0:
        raise ValueError("cannot synthesize an empty field batch")
    ntrunc = _default_synthesis_ntrunc(target_grid, ntrunc, arr.shape[1])
    add_native_runtime_dir()
    _ectrans = _import_native_backend()

    global _NATIVE_FAILURE, _NATIVE_LAST_OK
    try:
        out = _ectrans.synthesis(
            np.ascontiguousarray(arr, dtype=np.float64),
            _native_pl(target_grid),
            int(ntrunc),
        )
    except Exception as exc:
        _NATIVE_FAILURE = exc
        _NATIVE_LAST_OK = False
        raise FullposBackendError(
            "native FULLPOS/ECTRANS backend failed during spectral synthesis"
        ) from exc

    _NATIVE_FAILURE = None
    _NATIVE_LAST_OK = True
    if target_grid.is_reduced:
        return np.asarray(out, dtype=np.float64)
    return np.asarray(out, dtype=np.float64).reshape(
        (arr.shape[0], target_grid.nlat, target_grid.work_nlon)
    )


def spectral_wind_synthesis(
    vorticity: np.ndarray,
    divergence: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Synthesize U/V wind from vorticity/divergence spectral coefficients.

    This is a native ECTRANS ``nvordiv`` inverse-transform path. Inputs must be
    flattened coefficient batches with shape ``(nfield, nspec2)`` using the
    ECMWF/FULLPOS global real/imaginary coefficient layout accepted by
    :func:`spectral_synthesis_batch`. The returned arrays are ``(nfield,
    npoints)`` for reduced grids and ``(nfield, nlat, nlon)`` for regular
    Gaussian grids.
    """
    target_grid = _ensure_grid(grid)
    vor = np.asarray(vorticity)
    div = np.asarray(divergence)
    if vor.ndim != 2 or div.ndim != 2:
        raise ValueError("vorticity and divergence must have shape (nfield, nspec2)")
    if vor.shape != div.shape:
        raise ValueError("vorticity and divergence must have the same shape")
    if vor.shape[0] == 0:
        raise ValueError("cannot synthesize an empty vorticity/divergence batch")

    ntrunc = _default_coefficient_ntrunc(ntrunc, vor.shape[1])
    add_native_runtime_dir()
    _ectrans = _import_native_backend()

    global _NATIVE_FAILURE, _NATIVE_LAST_OK
    try:
        out = _ectrans.vordiv_synthesis(
            np.ascontiguousarray(vor, dtype=np.float64),
            np.ascontiguousarray(div, dtype=np.float64),
            _native_pl(target_grid),
            int(ntrunc),
        )
    except Exception as exc:
        _NATIVE_FAILURE = exc
        _NATIVE_LAST_OK = False
        raise FullposBackendError(
            "native FULLPOS/ECTRANS backend failed during wind synthesis"
        ) from exc

    _NATIVE_FAILURE = None
    _NATIVE_LAST_OK = True
    out = np.asarray(out, dtype=np.float64)
    nfield = vor.shape[0]
    u = out[:nfield]
    v = out[nfield:]
    if target_grid.is_reduced:
        return u, v
    shape = (nfield, target_grid.nlat, target_grid.work_nlon)
    return u.reshape(shape), v.reshape(shape)


def spectral_fit_chunks(
    values: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
    chunk_size: int | None = 64,
    missing_policy: str = "error",
) -> np.ndarray:
    """Chunk a flattened field batch and fit each chunk to coefficients."""
    source_grid = _ensure_grid(grid)
    arr = np.asarray(values)
    if arr.ndim != 2:
        raise ValueError("chunked values must have shape (nfield, npoints)")
    if arr.shape[0] == 0:
        raise ValueError("cannot fit an empty field batch")

    chunks = []
    for start, stop in _chunk_slices(arr.shape[0], chunk_size):
        chunks.append(
            spectral_fit_batch(
                np.ascontiguousarray(arr[start:stop]),
                grid=source_grid,
                ntrunc=ntrunc,
                missing_policy=missing_policy,
            )
        )
    return np.concatenate(chunks, axis=0)


def spectral_synthesis_chunks(
    coefficients: np.ndarray,
    *,
    grid: str | GaussianGrid,
    ntrunc: int | None = None,
    chunk_size: int | None = 64,
) -> np.ndarray:
    """Chunk a coefficient batch and synthesize each chunk to grid-point values."""
    target_grid = _ensure_grid(grid)
    arr = np.asarray(coefficients)
    if arr.ndim != 2:
        raise ValueError("chunked coefficients must have shape (nfield, nspec2)")
    if arr.shape[0] == 0:
        raise ValueError("cannot synthesize an empty coefficient batch")

    chunks = []
    for start, stop in _chunk_slices(arr.shape[0], chunk_size):
        native = spectral_synthesis_batch(
            np.ascontiguousarray(arr[start:stop]),
            grid=target_grid,
            ntrunc=ntrunc,
        )
        if native.ndim > 2:
            native = native.reshape((native.shape[0], -1))
        chunks.append(native)
    return np.concatenate(chunks, axis=0)


def native_backend_status() -> dict:
    """Return import/runtime status for the required native ECTRANS backend."""
    runtime = native_runtime_info()
    importable = False
    import_error = None
    try:
        add_native_runtime_dir()
        _import_native_backend()
        importable = True
    except Exception as exc:
        import_error = f"{type(exc).__name__}: {exc}"

    return {
        "backend": "native",
        "native_available": importable,
        "native_module": "fullpos._ectrans",
        **runtime,
        # Backward-compatible Windows-oriented aliases for existing diagnostics consumers.
        "native_dll_dir": runtime["native_runtime_dir"],
        "native_dll_dir_exists": runtime["native_runtime_dir_exists"],
        "required_native_dlls": runtime["required_native_libraries"],
        "required_native_dlls_present": runtime["required_native_libraries_present"],
        "last_native_ok": _NATIVE_LAST_OK,
        "last_native_failure": None
        if _NATIVE_FAILURE is None
        else f"{type(_NATIVE_FAILURE).__name__}: {_NATIVE_FAILURE}",
        "import_error": import_error,
    }


def _import_native_backend():
    global _NATIVE_FAILURE, _NATIVE_LAST_OK
    try:
        from . import _ectrans
    except Exception as exc:
        _NATIVE_FAILURE = exc
        _NATIVE_LAST_OK = False
        raise FullposBackendError(
            "native FULLPOS/ECTRANS backend is required but could not be imported"
        ) from exc
    return _ectrans


def _validate_batch(
    values: np.ndarray,
    source_grid: GaussianGrid,
    *,
    missing_policy: str,
) -> None:
    arr = np.asarray(values)
    if arr.ndim != 2:
        raise ValueError("batch values must have shape (nfield, npoints)")
    expected = source_grid.size
    if arr.shape[1] != expected:
        raise ValueError(f"batch fields have {arr.shape[1]} points, expected {expected}")
    policy = _normalize_missing_policy(missing_policy)
    if policy == "error":
        finite = np.isfinite(arr)
        if not bool(finite.all()):
            missing_count = int(arr.size - finite.sum())
            raise ValueError(
                "spectral regridding requires finite input; "
                f"found {missing_count} non-finite value(s) out of {arr.size}. "
                "Fields with GRIB bitmap or missing values, for example SST over land, "
                "need masked surface interpolation instead of global spectral interpolation. "
                "Pass missing_policy='ignore' only to preserve NaN-propagating behavior."
            )


def _normalize_missing_policy(missing_policy: str) -> str:
    policy = str(missing_policy).lower()
    if policy not in {"error", "ignore"}:
        raise ValueError("missing_policy must be 'error' or 'ignore'")
    return policy


def _chunk_slices(nfield: int, chunk_size: int | None):
    if chunk_size is None:
        yield 0, nfield
        return
    chunk = int(chunk_size)
    if chunk <= 0:
        raise ValueError("chunk_size must be a positive integer or None")
    for start in range(0, nfield, chunk):
        yield start, min(start + chunk, nfield)


def _native_pl(grid: GaussianGrid) -> np.ndarray:
    if grid.pl is None:
        return np.full(grid.nlat, grid.work_nlon, dtype=np.int32)
    return np.asarray(grid.pl, dtype=np.int32)


def _to_native_global(values: np.ndarray, grid: GaussianGrid) -> np.ndarray:
    arr = np.asarray(values)
    if grid.is_reduced:
        if arr.ndim != 1:
            raise ValueError(f"reduced {grid.name} input must be a 1D packed field")
        expected = grid.size
        if arr.size != expected:
            raise ValueError(f"packed field has {arr.size} points, expected {expected}")
        return np.ascontiguousarray(arr, dtype=np.float64)

    expected_shape = (grid.nlat, grid.work_nlon)
    if arr.shape != expected_shape:
        raise ValueError(f"regular {grid.name} input must have shape {expected_shape}, got {arr.shape}")
    return np.ascontiguousarray(arr.reshape(-1), dtype=np.float64)


def _ensure_grid(grid: str | GaussianGrid) -> GaussianGrid:
    if isinstance(grid, GaussianGrid):
        return grid
    return parse_grid(grid)


def _horizontal_shape(grid: GaussianGrid) -> tuple[int, ...]:
    if grid.is_reduced:
        return (grid.size,)
    return (grid.nlat, grid.work_nlon)


def _move_horizontal_to_end(
    arr: np.ndarray,
    expected_shape: tuple[int, ...],
    axis: int | tuple[int, ...],
) -> tuple[np.ndarray, tuple[int, ...]]:
    horizontal_ndim = len(expected_shape)
    if horizontal_ndim == 1:
        axes = (int(axis),)
    else:
        if not isinstance(axis, tuple):
            if int(axis) != -1:
                raise ValueError(
                    "regular Gaussian input uses two horizontal dimensions; "
                    "pass axis=(lat_axis, lon_axis)"
                )
            axes = tuple(range(arr.ndim - horizontal_ndim, arr.ndim))
        else:
            axes = tuple(int(v) for v in axis)
    axes = tuple(a + arr.ndim if a < 0 else a for a in axes)
    if len(axes) != horizontal_ndim or len(set(axes)) != horizontal_ndim:
        raise ValueError(f"axis must identify {horizontal_ndim} horizontal dimension(s)")
    if any(a < 0 or a >= arr.ndim for a in axes):
        raise ValueError("axis is out of bounds")
    actual_shape = tuple(arr.shape[a] for a in axes)
    if actual_shape != expected_shape:
        raise ValueError(f"horizontal shape must be {expected_shape}, got {actual_shape}")
    leading_axes = tuple(i for i in range(arr.ndim) if i not in axes)
    moved = np.transpose(arr, leading_axes + axes)
    leading_shape = moved.shape[: arr.ndim - horizontal_ndim]
    return moved, leading_shape


def _default_fit_ntrunc(grid: GaussianGrid, ntrunc: int | None) -> int:
    max_ntrunc = grid.n - 1
    if ntrunc is None:
        return max_ntrunc
    if ntrunc < 0 or ntrunc > max_ntrunc:
        raise ValueError(f"ntrunc must be between 0 and {max_ntrunc}, got {ntrunc}")
    return int(ntrunc)


def _default_synthesis_ntrunc(
    grid: GaussianGrid,
    ntrunc: int | None,
    coeff_count: int,
) -> int:
    if ntrunc is None:
        ntrunc = _infer_ntrunc_from_nspec2(coeff_count)
    ntrunc = _default_fit_ntrunc(grid, int(ntrunc))
    expected = _nspec2_from_ntrunc(ntrunc)
    if coeff_count != expected:
        raise ValueError(
            f"coefficients have {coeff_count} values per field, "
            f"but T{ntrunc} expects {expected}"
        )
    return ntrunc


def _default_coefficient_ntrunc(ntrunc: int | None, coeff_count: int) -> int:
    if ntrunc is None:
        ntrunc = _infer_ntrunc_from_nspec2(coeff_count)
    ntrunc = int(ntrunc)
    if ntrunc < 0:
        raise ValueError(f"ntrunc must be non-negative, got {ntrunc}")
    expected = _nspec2_from_ntrunc(ntrunc)
    if coeff_count != expected:
        raise ValueError(
            f"coefficients have {coeff_count} values per field, "
            f"but T{ntrunc} expects {expected}"
        )
    return ntrunc


def _infer_ntrunc_from_nspec2(coeff_count: int) -> int:
    if coeff_count <= 0:
        raise ValueError("coefficient count must be positive")
    root = int((np.sqrt(1 + 4 * coeff_count) - 3) / 2)
    if _nspec2_from_ntrunc(root) != coeff_count:
        raise ValueError("ntrunc is required when coefficient count is not triangular")
    return root


def _nspec2_from_ntrunc(ntrunc: int) -> int:
    return (int(ntrunc) + 1) * (int(ntrunc) + 2)


def _default_ntrunc(
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    ntrunc: int | None,
) -> int:
    max_ntrunc = min(source_grid.n - 1, target_grid.n - 1)
    if ntrunc is None:
        return max_ntrunc
    if ntrunc < 0 or ntrunc > max_ntrunc:
        raise ValueError(f"ntrunc must be between 0 and {max_ntrunc}, got {ntrunc}")
    return int(ntrunc)
