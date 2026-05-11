from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from ..metadata import append_history
from ..native import add_native_runtime_dir


@dataclass(frozen=True)
class PressureRequest:
    """Validated hybrid-to-pressure request prepared for a native backend."""

    levels: np.ndarray
    hybrid_dim: str
    ak: np.ndarray
    bk: np.ndarray
    surface_pressure: xr.DataArray
    variables: tuple[str, ...] | None


def interpolate_to_pressure(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    keep_attrs: bool = True,
):
    """Interpolate hybrid model-level fields to pressure levels with FULLPOS.

    The numerical path uses native ECMWF/OpenIFS FULLPOS vertical routines
    compiled into ``fullpos._vertical_native``. Python only validates xarray
    metadata, chunks leading dimensions, and restores coordinates/attributes.
    """
    request = prepare_pressure_request(
        values,
        levels=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
    )
    if isinstance(values, xr.Dataset):
        selected = request.variables or tuple(values.data_vars)
        outputs: dict[str, xr.DataArray] = {}
        used: set[str] = set()
        if "u" in selected and "v" in selected:
            u_out, v_out = _interpolate_wind_pair(
                values["u"],
                values["v"],
                request=request,
                chunks=chunks,
                keep_attrs=keep_attrs,
            )
            outputs["u"] = u_out
            outputs["v"] = v_out
            used.update(("u", "v"))
        for name in selected:
            if name in used:
                continue
            outputs[name] = _interpolate_data_array(
                values[name],
                request=request,
                chunks=chunks,
                variable_name=name,
                keep_attrs=keep_attrs,
            )
        out = xr.Dataset(outputs, attrs=dict(values.attrs) if keep_attrs else {})
        if keep_attrs:
            out.attrs = append_history(
                out.attrs,
                f"fullpos vertical_interpolate target=pressure levels={','.join(f'{x:g}' for x in request.levels)}",
            )
        return out
    return _interpolate_data_array(
        values,
        request=request,
        chunks=chunks,
        variable_name=str(values.name or ""),
        keep_attrs=keep_attrs,
    )


def pressure_capabilities() -> dict[str, str]:
    """Return the current pressure-level implementation status."""
    return {
        "status": "native",
        "native_backend": "FULLPOS",
        "input_levels": "hybrid model levels",
        "output_levels": "pressure levels",
    }


def _interpolate_data_array(
    obj: xr.DataArray,
    *,
    request: PressureRequest,
    chunks: dict[str, int] | None,
    variable_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    kernel = _vertical_native.pressure_ppt if _is_temperature_name(variable_name, obj) else _vertical_native.pressure_ppq
    return _apply_native_scalar_kernel(
        obj,
        request=request,
        chunks=chunks,
        keep_attrs=keep_attrs,
        kernel=kernel,
    )


def _interpolate_wind_pair(
    u: xr.DataArray,
    v: xr.DataArray,
    *,
    request: PressureRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> tuple[xr.DataArray, xr.DataArray]:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    if u.dims != v.dims or u.shape != v.shape:
        raise ValueError("u and v variables must have matching dimensions and shapes for FULLPOS PPUV")
    _validate_same_coords(u, v)
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("pressure" if dim == hybrid_dim else dim for dim in u.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else u.sizes[dim] for dim in u.dims)
    u_values = np.empty(out_shape, dtype=np.float64)
    v_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(u, hybrid_dim=hybrid_dim, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        ps_block = _surface_block(request.surface_pressure, selection)
        u_native, v_native = _vertical_native.pressure_ppuv(
            _flatten_hybrid_columns(u_block, hybrid_dim),
            _flatten_hybrid_columns(v_block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            request.levels,
        )
        target = _output_selection(selection, u, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        u_values[target] = _unflatten_output_columns(u_native, u_block, hybrid_dim, request.levels.size)
        v_values[target] = _unflatten_output_columns(v_native, v_block, hybrid_dim, request.levels.size)
    return (
        _wrap_pressure_output(u, u_values, out_dims, request.levels, keep_attrs=keep_attrs),
        _wrap_pressure_output(v, v_values, out_dims, request.levels, keep_attrs=keep_attrs),
    )


def _apply_native_scalar_kernel(
    obj: xr.DataArray,
    *,
    request: PressureRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    kernel,
) -> xr.DataArray:
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("pressure" if dim == hybrid_dim else dim for dim in obj.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else obj.sizes[dim] for dim in obj.dims)
    out_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(obj, hybrid_dim=hybrid_dim, chunks=chunks):
        block = obj.isel(selection) if selection else obj
        ps_block = _surface_block(request.surface_pressure, selection)
        native = kernel(
            _flatten_hybrid_columns(block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            request.levels,
        )
        target = _output_selection(selection, obj, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        out_values[target] = _unflatten_output_columns(native, block, hybrid_dim, request.levels.size)
    return _wrap_pressure_output(obj, out_values, out_dims, request.levels, keep_attrs=keep_attrs)


def _leading_chunk_selections(
    obj: xr.DataArray,
    *,
    hybrid_dim: str,
    chunks: dict[str, int] | None,
):
    leading_dims = [dim for dim in obj.dims if dim != hybrid_dim]
    chunked = []
    for dim in leading_dims:
        size = obj.sizes[dim]
        step = size if chunks is None or dim not in chunks else int(chunks[dim])
        starts = range(0, size, step)
        chunked.append((dim, [slice(start, min(start + step, size)) for start in starts]))
    if not chunked:
        yield {}
        return
    selections = [{}]
    for dim, slices in chunked:
        selections = [dict(sel, **{dim: slc}) for sel in selections for slc in slices]
    for selection in selections:
        yield selection


def _surface_block(ps: xr.DataArray, selection: dict[str, slice]) -> xr.DataArray:
    ps_selection = {dim: value for dim, value in selection.items() if dim in ps.dims}
    return ps.isel(ps_selection) if ps_selection else ps


def _flatten_hybrid_columns(obj: xr.DataArray, hybrid_dim: str) -> np.ndarray:
    other_dims = [dim for dim in obj.dims if dim != hybrid_dim]
    transposed = obj.transpose(*other_dims, hybrid_dim)
    arr = np.asarray(transposed.values, dtype=np.float64)
    return np.ascontiguousarray(arr.reshape((-1, obj.sizes[hybrid_dim])))


def _flatten_surface_columns(ps: xr.DataArray) -> np.ndarray:
    arr = np.asarray(ps.values, dtype=np.float64)
    return np.ascontiguousarray(arr.reshape(-1))


def _unflatten_output_columns(
    native: np.ndarray,
    source_block: xr.DataArray,
    hybrid_dim: str,
    nlevels: int,
) -> np.ndarray:
    other_dims = [dim for dim in source_block.dims if dim != hybrid_dim]
    shape_other = tuple(source_block.sizes[dim] for dim in other_dims)
    arr = np.asarray(native, dtype=np.float64).reshape(shape_other + (nlevels,))
    dims_other_pressure = tuple(other_dims) + ("pressure",)
    target_dims = tuple("pressure" if dim == hybrid_dim else dim for dim in source_block.dims)
    return xr.DataArray(arr, dims=dims_other_pressure).transpose(*target_dims).values


def _output_selection(
    selection: dict[str, slice],
    obj: xr.DataArray,
    *,
    hybrid_dim: str,
    nlevels: int,
) -> tuple:
    return tuple(slice(0, nlevels) if dim == hybrid_dim else selection.get(dim, slice(None)) for dim in obj.dims)


def _wrap_pressure_output(
    template: xr.DataArray,
    values: np.ndarray,
    dims: tuple[str, ...],
    levels: np.ndarray,
    *,
    keep_attrs: bool,
) -> xr.DataArray:
    coords = {}
    for old_dim, new_dim in zip(template.dims, dims):
        if new_dim == "pressure":
            coords[new_dim] = levels
        elif old_dim in template.coords and template.coords[old_dim].ndim == 1:
            coords[new_dim] = template.coords[old_dim].values
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["vertical_target"] = "pressure"
    attrs["vertical_backend"] = "FULLPOS"
    attrs["pressure_units"] = "Pa"
    if keep_attrs:
        attrs = append_history(
            attrs,
            f"fullpos vertical_interpolate target=pressure levels={','.join(f'{x:g}' for x in levels)}",
        )
    return xr.DataArray(values, dims=dims, coords=coords, attrs=attrs, name=template.name)


def _validate_same_coords(a: xr.DataArray, b: xr.DataArray) -> None:
    for dim in a.dims:
        if dim in a.coords and dim in b.coords and a.coords[dim].ndim == 1:
            if not np.array_equal(a.coords[dim].values, b.coords[dim].values):
                raise ValueError(f"u and v coordinate {dim!r} must match")


def _is_temperature_name(name: str, obj: xr.DataArray) -> bool:
    short = str(name or obj.attrs.get("GRIB_shortName", "") or obj.name or "").lower()
    return short in {"t", "temperature"}


def prepare_pressure_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
) -> PressureRequest:
    """Validate and normalize a hybrid-to-pressure request."""
    normalized_levels = _normalize_pressure_levels(levels)
    reference, selected = _validate_pressure_request(
        values,
        levels=normalized_levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
    )
    hybrid_dim = _find_hybrid_dim(reference)
    assert hybrid_dim is not None
    ak, bk = _resolve_hybrid_coefficients(
        reference,
        hybrid_coefficients=hybrid_coefficients,
    )
    expected_half_levels = reference.sizes[hybrid_dim] + 1
    ak, bk = _select_hybrid_coefficients_for_levels(
        reference,
        hybrid_dim=hybrid_dim,
        ak=ak,
        bk=bk,
        expected_half_levels=expected_half_levels,
    )
    if ak.size != expected_half_levels or bk.size != expected_half_levels:
        raise ValueError(
            "hybrid coefficients do not match the input hybrid dimension: "
            f"expected {expected_half_levels} half levels, got {ak.size} A and {bk.size} B values"
        )
    ps = _resolve_surface_pressure(reference, surface_pressure=surface_pressure)
    return PressureRequest(
        levels=normalized_levels,
        hybrid_dim=hybrid_dim,
        ak=ak,
        bk=bk,
        surface_pressure=ps,
        variables=selected,
    )


def _normalize_pressure_levels(levels) -> np.ndarray:
    arr = np.asarray(levels, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError("pressure interpolation requires at least one target level")
    if not np.isfinite(arr).all():
        raise ValueError("pressure levels must be finite")
    if np.any(arr <= 0.0):
        raise ValueError("pressure levels must be positive")
    return arr


def _validate_pressure_request(
    values,
    *,
    levels: np.ndarray,
    variables,
    chunks: dict[str, int] | None,
    surface_pressure,
) -> tuple[xr.DataArray, tuple[str, ...] | None]:
    _ = levels
    if isinstance(values, xr.Dataset):
        reference, selected = _validate_pressure_dataset(values, variables=variables, chunks=chunks)
    elif isinstance(values, xr.DataArray):
        _validate_pressure_data_array(values, chunks=chunks)
        reference = values
        selected = None
    else:
        raise TypeError(
            "pressure interpolation currently requires xarray DataArray or Dataset input"
        )

    if surface_pressure is not None and isinstance(values, (xr.Dataset, xr.DataArray)):
        if not isinstance(surface_pressure, xr.DataArray):
            raise TypeError("surface_pressure must be an xarray.DataArray when values is xarray")
    return reference, selected


def _validate_pressure_dataset(
    obj: xr.Dataset,
    *,
    variables,
    chunks: dict[str, int] | None,
) -> tuple[xr.DataArray, tuple[str, ...]]:
    selected = list(obj.data_vars) if variables is None else [str(v) for v in variables]
    missing = [name for name in selected if name not in obj.data_vars]
    if missing:
        raise KeyError(f"variables not found in dataset: {missing}")
    for name in selected:
        _validate_pressure_data_array(obj[name], chunks=chunks, variable_name=name)
    return obj[selected[0]], tuple(selected)


def _validate_pressure_data_array(
    obj: xr.DataArray,
    *,
    chunks: dict[str, int] | None,
    variable_name: str | None = None,
) -> None:
    hybrid_dim = _find_hybrid_dim(obj)
    if chunks is not None:
        unknown = set(chunks) - set(obj.dims)
        if unknown:
            raise ValueError(f"chunks contains dimensions not present in object: {sorted(unknown)}")
        for dim, size in chunks.items():
            if int(size) <= 0:
                raise ValueError("chunk sizes must be positive")
    if hybrid_dim is None:
        label = "data array" if variable_name is None else f"variable {variable_name!r}"
        raise ValueError(f"{label} does not have a supported hybrid/model-level dimension")


def _find_hybrid_dim(obj: xr.DataArray) -> str | None:
    for dim in ("hybrid", "level", "model_level", "ml"):
        if dim in obj.dims:
            return dim
    return None


def _resolve_hybrid_coefficients(
    obj: xr.DataArray,
    *,
    hybrid_coefficients,
) -> tuple[np.ndarray, np.ndarray]:
    if hybrid_coefficients is not None:
        return _coefficients_from_source(hybrid_coefficients)
    pv = obj.attrs.get("GRIB_pv")
    if pv is not None:
        return _coefficients_from_grib_pv(pv)
    raise ValueError(
        "hybrid coefficients are required; pass hybrid_coefficients=... or use a GRIB-backed field with GRIB_pv metadata"
    )


def _coefficients_from_grib_pv(pv) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(pv, dtype=np.float64).reshape(-1)
    if arr.size == 0 or arr.size % 2 != 0:
        raise ValueError("GRIB_pv must contain an even number of A/B half-level coefficients")
    half = arr.size // 2
    return arr[:half].copy(), arr[half:].copy()


def _select_hybrid_coefficients_for_levels(
    obj: xr.DataArray,
    *,
    hybrid_dim: str,
    ak: np.ndarray,
    bk: np.ndarray,
    expected_half_levels: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Select matching half-level coefficients for a contiguous level subset."""
    if ak.size == expected_half_levels and bk.size == expected_half_levels:
        return ak, bk
    if ak.size != bk.size or hybrid_dim not in obj.coords:
        return ak, bk
    coord = np.asarray(obj.coords[hybrid_dim].values)
    if coord.size + 1 != expected_half_levels:
        return ak, bk
    try:
        levels = coord.astype(np.int64)
    except (TypeError, ValueError):
        return ak, bk
    if not np.allclose(coord, levels):
        return ak, bk
    if levels.size == 0 or not np.array_equal(np.diff(levels), np.ones(levels.size - 1, dtype=np.int64)):
        return ak, bk

    candidates = []
    if levels.min() >= 1:
        candidates.append((int(levels.min()) - 1, int(levels.max()) + 1))
    if levels.min() >= 0:
        candidates.append((int(levels.min()), int(levels.max()) + 2))
    for start, stop in candidates:
        if start >= 0 and stop <= ak.size and stop - start == expected_half_levels:
            return ak[start:stop].copy(), bk[start:stop].copy()
    return ak, bk


def _midlevel_coefficients_from_half_levels(
    ak: np.ndarray,
    bk: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Collapse hybrid half-level coefficients to model-level midpoints."""
    ak = np.asarray(ak, dtype=np.float64).reshape(-1)
    bk = np.asarray(bk, dtype=np.float64).reshape(-1)
    if ak.shape != bk.shape:
        raise ValueError(
            "hybrid half-level coefficients must have matching shapes: "
            f"{ak.shape} vs {bk.shape}"
        )
    if ak.size < 2:
        raise ValueError("at least two half levels are required to form one model level")
    return 0.5 * (ak[:-1] + ak[1:]), 0.5 * (bk[:-1] + bk[1:])


def _infer_reference_pressure_from_ak(ak: np.ndarray, *, units: str | None = None) -> float:
    """Infer whether hybrid A coefficients use normalized or pressure units."""
    if units is not None:
        normalized = units.strip().lower()
        if normalized in {"pa", "pascal", "pascals"}:
            return 1.0
    ak = np.asarray(ak, dtype=np.float64).reshape(-1)
    if ak.size and float(np.nanmax(np.abs(ak))) > 10.0:
        return 1.0
    return 100000.0


def _coefficients_from_source(source) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(source, xr.Dataset):
        ds = source
    else:
        path = Path(source)
        ds = xr.open_dataset(path)
    if "hyai" in ds and "hybi" in ds:
        return (
            np.asarray(ds["hyai"].values, dtype=np.float64).reshape(-1),
            np.asarray(ds["hybi"].values, dtype=np.float64).reshape(-1),
        )
    if "a" in ds and "b" in ds:
        return (
            np.asarray(ds["a"].values, dtype=np.float64).reshape(-1),
            np.asarray(ds["b"].values, dtype=np.float64).reshape(-1),
        )
    raise ValueError("hybrid_coefficients must provide hyai/hybi or a/b variables")


def _resolve_surface_pressure(
    reference: xr.DataArray,
    *,
    surface_pressure,
) -> xr.DataArray:
    ps = surface_pressure
    if ps is None:
        raise ValueError(
            "surface_pressure is required for pressure interpolation; pass sp/lnsp as an xarray.DataArray"
        )
    ps = _normalize_surface_pressure(ps)
    return _align_surface_pressure(reference, ps)


def _normalize_surface_pressure(ps: xr.DataArray) -> xr.DataArray:
    name = str(ps.name or "").lower()
    short_name = str(ps.attrs.get("GRIB_shortName", "")).lower()
    if name == "lnsp" or short_name == "lnsp":
        out = xr.apply_ufunc(np.exp, ps)
        out.name = "sp"
        out.attrs = dict(ps.attrs)
        out.attrs["GRIB_shortName"] = "sp"
        out.attrs["long_name"] = "surface pressure"
        out.attrs["units"] = "Pa"
        return out
    return ps


def _align_surface_pressure(
    reference: xr.DataArray,
    surface_pressure: xr.DataArray,
) -> xr.DataArray:
    ref = reference
    ps = surface_pressure
    if "time" in ref.dims:
        if "time" not in ps.dims:
            raise ValueError("surface_pressure must include a time dimension when the input field has time")
        try:
            ps = ps.sel(time=ref["time"])
        except Exception as exc:
            raise ValueError("surface_pressure time coordinate does not cover the input field times") from exc

    ref_horizontal_dims = tuple(dim for dim in ref.dims if dim != _find_hybrid_dim(ref) and dim != "time")
    ps_dims = tuple(dim for dim in ps.dims if dim != "time")
    if ps_dims != ref_horizontal_dims:
        raise ValueError(
            "surface_pressure horizontal dimensions must match the input field after removing the hybrid dimension: "
            f"expected {ref_horizontal_dims}, got {ps_dims}"
        )
    for dim in ps_dims:
        if ps.sizes[dim] != ref.sizes[dim]:
            raise ValueError(
                f"surface_pressure dimension {dim!r} has size {ps.sizes[dim]}, expected {ref.sizes[dim]}"
            )
    return ps.transpose(*([dim for dim in ref.dims if dim != _find_hybrid_dim(ref)]))
