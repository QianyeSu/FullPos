from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

from ..metadata import append_history
from ..native import add_native_runtime_dir
from .pressure import (
    _find_hybrid_dim,
    _flatten_hybrid_columns,
    _flatten_surface_columns,
    _is_temperature_name,
    _leading_chunk_selections,
    _output_selection,
    _resolve_hybrid_coefficients,
    _resolve_surface_pressure,
    _select_hybrid_coefficients_for_levels,
    _surface_block,
    _validate_pressure_data_array,
    _validate_same_coords,
)


@dataclass(frozen=True)
class EtaRequest:
    """Validated hybrid-to-eta request for native FULLPOS ``PPLETA``."""

    levels: np.ndarray
    hybrid_dim: str
    ak: np.ndarray
    bk: np.ndarray
    surface_pressure: xr.DataArray
    variables: tuple[str, ...] | None


def interpolate_to_eta(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    keep_attrs: bool = True,
):
    """Interpolate fields to FULLPOS eta/model-level indexes.

    The target pressures are computed by native OpenIFS/FULLPOS ``PPLETA``
    and ``GPHPRE`` from integer eta level indexes. The resulting per-column
    pressures are then passed to the native ``PPQ``/``PPUV``/``PPT`` kernels.
    Python only validates inputs and shapes around that native boundary.
    """
    request = prepare_eta_request(
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
                f"fullpos vertical_interpolate target=eta levels={','.join(str(int(x)) for x in request.levels)}",
            )
        return out
    return _interpolate_data_array(
        values,
        request=request,
        chunks=chunks,
        variable_name=str(values.name or ""),
        keep_attrs=keep_attrs,
    )


def prepare_eta_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
) -> EtaRequest:
    """Normalize an eta request before the native FULLPOS pressure lookup."""
    normalized_levels = _normalize_eta_levels(levels)
    reference, selected = _validate_eta_request(values, variables=variables, chunks=chunks)
    hybrid_dim = _find_hybrid_dim(reference)
    assert hybrid_dim is not None
    ak, bk = _resolve_hybrid_coefficients(reference, hybrid_coefficients=hybrid_coefficients)
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
    if np.any(normalized_levels > reference.sizes[hybrid_dim]):
        raise ValueError(
            "eta levels are FULLPOS integer model-level indexes and must be in "
            f"1..{reference.sizes[hybrid_dim]}"
        )
    ps = _resolve_surface_pressure(reference, surface_pressure=surface_pressure)
    return EtaRequest(
        levels=normalized_levels,
        hybrid_dim=hybrid_dim,
        ak=ak,
        bk=bk,
        surface_pressure=ps,
        variables=selected,
    )


def _validate_eta_request(
    values,
    *,
    variables,
    chunks: dict[str, int] | None,
) -> tuple[xr.DataArray, tuple[str, ...] | None]:
    """Validate xarray inputs before native eta interpolation."""
    if isinstance(values, xr.Dataset):
        selected = list(values.data_vars) if variables is None else [str(v) for v in variables]
        missing = [name for name in selected if name not in values.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        for name in selected:
            _validate_pressure_data_array(values[name], chunks=chunks, variable_name=name)
        return values[selected[0]], tuple(selected)
    if isinstance(values, xr.DataArray):
        _validate_pressure_data_array(values, chunks=chunks)
        return values, None
    raise TypeError("eta interpolation currently requires xarray DataArray or Dataset input")


def _normalize_eta_levels(levels) -> np.ndarray:
    """Normalize eta/model-level indexes before native dispatch."""
    arr = np.asarray(levels, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError("eta interpolation requires at least one target level")
    if not np.isfinite(arr).all():
        raise ValueError("eta levels must be finite")
    rounded = np.rint(arr)
    if not np.allclose(arr, rounded, rtol=0.0, atol=1.0e-9):
        raise ValueError("eta levels must be integer FULLPOS model-level indexes")
    if np.any(rounded < 1):
        raise ValueError("eta levels must be positive 1-based model-level indexes")
    return rounded.astype(np.int32)


def _interpolate_data_array(
    obj: xr.DataArray,
    *,
    request: EtaRequest,
    chunks: dict[str, int] | None,
    variable_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    """Dispatch a scalar field to the native eta pressure kernel."""
    add_native_runtime_dir()
    from fullpos import _vertical_native

    kernel = _vertical_native.column_pressure_ppt if _is_temperature_name(variable_name, obj) else _vertical_native.column_pressure_ppq
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
    request: EtaRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Dispatch a wind pair to the native eta PPUV kernel."""
    add_native_runtime_dir()
    from fullpos import _vertical_native

    if u.dims != v.dims or u.shape != v.shape:
        raise ValueError("u and v variables must have matching dimensions and shapes for FULLPOS PPUV")
    _validate_same_coords(u, v)
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("eta" if dim == hybrid_dim else dim for dim in u.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else u.sizes[dim] for dim in u.dims)
    u_values = np.empty(out_shape, dtype=np.float64)
    v_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(u, hybrid_dim=hybrid_dim, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        ps_block = _surface_block(request.surface_pressure, selection)
        target_pressures = _eta_target_pressures(request=request, ps_block=ps_block)
        u_native, v_native = _vertical_native.column_pressure_ppuv(
            _flatten_hybrid_columns(u_block, hybrid_dim),
            _flatten_hybrid_columns(v_block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            target_pressures,
        )
        target = _output_selection(selection, u, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        u_values[target] = _unflatten_eta_output_columns(u_native, u_block, hybrid_dim, request.levels.size)
        v_values[target] = _unflatten_eta_output_columns(v_native, v_block, hybrid_dim, request.levels.size)
    return (
        _wrap_eta_output(u, u_values, out_dims, request.levels, keep_attrs=keep_attrs),
        _wrap_eta_output(v, v_values, out_dims, request.levels, keep_attrs=keep_attrs),
    )


def _apply_native_scalar_kernel(
    obj: xr.DataArray,
    *,
    request: EtaRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    kernel,
) -> xr.DataArray:
    """Apply a native scalar eta kernel block-by-block."""
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("eta" if dim == hybrid_dim else dim for dim in obj.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else obj.sizes[dim] for dim in obj.dims)
    out_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(obj, hybrid_dim=hybrid_dim, chunks=chunks):
        block = obj.isel(selection) if selection else obj
        ps_block = _surface_block(request.surface_pressure, selection)
        target_pressures = _eta_target_pressures(request=request, ps_block=ps_block)
        native = kernel(
            _flatten_hybrid_columns(block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            target_pressures,
        )
        target = _output_selection(selection, obj, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        out_values[target] = _unflatten_eta_output_columns(native, block, hybrid_dim, request.levels.size)
    return _wrap_eta_output(obj, out_values, out_dims, request.levels, keep_attrs=keep_attrs)


def _eta_target_pressures(
    *,
    request: EtaRequest,
    ps_block: xr.DataArray,
) -> np.ndarray:
    """Compute native target pressures for eta model-level indexes."""
    from fullpos import _vertical_native

    return _vertical_native.eta_pressures(
        request.ak,
        request.bk,
        _flatten_surface_columns(ps_block),
        request.levels.astype(np.float64),
    )


def _unflatten_eta_output_columns(
    native: np.ndarray,
    source_block: xr.DataArray,
    hybrid_dim: str,
    nlevels: int,
) -> np.ndarray:
    """Restore native eta output columns to the original xarray layout."""
    other_dims = [dim for dim in source_block.dims if dim != hybrid_dim]
    shape_other = tuple(source_block.sizes[dim] for dim in other_dims)
    arr = np.asarray(native, dtype=np.float64).reshape(shape_other + (nlevels,))
    dims_other_eta = tuple(other_dims) + ("eta",)
    target_dims = tuple("eta" if dim == hybrid_dim else dim for dim in source_block.dims)
    return xr.DataArray(arr, dims=dims_other_eta).transpose(*target_dims).values


def _wrap_eta_output(
    template: xr.DataArray,
    values: np.ndarray,
    dims: tuple[str, ...],
    levels: np.ndarray,
    *,
    keep_attrs: bool,
) -> xr.DataArray:
    """Attach Python-side metadata after native eta interpolation."""
    coords = {}
    dim_map = dict(zip(template.dims, dims))
    for old_dim, new_dim in zip(template.dims, dims):
        if new_dim == "eta":
            coords[new_dim] = levels
        elif old_dim in template.coords and template.coords[old_dim].ndim == 1:
            coords[new_dim] = template.coords[old_dim].values
    for name, coord in template.coords.items():
        if name in coords or _coord_uses_hybrid_dim(coord, template.dims, dims):
            continue
        if all(dim in dim_map for dim in coord.dims):
            coords[name] = (
                tuple(dim_map[dim] for dim in coord.dims),
                coord.values,
                dict(coord.attrs),
            )
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["vertical_target"] = "eta"
    attrs["vertical_backend"] = "FULLPOS"
    attrs["eta_level_kind"] = "FULLPOS integer model-level index"
    attrs["vertical_native_path"] = "PPLETA/GPHPRE + PPINIT/PPFLEV/PPQ-PPUV-PPT"
    if keep_attrs:
        attrs = append_history(
            attrs,
            f"fullpos vertical_interpolate target=eta levels={','.join(str(int(x)) for x in levels)}",
        )
    return xr.DataArray(values, dims=dims, coords=coords, attrs=attrs, name=template.name)


def _coord_uses_hybrid_dim(
    coord: xr.DataArray,
    old_dims: tuple[str, ...],
    new_dims: tuple[str, ...],
) -> bool:
    """Return True for coordinates tied to the original hybrid dimension."""
    hybrid_dims = {old for old, new in zip(old_dims, new_dims) if new == "eta"}
    return any(dim in hybrid_dims for dim in coord.dims)
