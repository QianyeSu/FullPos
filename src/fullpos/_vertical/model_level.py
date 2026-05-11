from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

from ..metadata import append_history
from ..native import add_native_runtime_dir
from .pressure import (
    _coefficients_from_source,
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
    _unflatten_output_columns,
    _validate_pressure_data_array,
    _validate_same_coords,
)


@dataclass(frozen=True)
class ModelLevelRequest:
    """Validated hybrid-to-model-level request prepared for native FULLPOS."""

    hybrid_dim: str
    source_ak: np.ndarray
    source_bk: np.ndarray
    target_ak: np.ndarray
    target_bk: np.ndarray
    surface_pressure: xr.DataArray
    variables: tuple[str, ...] | None


def interpolate_to_model_levels(
    values,
    *,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    target_hybrid_coefficients=None,
    keep_attrs: bool = True,
):
    """Interpolate hybrid fields to native FULLPOS model-level targets.

    Target model levels are represented by target hybrid half-level
    coefficients. The native FULLPOS wrapper computes per-column target
    full-level pressures and then runs the PP-chain interpolation.
    """
    request = prepare_model_level_request(
        values,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        target_hybrid_coefficients=target_hybrid_coefficients,
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
            out.attrs = append_history(out.attrs, "fullpos vertical_interpolate target=model_level")
        return out
    return _interpolate_data_array(
        values,
        request=request,
        chunks=chunks,
        variable_name=str(values.name or ""),
        keep_attrs=keep_attrs,
    )


def prepare_model_level_request(
    values,
    *,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    target_hybrid_coefficients=None,
) -> ModelLevelRequest:
    """Validate and normalize a model-level interpolation request."""
    reference, selected = _validate_model_level_request(values, variables=variables, chunks=chunks)
    hybrid_dim = _find_hybrid_dim(reference)
    assert hybrid_dim is not None
    source_ak, source_bk = _resolve_hybrid_coefficients(
        reference,
        hybrid_coefficients=hybrid_coefficients,
    )
    expected_half_levels = reference.sizes[hybrid_dim] + 1
    source_ak, source_bk = _select_hybrid_coefficients_for_levels(
        reference,
        hybrid_dim=hybrid_dim,
        ak=source_ak,
        bk=source_bk,
        expected_half_levels=expected_half_levels,
    )
    if source_ak.size != expected_half_levels or source_bk.size != expected_half_levels:
        raise ValueError(
            "source hybrid coefficients do not match the input hybrid dimension: "
            f"expected {expected_half_levels} half levels, got {source_ak.size} A and {source_bk.size} B values"
        )
    if target_hybrid_coefficients is None:
        target_ak, target_bk = source_ak.copy(), source_bk.copy()
    else:
        target_ak, target_bk = _coefficients_from_source(target_hybrid_coefficients)
    if target_ak.size != target_bk.size:
        raise ValueError("target hybrid A/B coefficient arrays must have matching lengths")
    if target_ak.size < 2:
        raise ValueError("target hybrid coefficients must contain at least two half levels")
    ps = _resolve_surface_pressure(reference, surface_pressure=surface_pressure)
    return ModelLevelRequest(
        hybrid_dim=hybrid_dim,
        source_ak=source_ak,
        source_bk=source_bk,
        target_ak=target_ak,
        target_bk=target_bk,
        surface_pressure=ps,
        variables=selected,
    )


def _validate_model_level_request(
    values,
    *,
    variables,
    chunks: dict[str, int] | None,
) -> tuple[xr.DataArray, tuple[str, ...] | None]:
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
    raise TypeError("model-level interpolation currently requires xarray DataArray or Dataset input")


def _interpolate_data_array(
    obj: xr.DataArray,
    *,
    request: ModelLevelRequest,
    chunks: dict[str, int] | None,
    variable_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    kernel = _vertical_native.hybrid_pressure_ppt if _is_temperature_name(variable_name, obj) else _vertical_native.hybrid_pressure_ppq
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
    request: ModelLevelRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> tuple[xr.DataArray, xr.DataArray]:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    if u.dims != v.dims or u.shape != v.shape:
        raise ValueError("u and v variables must have matching dimensions and shapes for FULLPOS PPUV")
    _validate_same_coords(u, v)
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("model_level" if dim == hybrid_dim else dim for dim in u.dims)
    nout = request.target_ak.size - 1
    out_shape = tuple(nout if dim == hybrid_dim else u.sizes[dim] for dim in u.dims)
    u_values = np.empty(out_shape, dtype=np.float64)
    v_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(u, hybrid_dim=hybrid_dim, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        ps_block = _surface_block(request.surface_pressure, selection)
        u_native, v_native = _vertical_native.hybrid_pressure_ppuv(
            _flatten_hybrid_columns(u_block, hybrid_dim),
            _flatten_hybrid_columns(v_block, hybrid_dim),
            request.source_ak,
            request.source_bk,
            _flatten_surface_columns(ps_block),
            request.target_ak,
            request.target_bk,
        )
        target = _output_selection(selection, u, hybrid_dim=hybrid_dim, nlevels=nout)
        u_values[target] = _unflatten_output_columns(u_native, u_block, hybrid_dim, nout)
        v_values[target] = _unflatten_output_columns(v_native, v_block, hybrid_dim, nout)
    return (
        _wrap_model_level_output(u, u_values, out_dims, nout, keep_attrs=keep_attrs),
        _wrap_model_level_output(v, v_values, out_dims, nout, keep_attrs=keep_attrs),
    )


def _apply_native_scalar_kernel(
    obj: xr.DataArray,
    *,
    request: ModelLevelRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    kernel,
) -> xr.DataArray:
    hybrid_dim = request.hybrid_dim
    nout = request.target_ak.size - 1
    out_dims = tuple("model_level" if dim == hybrid_dim else dim for dim in obj.dims)
    out_shape = tuple(nout if dim == hybrid_dim else obj.sizes[dim] for dim in obj.dims)
    out_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(obj, hybrid_dim=hybrid_dim, chunks=chunks):
        block = obj.isel(selection) if selection else obj
        ps_block = _surface_block(request.surface_pressure, selection)
        native = kernel(
            _flatten_hybrid_columns(block, hybrid_dim),
            request.source_ak,
            request.source_bk,
            _flatten_surface_columns(ps_block),
            request.target_ak,
            request.target_bk,
        )
        target = _output_selection(selection, obj, hybrid_dim=hybrid_dim, nlevels=nout)
        out_values[target] = _unflatten_output_columns(native, block, hybrid_dim, nout)
    return _wrap_model_level_output(obj, out_values, out_dims, nout, keep_attrs=keep_attrs)


def _target_full_level_pressures(
    request: ModelLevelRequest,
    ps: xr.DataArray,
) -> np.ndarray:
    surface = _flatten_surface_columns(ps)
    half = request.target_ak[None, :] + surface[:, None] * request.target_bk[None, :]
    full = 0.5 * (half[:, :-1] + half[:, 1:])
    if not np.isfinite(full).all() or np.any(full <= 0.0):
        raise ValueError("target hybrid coefficients produce non-positive or non-finite full-level pressures")
    return np.asfortranarray(full, dtype=np.float64)


def _wrap_model_level_output(
    template: xr.DataArray,
    values: np.ndarray,
    dims: tuple[str, ...],
    nlevels: int,
    *,
    keep_attrs: bool,
) -> xr.DataArray:
    coords = {}
    for old_dim, new_dim in zip(template.dims, dims):
        if new_dim == "model_level":
            coords[new_dim] = np.arange(1, nlevels + 1, dtype=np.int32)
        elif old_dim in template.coords and template.coords[old_dim].ndim == 1:
            coords[new_dim] = template.coords[old_dim].values
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["vertical_target"] = "model_level"
    attrs["vertical_backend"] = "FULLPOS"
    attrs["vertical_native_path"] = "PPINIT/PPFLEV/PPQ-PPUV-PPT"
    if keep_attrs:
        attrs = append_history(attrs, "fullpos vertical_interpolate target=model_level")
    return xr.DataArray(values, dims=dims, coords=coords, attrs=attrs, name=template.name)
