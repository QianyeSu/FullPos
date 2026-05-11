from __future__ import annotations

import numpy as np
import xarray as xr

from ..metadata import append_history
from ..native import add_native_runtime_dir
from .potential_temperature import (
    PotentialTemperatureRequest,
    _resolve_temperature,
    _validate_theta_request,
)
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
    _unflatten_output_columns,
    _validate_same_coords,
)


TemperatureRequest = PotentialTemperatureRequest


def interpolate_to_temperature(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    keep_attrs: bool = True,
):
    """Interpolate fields to temperature surfaces with native FULLPOS code.

    Python only prepares xarray inputs and target pressures; the actual
    interpolation is performed by the native FULLPOS/Fortran kernels.
    """
    request = prepare_temperature_request(
        values,
        levels=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
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
                f"fullpos vertical_interpolate target=temperature levels={','.join(f'{x:g}' for x in request.levels)}",
            )
        return out
    return _interpolate_data_array(
        values,
        request=request,
        chunks=chunks,
        variable_name=str(values.name or ""),
        keep_attrs=keep_attrs,
    )


def prepare_temperature_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
) -> TemperatureRequest:
    """Validate and normalize a temperature-surface request for native FULLPOS.

    The returned request is Python-side metadata that feeds the native
    FULLPOS pressure lookup and interpolation kernels.
    """
    normalized_levels = _normalize_temperature_levels(levels)
    reference, selected = _validate_theta_request(values, variables=variables, chunks=chunks)
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
    ps = _resolve_surface_pressure(reference, surface_pressure=surface_pressure)
    temp = _resolve_temperature(values, reference, temperature=temperature, chunks=chunks)
    return TemperatureRequest(
        levels=normalized_levels,
        hybrid_dim=hybrid_dim,
        ak=ak,
        bk=bk,
        surface_pressure=ps,
        temperature=temp,
        variables=selected,
    )


def _normalize_temperature_levels(levels) -> np.ndarray:
    """Normalize temperature targets before they are handed to native FULLPOS."""
    arr = np.asarray(levels, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError("temperature interpolation requires at least one target temperature level")
    if not np.isfinite(arr).all():
        raise ValueError("temperature levels must be finite")
    if np.any(arr <= 0.0):
        raise ValueError("temperature levels must be positive")
    return arr


def _interpolate_data_array(
    obj: xr.DataArray,
    *,
    request: TemperatureRequest,
    chunks: dict[str, int] | None,
    variable_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    """Run the native scalar pressure kernel for one temperature target field."""
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
    request: TemperatureRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Run the native PPUV kernel for one wind pair on temperature surfaces."""
    add_native_runtime_dir()
    from fullpos import _vertical_native

    if u.dims != v.dims or u.shape != v.shape:
        raise ValueError("u and v variables must have matching dimensions and shapes for FULLPOS PPUV")
    _validate_same_coords(u, v)
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("temperature" if dim == hybrid_dim else dim for dim in u.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else u.sizes[dim] for dim in u.dims)
    u_values = np.empty(out_shape, dtype=np.float64)
    v_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(u, hybrid_dim=hybrid_dim, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        temp_block = request.temperature.isel(selection) if selection else request.temperature
        ps_block = _surface_block(request.surface_pressure, selection)
        target_pressures = _temperature_target_pressures(temp_block, request=request, ps_block=ps_block)
        u_native, v_native = _vertical_native.column_pressure_ppuv(
            _flatten_hybrid_columns(u_block, hybrid_dim),
            _flatten_hybrid_columns(v_block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            target_pressures,
        )
        target = _output_selection(selection, u, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        u_values[target] = _unflatten_output_columns(u_native, u_block, hybrid_dim, request.levels.size)
        v_values[target] = _unflatten_output_columns(v_native, v_block, hybrid_dim, request.levels.size)
    return (
        _wrap_temperature_output(u, u_values, out_dims, request.levels, keep_attrs=keep_attrs),
        _wrap_temperature_output(v, v_values, out_dims, request.levels, keep_attrs=keep_attrs),
    )


def _apply_native_scalar_kernel(
    obj: xr.DataArray,
    *,
    request: TemperatureRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    kernel,
) -> xr.DataArray:
    """Apply a native scalar pressure kernel block-by-block along leading dims."""
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("temperature" if dim == hybrid_dim else dim for dim in obj.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else obj.sizes[dim] for dim in obj.dims)
    out_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(obj, hybrid_dim=hybrid_dim, chunks=chunks):
        block = obj.isel(selection) if selection else obj
        temp_block = request.temperature.isel(selection) if selection else request.temperature
        ps_block = _surface_block(request.surface_pressure, selection)
        target_pressures = _temperature_target_pressures(temp_block, request=request, ps_block=ps_block)
        native = kernel(
            _flatten_hybrid_columns(block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            target_pressures,
        )
        target = _output_selection(selection, obj, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        out_values[target] = _unflatten_output_columns(native, block, hybrid_dim, request.levels.size)
    return _wrap_temperature_output(obj, out_values, out_dims, request.levels, keep_attrs=keep_attrs)


def _temperature_target_pressures(
    temperature: xr.DataArray,
    *,
    request: TemperatureRequest,
    ps_block: xr.DataArray,
) -> np.ndarray:
    """Compute native target pressures for temperature surfaces."""
    from fullpos import _vertical_native

    return _vertical_native.temperature_pressures(
        _flatten_hybrid_columns(temperature, request.hybrid_dim),
        request.ak,
        request.bk,
        _flatten_surface_columns(ps_block),
        request.levels,
    )


def _wrap_temperature_output(
    template: xr.DataArray,
    values: np.ndarray,
    dims: tuple[str, ...],
    levels: np.ndarray,
    *,
    keep_attrs: bool,
) -> xr.DataArray:
    """Attach Python-side metadata after native FULLPOS temperature interpolation."""
    coords = {}
    for old_dim, new_dim in zip(template.dims, dims):
        if new_dim == "temperature":
            coords[new_dim] = levels
        elif old_dim in template.coords and template.coords[old_dim].ndim == 1:
            coords[new_dim] = template.coords[old_dim].values
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["vertical_target"] = "temperature"
    attrs["vertical_backend"] = "FULLPOS"
    attrs["temperature_units"] = "K"
    attrs["vertical_native_path"] = "PPLTW/FPPS + PPINIT/PPFLEV/PPQ-PPUV-PPT"
    if keep_attrs:
        attrs = append_history(
            attrs,
            f"fullpos vertical_interpolate target=temperature levels={','.join(f'{x:g}' for x in levels)}",
        )
    return xr.DataArray(values, dims=dims, coords=coords, attrs=attrs, name=template.name)
