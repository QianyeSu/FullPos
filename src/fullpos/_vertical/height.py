from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

from ..metadata import append_history
from ..native import add_native_runtime_dir
from .potential_temperature import _resolve_temperature
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
    _validate_pressure_data_array,
    _validate_same_coords,
)


@dataclass(frozen=True)
class HeightTargetRequest:
    """Validated height-target request for native FULLPOS ``FPPS``."""

    levels: np.ndarray
    target_levels_m: np.ndarray
    target_name: str
    hybrid_dim: str
    ak: np.ndarray
    bk: np.ndarray
    surface_pressure: xr.DataArray
    temperature: xr.DataArray
    orography_geopotential: xr.DataArray
    specific_humidity: xr.DataArray | None
    variables: tuple[str, ...] | None


def interpolate_to_height_above_orography(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
    keep_attrs: bool = True,
):
    """Interpolate fields to heights above local orography with FULLPOS.

    Target pressures are computed by native OpenIFS/FULLPOS ``GPHPRE``,
    ``GPGEO``, and ``FPPS``. The resulting per-column pressures are then
    passed to the native ``PPQ``/``PPUV``/``PPT`` interpolation kernels.
    ``levels`` are heights in metres above the supplied surface geopotential.
    """
    return _interpolate_to_height(
        values,
        target_name="height_above_orography",
        levels=levels,
        target_levels_m=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
        keep_attrs=keep_attrs,
    )


def interpolate_to_height_above_sea(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
    keep_attrs: bool = True,
):
    """Interpolate fields to absolute heights above mean sea level with FULLPOS.

    This follows the FULLPOS ``CDCONF='F'`` pressure lookup: target
    geopotential is ``g * level``. The input surface geopotential is still
    required so ``GPGEO`` and ``FPPS`` can reconstruct the source column
    geometry.
    """
    return _interpolate_to_height(
        values,
        target_name="height_above_sea",
        levels=levels,
        target_levels_m=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
        keep_attrs=keep_attrs,
    )


def interpolate_to_flight_level(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
    keep_attrs: bool = True,
):
    """Interpolate fields to FULLPOS flight levels.

    ``levels`` are flight-level numbers, so ``FL350`` is passed as ``350``.
    The native pressure lookup follows the FULLPOS flight-level branch through
    ``FPPS`` after converting the level numbers to metres above mean sea level.
    """
    flight_levels = _normalize_flight_levels(levels)
    return _interpolate_to_height(
        values,
        target_name="flight_level",
        levels=flight_levels,
        target_levels_m=flight_levels * 30.48,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
        keep_attrs=keep_attrs,
    )


def _interpolate_to_height(
    values,
    *,
    target_name: str,
    levels,
    target_levels_m,
    variables,
    chunks: dict[str, int] | None,
    surface_pressure,
    hybrid_coefficients,
    temperature,
    orography_geopotential,
    specific_humidity,
    keep_attrs: bool,
):
    """Prepare a height request and dispatch to the native FULLPOS kernels."""
    request = _prepare_height_request(
        values,
        target_name=target_name,
        levels=levels,
        target_levels_m=target_levels_m,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
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
                f"fullpos vertical_interpolate target={request.target_name} levels={','.join(f'{x:g}' for x in request.levels)}",
            )
        return out
    return _interpolate_data_array(
        values,
        request=request,
        chunks=chunks,
        variable_name=str(values.name or ""),
        keep_attrs=keep_attrs,
    )


def prepare_height_above_orography_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
) -> HeightTargetRequest:
    """Validate a height-above-orography request before native FULLPOS runs."""
    return _prepare_height_request(
        values,
        target_name="height_above_orography",
        levels=levels,
        target_levels_m=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
    )


def prepare_height_above_sea_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
) -> HeightTargetRequest:
    """Validate a height-above-sea request before native FULLPOS runs."""
    return _prepare_height_request(
        values,
        target_name="height_above_sea",
        levels=levels,
        target_levels_m=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
    )


def prepare_flight_level_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
) -> HeightTargetRequest:
    """Validate a flight-level request before native FULLPOS runs."""
    flight_levels = _normalize_flight_levels(levels)
    return _prepare_height_request(
        values,
        target_name="flight_level",
        levels=flight_levels,
        target_levels_m=flight_levels * 30.48,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        temperature=temperature,
        orography_geopotential=orography_geopotential,
        specific_humidity=specific_humidity,
    )


def _prepare_height_request(
    values,
    *,
    target_name: str,
    levels,
    target_levels_m,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    temperature=None,
    orography_geopotential=None,
    specific_humidity=None,
) -> HeightTargetRequest:
    """Normalize height-target inputs into the Python request for native FULLPOS."""
    if target_name not in {"height_above_orography", "height_above_sea", "flight_level"}:
        raise ValueError(f"unsupported height target {target_name!r}")
    normalized_levels = _normalize_height_levels(levels)
    normalized_target_levels_m = np.asarray(target_levels_m, dtype=np.float64).reshape(-1)
    if normalized_target_levels_m.shape != normalized_levels.shape:
        raise ValueError("target_levels_m must have the same shape as levels")
    if not np.isfinite(normalized_target_levels_m).all():
        raise ValueError("target_levels_m must be finite")
    if np.any(normalized_target_levels_m < 0.0):
        raise ValueError("target_levels_m must be non-negative")
    reference, selected = _validate_height_request(values, variables=variables, chunks=chunks)
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
    orog = _resolve_orography_geopotential(
        values,
        reference,
        ps,
        orography_geopotential=orography_geopotential,
    )
    q = _resolve_specific_humidity(
        values,
        reference,
        specific_humidity=specific_humidity,
        chunks=chunks,
    )
    return HeightTargetRequest(
        levels=normalized_levels,
        target_levels_m=normalized_target_levels_m,
        target_name=target_name,
        hybrid_dim=hybrid_dim,
        ak=ak,
        bk=bk,
        surface_pressure=ps,
        temperature=temp,
        orography_geopotential=orog,
        specific_humidity=q,
        variables=selected,
    )


def _validate_height_request(
    values,
    *,
    variables,
    chunks: dict[str, int] | None,
) -> tuple[xr.DataArray, tuple[str, ...] | None]:
    """Validate xarray inputs before height-target pressure lookup."""
    if isinstance(values, xr.Dataset):
        if variables is None:
            selected = [name for name, obj in values.data_vars.items() if _find_hybrid_dim(obj) is not None]
        else:
            selected = [str(v) for v in variables]
        if not selected:
            raise ValueError("height-above-orography interpolation requires at least one hybrid-level variable")
        missing = [name for name in selected if name not in values.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        for name in selected:
            _validate_pressure_data_array(values[name], chunks=chunks, variable_name=name)
        return values[selected[0]], tuple(selected)
    if isinstance(values, xr.DataArray):
        _validate_pressure_data_array(values, chunks=chunks)
        return values, None
    raise TypeError("height-above-orography interpolation currently requires xarray DataArray or Dataset input")


def _normalize_height_levels(levels) -> np.ndarray:
    """Normalize above-orography or above-sea target heights in metres."""
    arr = np.asarray(levels, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError("height-above-orography interpolation requires at least one target height")
    if not np.isfinite(arr).all():
        raise ValueError("height-above-orography levels must be finite")
    if np.any(arr < 0.0):
        raise ValueError("height-above-orography levels must be non-negative")
    return arr


def _normalize_flight_levels(levels) -> np.ndarray:
    """Normalize flight-level numbers before conversion to metres."""
    arr = np.asarray(levels, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError("flight-level interpolation requires at least one target level")
    if not np.isfinite(arr).all():
        raise ValueError("flight-level levels must be finite")
    if np.any(arr < 0.0):
        raise ValueError("flight-level levels must be non-negative")
    return arr


def _resolve_orography_geopotential(
    values,
    reference: xr.DataArray,
    surface_pressure: xr.DataArray,
    *,
    orography_geopotential,
) -> xr.DataArray:
    """Resolve the surface geopotential needed by native height targets."""
    if orography_geopotential is None:
        if isinstance(values, xr.Dataset):
            for name in ("orography_geopotential", "surface_geopotential", "z", "fis", "geopotential", "orog"):
                if name in values and _find_hybrid_dim(values[name]) is None:
                    orography_geopotential = values[name]
                    break
        if orography_geopotential is None:
            raise ValueError(
                "orography_geopotential is required for height targets; "
                "pass ECMWF surface geopotential 'z' in m2 s-2"
            )
    if not isinstance(orography_geopotential, xr.DataArray):
        raise TypeError("orography_geopotential must be an xarray.DataArray")
    return _align_surface_like(
        reference,
        surface_pressure,
        orography_geopotential,
        field_name="orography_geopotential",
    )


def _resolve_specific_humidity(
    values,
    reference: xr.DataArray,
    *,
    specific_humidity,
    chunks: dict[str, int] | None,
) -> xr.DataArray | None:
    """Resolve optional specific humidity for moist native height lookup."""
    if specific_humidity is None:
        if isinstance(values, xr.Dataset) and "q" in values and _find_hybrid_dim(values["q"]) is not None:
            specific_humidity = values["q"]
        elif isinstance(values, xr.DataArray) and str(values.name or "").lower() in {"q", "specific_humidity"}:
            specific_humidity = values
        else:
            return None
    if not isinstance(specific_humidity, xr.DataArray):
        raise TypeError("specific_humidity must be an xarray.DataArray")
    _validate_pressure_data_array(specific_humidity, chunks=chunks, variable_name="specific_humidity")
    if specific_humidity.dims != reference.dims or specific_humidity.shape != reference.shape:
        raise ValueError("specific_humidity must have the same dimensions and shape as the interpolated fields")
    _validate_same_coords(reference, specific_humidity)
    return specific_humidity


def _align_surface_like(
    reference: xr.DataArray,
    template: xr.DataArray,
    field: xr.DataArray,
    *,
    field_name: str,
) -> xr.DataArray:
    """Align a surface field to the horizontal grid used by native FULLPOS."""
    target_dims = template.dims
    unknown = set(field.dims) - set(target_dims)
    if unknown:
        raise ValueError(f"{field_name} contains dimensions not present in the surface grid: {sorted(unknown)}")
    aligned = field
    for dim in field.dims:
        if dim in template.coords and dim in aligned.coords and template.coords[dim].ndim == 1 and aligned.coords[dim].ndim == 1:
            try:
                aligned = aligned.sel({dim: template.coords[dim]})
            except Exception as exc:
                raise ValueError(f"{field_name} coordinate {dim!r} does not cover the input field") from exc
        elif aligned.sizes[dim] != template.sizes[dim]:
            raise ValueError(
                f"{field_name} dimension {dim!r} has size {aligned.sizes[dim]}, expected {template.sizes[dim]}"
            )

    for dim in target_dims:
        if dim in aligned.dims:
            continue
        coord = template.coords[dim].values if dim in template.coords and template.coords[dim].ndim == 1 else np.arange(template.sizes[dim])
        aligned = aligned.expand_dims({dim: coord})

    aligned = aligned.transpose(*target_dims)
    for dim in target_dims:
        if aligned.sizes[dim] != template.sizes[dim]:
            raise ValueError(
                f"{field_name} dimension {dim!r} has size {aligned.sizes[dim]}, expected {template.sizes[dim]}"
            )
    return aligned


def _interpolate_data_array(
    obj: xr.DataArray,
    *,
    request: HeightTargetRequest,
    chunks: dict[str, int] | None,
    variable_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    """Dispatch a scalar field to the native height interpolation kernels."""
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
    request: HeightTargetRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> tuple[xr.DataArray, xr.DataArray]:
    """Dispatch a wind pair to the native PPUV kernel for height targets."""
    add_native_runtime_dir()
    from fullpos import _vertical_native

    if u.dims != v.dims or u.shape != v.shape:
        raise ValueError("u and v variables must have matching dimensions and shapes for FULLPOS PPUV")
    _validate_same_coords(u, v)
    hybrid_dim = request.hybrid_dim
    out_dims = tuple(request.target_name if dim == hybrid_dim else dim for dim in u.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else u.sizes[dim] for dim in u.dims)
    u_values = np.empty(out_shape, dtype=np.float64)
    v_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(u, hybrid_dim=hybrid_dim, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        target_pressures = _height_target_pressures(request=request, selection=selection)
        u_native, v_native = _vertical_native.column_pressure_ppuv(
            _flatten_hybrid_columns(u_block, hybrid_dim),
            _flatten_hybrid_columns(v_block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(_surface_block(request.surface_pressure, selection)),
            target_pressures,
        )
        target = _output_selection(selection, u, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        u_values[target] = _unflatten_output_columns(u_native, u_block, hybrid_dim, request.levels.size)
        v_values[target] = _unflatten_output_columns(v_native, v_block, hybrid_dim, request.levels.size)
    return (
        _wrap_height_output(u, u_values, out_dims, request.levels, target_name=request.target_name, keep_attrs=keep_attrs),
        _wrap_height_output(v, v_values, out_dims, request.levels, target_name=request.target_name, keep_attrs=keep_attrs),
    )


def _apply_native_scalar_kernel(
    obj: xr.DataArray,
    *,
    request: HeightTargetRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    kernel,
) -> xr.DataArray:
    """Apply a native scalar height kernel block-by-block."""
    hybrid_dim = request.hybrid_dim
    out_dims = tuple(request.target_name if dim == hybrid_dim else dim for dim in obj.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else obj.sizes[dim] for dim in obj.dims)
    out_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(obj, hybrid_dim=hybrid_dim, chunks=chunks):
        block = obj.isel(selection) if selection else obj
        target_pressures = _height_target_pressures(request=request, selection=selection)
        native = kernel(
            _flatten_hybrid_columns(block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(_surface_block(request.surface_pressure, selection)),
            target_pressures,
        )
        target = _output_selection(selection, obj, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        out_values[target] = _unflatten_output_columns(native, block, hybrid_dim, request.levels.size)
    return _wrap_height_output(obj, out_values, out_dims, request.levels, target_name=request.target_name, keep_attrs=keep_attrs)


def _height_target_pressures(
    *,
    request: HeightTargetRequest,
    selection: dict[str, slice],
) -> np.ndarray:
    """Compute native target pressures for a height or flight-level surface."""
    from fullpos import _vertical_native

    temperature = request.temperature.isel(selection) if selection else request.temperature
    ps = _surface_block(request.surface_pressure, selection)
    orog = _surface_block(request.orography_geopotential, selection)
    q = None
    if request.specific_humidity is not None:
        q = request.specific_humidity.isel(selection) if selection else request.specific_humidity

    pressure_func = (
        _vertical_native.height_above_sea_pressures
        if request.target_name in {"height_above_sea", "flight_level"}
        else _vertical_native.height_above_orography_pressures
    )

    return pressure_func(
        _flatten_hybrid_columns(temperature, request.hybrid_dim),
        None if q is None else _flatten_hybrid_columns(q, request.hybrid_dim),
        request.ak,
        request.bk,
        _flatten_surface_columns(ps),
        _flatten_surface_columns(orog),
        request.target_levels_m,
    )


def _wrap_height_output(
    template: xr.DataArray,
    values: np.ndarray,
    dims: tuple[str, ...],
    levels: np.ndarray,
    *,
    target_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    """Attach Python-side metadata after native height interpolation."""
    coords = {}
    for old_dim, new_dim in zip(template.dims, dims):
        if new_dim == target_name:
            coords[new_dim] = levels
        elif old_dim in template.coords and template.coords[old_dim].ndim == 1:
            coords[new_dim] = template.coords[old_dim].values
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["vertical_target"] = target_name
    attrs["vertical_backend"] = "FULLPOS"
    attrs["height_units"] = "FL" if target_name == "flight_level" else "m"
    if target_name == "flight_level":
        attrs["flight_level_units"] = "hundreds of feet"
    attrs["orography_geopotential_units"] = "m2 s-2"
    attrs["vertical_native_path"] = "GPHPRE/GPGEO/FPPS + PPINIT/PPFLEV/PPQ-PPUV-PPT"
    if keep_attrs:
        attrs = append_history(
            attrs,
            f"fullpos vertical_interpolate target={target_name} levels={','.join(f'{x:g}' for x in levels)}",
        )
    return xr.DataArray(values, dims=dims, coords=coords, attrs=attrs, name=template.name)
