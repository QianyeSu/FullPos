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

_EARTH_OMEGA = 7.2921150e-5


@dataclass(frozen=True)
class PotentialVorticityRequest:
    """Validated PV-surface interpolation request for native FULLPOS."""

    levels: np.ndarray
    hybrid_dim: str
    ak: np.ndarray
    bk: np.ndarray
    surface_pressure: xr.DataArray
    potential_vorticity: xr.DataArray
    coriolis: xr.DataArray
    variables: tuple[str, ...] | None


def interpolate_to_potential_vorticity(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    potential_vorticity=None,
    coriolis=None,
    keep_attrs: bool = True,
):
    """Interpolate fields to iso-PV surfaces with native FULLPOS.

    This target currently uses the native FULLPOS ``PPLTP`` locator on a
    provided full-level potential-vorticity field, then interpolates the
    requested variables with the same native ``PPQ``/``PPUV``/``PPT`` kernels
    used by the other vertical targets.
    """
    request = prepare_potential_vorticity_request(
        values,
        levels=levels,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
        potential_vorticity=potential_vorticity,
        coriolis=coriolis,
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
                f"fullpos vertical_interpolate target=potential_vorticity levels={','.join(f'{x:g}' for x in request.levels)}",
            )
        return out
    return _interpolate_data_array(
        values,
        request=request,
        chunks=chunks,
        variable_name=str(values.name or ""),
        keep_attrs=keep_attrs,
    )


def prepare_potential_vorticity_request(
    values,
    *,
    levels,
    variables=None,
    chunks: dict[str, int] | None = None,
    surface_pressure=None,
    hybrid_coefficients=None,
    potential_vorticity=None,
    coriolis=None,
) -> PotentialVorticityRequest:
    """Validate and normalize an iso-PV interpolation request."""
    normalized_levels = _normalize_pv_levels(levels)
    reference, selected, pv_reference = _validate_pv_request(
        values,
        variables=variables,
        chunks=chunks,
        potential_vorticity=potential_vorticity,
    )
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
    pv = pv_reference
    corio = _resolve_coriolis(reference, hybrid_dim=hybrid_dim, coriolis=coriolis)
    return PotentialVorticityRequest(
        levels=normalized_levels,
        hybrid_dim=hybrid_dim,
        ak=ak,
        bk=bk,
        surface_pressure=ps,
        potential_vorticity=pv,
        coriolis=corio,
        variables=selected,
    )


def _validate_pv_request(
    values,
    *,
    variables,
    chunks: dict[str, int] | None,
    potential_vorticity,
) -> tuple[xr.DataArray, tuple[str, ...] | None, xr.DataArray]:
    if isinstance(values, xr.Dataset):
        selected = list(values.data_vars) if variables is None else [str(v) for v in variables]
        missing = [name for name in selected if name not in values.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        for name in selected:
            _validate_pressure_data_array(values[name], chunks=chunks, variable_name=name)
        reference = values[selected[0]]
    elif isinstance(values, xr.DataArray):
        _validate_pressure_data_array(values, chunks=chunks)
        selected = None
        reference = values
    else:
        raise TypeError("potential-vorticity interpolation currently requires xarray DataArray or Dataset input")
    pv = _resolve_pv_field(values, reference, potential_vorticity=potential_vorticity, chunks=chunks)
    return reference, selected, pv


def _resolve_pv_field(
    values,
    reference: xr.DataArray,
    *,
    potential_vorticity,
    chunks: dict[str, int] | None,
) -> xr.DataArray:
    if potential_vorticity is None:
        if isinstance(values, xr.Dataset):
            for name in ("potential_vorticity", "pv"):
                if name in values:
                    potential_vorticity = values[name]
                    break
        elif isinstance(values, xr.DataArray) and str(values.name or "").lower() in {"potential_vorticity", "pv"}:
            potential_vorticity = values
    if potential_vorticity is None:
        raise ValueError(
            "potential_vorticity is required to locate iso-PV surfaces; pass potential_vorticity=... "
            "or provide a Dataset containing a 'pv' or 'potential_vorticity' variable"
        )
    if not isinstance(potential_vorticity, xr.DataArray):
        raise TypeError("potential_vorticity must be an xarray.DataArray")
    _validate_pressure_data_array(potential_vorticity, chunks=chunks, variable_name="potential_vorticity")
    if potential_vorticity.dims != reference.dims or potential_vorticity.shape != reference.shape:
        raise ValueError("potential_vorticity must have the same dimensions and shape as the interpolated fields")
    _validate_same_coords(reference, potential_vorticity)
    return potential_vorticity


def _resolve_coriolis(
    reference: xr.DataArray,
    *,
    hybrid_dim: str,
    coriolis,
) -> xr.DataArray:
    if coriolis is None:
        coriolis = _infer_coriolis_from_coords(reference)
        if coriolis is None:
            raise ValueError(
                "coriolis is required when it cannot be inferred from latitude coordinates; "
                "pass coriolis=... or provide lat/latitude coordinates"
            )
    if not isinstance(coriolis, xr.DataArray):
        raise TypeError("coriolis must be an xarray.DataArray")
    if hybrid_dim in coriolis.dims:
        raise ValueError("coriolis must not include the hybrid dimension")
    base = reference.isel({hybrid_dim: 0}, drop=True)
    try:
        aligned = coriolis.broadcast_like(base)
    except Exception as exc:
        raise ValueError("coriolis dimensions must be broadcastable to the input horizontal grid") from exc
    return aligned


def _infer_coriolis_from_coords(reference: xr.DataArray) -> xr.DataArray | None:
    lat = None
    for name in ("latitude", "lat"):
        if name in reference.coords:
            lat = reference.coords[name]
            break
    if lat is None:
        return None
    values = 2.0 * _EARTH_OMEGA * np.sin(np.deg2rad(np.asarray(lat.values, dtype=np.float64)))
    return xr.DataArray(values, dims=lat.dims, coords=lat.coords, name="coriolis")


def _normalize_pv_levels(levels) -> np.ndarray:
    arr = np.asarray(levels, dtype=np.float64).reshape(-1)
    if arr.size == 0:
        raise ValueError("potential-vorticity interpolation requires at least one target level")
    if not np.isfinite(arr).all():
        raise ValueError("potential-vorticity levels must be finite")
    if np.any(arr <= 0.0):
        raise ValueError("potential-vorticity levels must be positive")
    return arr


def _interpolate_data_array(
    obj: xr.DataArray,
    *,
    request: PotentialVorticityRequest,
    chunks: dict[str, int] | None,
    variable_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    kernel = _vertical_native.column_pressure_ppt if str(variable_name or "").lower() in {"t", "temperature"} else _vertical_native.column_pressure_ppq
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
    request: PotentialVorticityRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> tuple[xr.DataArray, xr.DataArray]:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    if u.dims != v.dims or u.shape != v.shape:
        raise ValueError("u and v variables must have matching dimensions and shapes for FULLPOS PPUV")
    _validate_same_coords(u, v)
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("potential_vorticity" if dim == hybrid_dim else dim for dim in u.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else u.sizes[dim] for dim in u.dims)
    u_values = np.empty(out_shape, dtype=np.float64)
    v_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(u, hybrid_dim=hybrid_dim, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        ps_block = _surface_block(request.surface_pressure, selection)
        pv_block = request.potential_vorticity.isel(selection) if selection else request.potential_vorticity
        corio_block = request.coriolis.isel(selection) if selection else request.coriolis
        target_pressures = _pv_target_pressures(pv_block, request=request, ps_block=ps_block, coriolis_block=corio_block)
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
        _wrap_pv_output(u, u_values, out_dims, request.levels, keep_attrs=keep_attrs),
        _wrap_pv_output(v, v_values, out_dims, request.levels, keep_attrs=keep_attrs),
    )


def _apply_native_scalar_kernel(
    obj: xr.DataArray,
    *,
    request: PotentialVorticityRequest,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    kernel,
) -> xr.DataArray:
    hybrid_dim = request.hybrid_dim
    out_dims = tuple("potential_vorticity" if dim == hybrid_dim else dim for dim in obj.dims)
    out_shape = tuple(request.levels.size if dim == hybrid_dim else obj.sizes[dim] for dim in obj.dims)
    out_values = np.empty(out_shape, dtype=np.float64)
    for selection in _leading_chunk_selections(obj, hybrid_dim=hybrid_dim, chunks=chunks):
        block = obj.isel(selection) if selection else obj
        ps_block = _surface_block(request.surface_pressure, selection)
        pv_block = request.potential_vorticity.isel(selection) if selection else request.potential_vorticity
        corio_block = request.coriolis.isel(selection) if selection else request.coriolis
        target_pressures = _pv_target_pressures(pv_block, request=request, ps_block=ps_block, coriolis_block=corio_block)
        native = kernel(
            _flatten_hybrid_columns(block, hybrid_dim),
            request.ak,
            request.bk,
            _flatten_surface_columns(ps_block),
            target_pressures,
        )
        target = _output_selection(selection, obj, hybrid_dim=hybrid_dim, nlevels=request.levels.size)
        out_values[target] = _unflatten_output_columns(native, block, hybrid_dim, request.levels.size)
    return _wrap_pv_output(obj, out_values, out_dims, request.levels, keep_attrs=keep_attrs)


def _pv_target_pressures(
    pv: xr.DataArray,
    *,
    request: PotentialVorticityRequest,
    ps_block: xr.DataArray,
    coriolis_block: xr.DataArray,
) -> np.ndarray:
    from fullpos import _vertical_native

    return _vertical_native.potential_vorticity_pressures(
        _flatten_hybrid_columns(pv, request.hybrid_dim),
        request.ak,
        request.bk,
        _flatten_surface_columns(ps_block),
        _flatten_surface_columns(coriolis_block),
        request.levels,
    )


def _wrap_pv_output(
    template: xr.DataArray,
    values: np.ndarray,
    dims: tuple[str, ...],
    levels: np.ndarray,
    *,
    keep_attrs: bool,
) -> xr.DataArray:
    coords = {}
    for old_dim, new_dim in zip(template.dims, dims):
        if new_dim == "potential_vorticity":
            coords[new_dim] = levels
        elif old_dim in template.coords and template.coords[old_dim].ndim == 1:
            coords[new_dim] = template.coords[old_dim].values
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["vertical_target"] = "potential_vorticity"
    attrs["vertical_backend"] = "FULLPOS"
    attrs["vertical_native_path"] = "PPLTP + PPINIT/PPFLEV/PPQ-PPUV-PPT"
    if keep_attrs:
        attrs = append_history(
            attrs,
            f"fullpos vertical_interpolate target=potential_vorticity levels={','.join(f'{x:g}' for x in levels)}",
        )
    return xr.DataArray(values, dims=dims, coords=coords, attrs=attrs, name=template.name)
