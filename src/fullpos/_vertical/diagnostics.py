from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import xarray as xr

from ..grids import GaussianGrid, gaussian_latitudes, infer_grid_name_from_attrs, infer_grid_name_from_shape, parse_grid
from ..metadata import append_history
from ..native import add_native_runtime_dir
from .pressure import (
    _find_hybrid_dim,
    _flatten_hybrid_columns,
    _flatten_surface_columns,
    _resolve_hybrid_coefficients,
    _select_hybrid_coefficients_for_levels,
    _leading_chunk_selections,
    _output_selection,
    _unflatten_output_columns,
    _validate_pressure_data_array,
    _validate_same_coords,
)


@dataclass(frozen=True)
class PotentialVorticityDiagnostic:
    """Model-level outputs from native FULLPOS ``GPPVO``."""

    potential_vorticity: xr.DataArray
    potential_temperature: xr.DataArray


def diagnose_potential_vorticity(
    *,
    u: xr.DataArray,
    v: xr.DataArray,
    temperature: xr.DataArray,
    surface_pressure: xr.DataArray,
    coriolis: xr.DataArray | None = None,
    relative_vorticity: xr.DataArray | None = None,
    temperature_meridional_gradient: xr.DataArray | None = None,
    temperature_zonal_gradient: xr.DataArray | None = None,
    surface_pressure_meridional_gradient: xr.DataArray | None = None,
    surface_pressure_zonal_gradient: xr.DataArray | None = None,
    kappa: xr.DataArray | None = None,
    specific_humidity: xr.DataArray | None = None,
    source_grid: str | GaussianGrid | None = None,
    ntrunc: int | None = None,
    chunks: dict[str, int] | None = None,
    hybrid_coefficients=None,
    keep_attrs: bool = True,
) -> PotentialVorticityDiagnostic:
    """Diagnose model-level PV and theta with native FULLPOS ``GPPVO``.

    The numerical chain is native: ECTRANS prepares relative vorticity and
    scalar horizontal derivatives when they are not supplied, FULLPOS
    ``GPRCP`` computes ``kappa`` from specific humidity when needed, and
    FULLPOS ``GPPVO`` computes PV/theta.
    """
    _validate_pv_core_inputs(
        u=u,
        v=v,
        temperature=temperature,
    )
    hybrid_dim = _find_hybrid_dim(temperature)
    assert hybrid_dim is not None
    grid: GaussianGrid | None = None
    horizontal_dims: tuple[str, ...] | None = None
    if _needs_ectrans_diagnostics(
        relative_vorticity=relative_vorticity,
        temperature_meridional_gradient=temperature_meridional_gradient,
        temperature_zonal_gradient=temperature_zonal_gradient,
        surface_pressure_meridional_gradient=surface_pressure_meridional_gradient,
        surface_pressure_zonal_gradient=surface_pressure_zonal_gradient,
    ) or coriolis is None:
        grid = _resolve_grid(temperature, source_grid=source_grid)
        horizontal_dims = _horizontal_dims(temperature, grid)

    ak, bk = _resolve_hybrid_coefficients(temperature, hybrid_coefficients=hybrid_coefficients)
    expected_half_levels = temperature.sizes[hybrid_dim] + 1
    ak, bk = _select_hybrid_coefficients_for_levels(
        temperature,
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
    surface_pressure = _broadcast_surface_vector(surface_pressure, temperature, hybrid_dim=hybrid_dim, name="surface_pressure")
    if coriolis is None:
        assert grid is not None and horizontal_dims is not None
        coriolis = _auto_coriolis(temperature, grid=grid, horizontal_dims=horizontal_dims, hybrid_dim=hybrid_dim)
    else:
        coriolis = _broadcast_surface_vector(coriolis, temperature, hybrid_dim=hybrid_dim, name="coriolis")

    if _needs_ectrans_diagnostics(
        relative_vorticity=relative_vorticity,
        temperature_meridional_gradient=temperature_meridional_gradient,
        temperature_zonal_gradient=temperature_zonal_gradient,
        surface_pressure_meridional_gradient=surface_pressure_meridional_gradient,
        surface_pressure_zonal_gradient=surface_pressure_zonal_gradient,
    ):
        assert grid is not None and horizontal_dims is not None
        prepared = _prepare_ectrans_diagnostics(
            u=u,
            v=v,
            temperature=temperature,
            surface_pressure=surface_pressure,
            grid=grid,
            horizontal_dims=horizontal_dims,
            hybrid_dim=hybrid_dim,
            ntrunc=ntrunc,
            chunks=chunks,
            keep_attrs=keep_attrs,
        )
        if relative_vorticity is None:
            relative_vorticity = prepared.relative_vorticity
        if temperature_meridional_gradient is None:
            temperature_meridional_gradient = prepared.temperature_meridional_gradient
        if temperature_zonal_gradient is None:
            temperature_zonal_gradient = prepared.temperature_zonal_gradient
        if surface_pressure_meridional_gradient is None:
            surface_pressure_meridional_gradient = prepared.surface_pressure_meridional_gradient
        if surface_pressure_zonal_gradient is None:
            surface_pressure_zonal_gradient = prepared.surface_pressure_zonal_gradient

    assert relative_vorticity is not None
    assert temperature_meridional_gradient is not None
    assert temperature_zonal_gradient is not None
    assert surface_pressure_meridional_gradient is not None
    assert surface_pressure_zonal_gradient is not None

    _validate_pv_prepared_inputs(
        u=u,
        relative_vorticity=relative_vorticity,
        temperature_meridional_gradient=temperature_meridional_gradient,
        temperature_zonal_gradient=temperature_zonal_gradient,
    )
    spm = _broadcast_surface_vector(
        surface_pressure_meridional_gradient,
        temperature,
        hybrid_dim=hybrid_dim,
        name="surface_pressure_meridional_gradient",
    )
    spl = _broadcast_surface_vector(
        surface_pressure_zonal_gradient,
        temperature,
        hybrid_dim=hybrid_dim,
        name="surface_pressure_zonal_gradient",
    )
    kappa = _resolve_kappa(
        kappa,
        specific_humidity=specific_humidity,
        reference=temperature,
        hybrid_dim=hybrid_dim,
        chunks=chunks,
        keep_attrs=keep_attrs,
    )

    add_native_runtime_dir()
    from fullpos import _vertical_native

    pv_native, theta_native = _vertical_native.diagnose_potential_vorticity(
        _flatten_hybrid_columns(u, hybrid_dim),
        _flatten_hybrid_columns(v, hybrid_dim),
        _flatten_hybrid_columns(temperature, hybrid_dim),
        _flatten_hybrid_columns(relative_vorticity, hybrid_dim),
        _flatten_hybrid_columns(temperature_meridional_gradient, hybrid_dim),
        _flatten_hybrid_columns(temperature_zonal_gradient, hybrid_dim),
        _flatten_surface_columns(spm),
        _flatten_surface_columns(spl),
        _flatten_hybrid_columns(kappa, hybrid_dim),
        ak,
        bk,
        _flatten_surface_columns(surface_pressure),
        _flatten_surface_columns(coriolis),
    )
    pv = _wrap_diagnostic_output(
        temperature,
        _unflatten_output_columns(pv_native, temperature, hybrid_dim, temperature.sizes[hybrid_dim]),
        name="potential_vorticity",
        attrs_name="potential_vorticity",
        keep_attrs=keep_attrs,
    )
    theta = _wrap_diagnostic_output(
        temperature,
        _unflatten_output_columns(theta_native, temperature, hybrid_dim, temperature.sizes[hybrid_dim]),
        name="potential_temperature",
        attrs_name="potential_temperature",
        keep_attrs=keep_attrs,
    )
    return PotentialVorticityDiagnostic(potential_vorticity=pv, potential_temperature=theta)


@dataclass(frozen=True)
class _PreparedPvFields:
    relative_vorticity: xr.DataArray
    temperature_meridional_gradient: xr.DataArray
    temperature_zonal_gradient: xr.DataArray
    surface_pressure_meridional_gradient: xr.DataArray
    surface_pressure_zonal_gradient: xr.DataArray


def _validate_pv_core_inputs(
    *,
    u: xr.DataArray,
    v: xr.DataArray,
    temperature: xr.DataArray,
) -> None:
    for name, obj in {
        "u": u,
        "v": v,
        "temperature": temperature,
    }.items():
        if not isinstance(obj, xr.DataArray):
            raise TypeError(f"{name} must be an xarray.DataArray")
        _validate_pressure_data_array(obj, chunks=None, variable_name=name)
    for name, obj in {"v": v, "temperature": temperature}.items():
        if obj.dims != u.dims or obj.shape != u.shape:
            raise ValueError(f"{name} must have the same dimensions and shape as u")
        _validate_same_coords(u, obj)


def _validate_pv_prepared_inputs(
    *,
    u: xr.DataArray,
    relative_vorticity: xr.DataArray,
    temperature_meridional_gradient: xr.DataArray,
    temperature_zonal_gradient: xr.DataArray,
) -> None:
    for name, obj in {
        "relative_vorticity": relative_vorticity,
        "temperature_meridional_gradient": temperature_meridional_gradient,
        "temperature_zonal_gradient": temperature_zonal_gradient,
    }.items():
        if not isinstance(obj, xr.DataArray):
            raise TypeError(f"{name} must be an xarray.DataArray")
        _validate_pressure_data_array(obj, chunks=None, variable_name=name)
        if obj.dims != u.dims or obj.shape != u.shape:
            raise ValueError(f"{name} must have the same dimensions and shape as u")
        _validate_same_coords(u, obj)


def _needs_ectrans_diagnostics(
    *,
    relative_vorticity,
    temperature_meridional_gradient,
    temperature_zonal_gradient,
    surface_pressure_meridional_gradient,
    surface_pressure_zonal_gradient,
) -> bool:
    return any(
        obj is None
        for obj in (
            relative_vorticity,
            temperature_meridional_gradient,
            temperature_zonal_gradient,
            surface_pressure_meridional_gradient,
            surface_pressure_zonal_gradient,
        )
    )


def _prepare_ectrans_diagnostics(
    *,
    u: xr.DataArray,
    v: xr.DataArray,
    temperature: xr.DataArray,
    surface_pressure: xr.DataArray,
    grid: GaussianGrid,
    horizontal_dims: tuple[str, ...],
    hybrid_dim: str,
    ntrunc: int | None,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> _PreparedPvFields:
    _validate_spectral_chunks(temperature, hybrid_dim=hybrid_dim, horizontal_dims=horizontal_dims, chunks=chunks)
    ntrunc = _normalize_ntrunc(grid, ntrunc)
    rel_vort_values = np.empty(temperature.shape, dtype=np.float64)
    tm_values = np.empty(temperature.shape, dtype=np.float64)
    tl_values = np.empty(temperature.shape, dtype=np.float64)
    spm_values = np.empty(surface_pressure.shape, dtype=np.float64)
    spl_values = np.empty(surface_pressure.shape, dtype=np.float64)

    for selection in _spectral_chunk_selections(temperature, hybrid_dim=hybrid_dim, horizontal_dims=horizontal_dims, chunks=chunks):
        u_block = u.isel(selection) if selection else u
        v_block = v.isel(selection) if selection else v
        t_block = temperature.isel(selection) if selection else temperature
        sp_block = _surface_selection(surface_pressure, selection)
        block = _run_ectrans_diagnostics_block(
            u_block,
            v_block,
            t_block,
            sp_block,
            grid=grid,
            horizontal_dims=horizontal_dims,
            hybrid_dim=hybrid_dim,
            ntrunc=ntrunc,
        )
        rel_vort_values[_array_selection(selection, temperature)] = block.relative_vorticity.values
        tm_values[_array_selection(selection, temperature)] = block.temperature_meridional_gradient.values
        tl_values[_array_selection(selection, temperature)] = block.temperature_zonal_gradient.values
        spm_values[_array_selection(selection, surface_pressure)] = block.surface_pressure_meridional_gradient.values
        spl_values[_array_selection(selection, surface_pressure)] = block.surface_pressure_zonal_gradient.values

    return _PreparedPvFields(
        relative_vorticity=_wrap_prepared_like(
            temperature,
            rel_vort_values,
            name="relative_vorticity",
            diagnostic="relative_vorticity",
            keep_attrs=keep_attrs,
        ),
        temperature_meridional_gradient=_wrap_prepared_like(
            temperature,
            tm_values,
            name="temperature_meridional_gradient",
            diagnostic="temperature_meridional_gradient",
            keep_attrs=keep_attrs,
        ),
        temperature_zonal_gradient=_wrap_prepared_like(
            temperature,
            tl_values,
            name="temperature_zonal_gradient",
            diagnostic="temperature_zonal_gradient",
            keep_attrs=keep_attrs,
        ),
        surface_pressure_meridional_gradient=_wrap_prepared_like(
            surface_pressure,
            spm_values,
            name="surface_pressure_meridional_gradient",
            diagnostic="surface_pressure_meridional_gradient",
            keep_attrs=keep_attrs,
        ),
        surface_pressure_zonal_gradient=_wrap_prepared_like(
            surface_pressure,
            spl_values,
            name="surface_pressure_zonal_gradient",
            diagnostic="surface_pressure_zonal_gradient",
            keep_attrs=keep_attrs,
        ),
    )


def _run_ectrans_diagnostics_block(
    u: xr.DataArray,
    v: xr.DataArray,
    temperature: xr.DataArray,
    surface_pressure: xr.DataArray,
    *,
    grid: GaussianGrid,
    horizontal_dims: tuple[str, ...],
    hybrid_dim: str,
    ntrunc: int,
) -> _PreparedPvFields:
    add_native_runtime_dir()
    from fullpos import _ectrans

    leading_dims = [dim for dim in temperature.dims if dim != hybrid_dim and dim not in horizontal_dims]
    horizontal_shape = tuple(temperature.sizes[dim] for dim in horizontal_dims)
    leading_shape = tuple(temperature.sizes[dim] for dim in leading_dims)
    nlead = int(np.prod(leading_shape, dtype=np.int64)) if leading_shape else 1
    nlev = temperature.sizes[hybrid_dim]
    point_count = int(np.prod(horizontal_shape, dtype=np.int64))
    if point_count != grid.size:
        raise ValueError(f"horizontal dimensions contain {point_count} points, but {grid.name} expects {grid.size}")

    u_flat = _flatten_hybrid_grid(u, leading_dims=leading_dims, hybrid_dim=hybrid_dim, horizontal_dims=horizontal_dims)
    v_flat = _flatten_hybrid_grid(v, leading_dims=leading_dims, hybrid_dim=hybrid_dim, horizontal_dims=horizontal_dims)
    t_flat = _flatten_hybrid_grid(temperature, leading_dims=leading_dims, hybrid_dim=hybrid_dim, horizontal_dims=horizontal_dims)
    sp_flat = _flatten_surface_grid(surface_pressure, leading_dims=leading_dims, horizontal_dims=horizontal_dims)
    if np.any(sp_flat <= 0.0):
        raise ValueError("surface_pressure must be positive to derive ln(surface_pressure) gradients")
    scalars = np.concatenate([t_flat, np.log(sp_flat)], axis=0)

    vort_flat, _div_flat, scalar_ns_flat, scalar_ew_flat = _ectrans.vector_scalar_diagnostics(
        np.ascontiguousarray(u_flat, dtype=np.float64),
        np.ascontiguousarray(v_flat, dtype=np.float64),
        np.ascontiguousarray(scalars, dtype=np.float64),
        _native_pl(grid),
        int(ntrunc),
    )
    scalar_count_3d = nlead * nlev
    tm_flat = np.asarray(scalar_ns_flat[:scalar_count_3d], dtype=np.float64)
    tl_flat = np.asarray(scalar_ew_flat[:scalar_count_3d], dtype=np.float64)
    spm_flat = np.asarray(scalar_ns_flat[scalar_count_3d:], dtype=np.float64) * sp_flat
    spl_flat = np.asarray(scalar_ew_flat[scalar_count_3d:], dtype=np.float64) * sp_flat

    return _PreparedPvFields(
        relative_vorticity=_unflatten_hybrid_grid(
            np.asarray(vort_flat, dtype=np.float64),
            template=temperature,
            leading_dims=leading_dims,
            hybrid_dim=hybrid_dim,
            horizontal_dims=horizontal_dims,
            leading_shape=leading_shape,
            horizontal_shape=horizontal_shape,
        ),
        temperature_meridional_gradient=_unflatten_hybrid_grid(
            tm_flat,
            template=temperature,
            leading_dims=leading_dims,
            hybrid_dim=hybrid_dim,
            horizontal_dims=horizontal_dims,
            leading_shape=leading_shape,
            horizontal_shape=horizontal_shape,
        ),
        temperature_zonal_gradient=_unflatten_hybrid_grid(
            tl_flat,
            template=temperature,
            leading_dims=leading_dims,
            hybrid_dim=hybrid_dim,
            horizontal_dims=horizontal_dims,
            leading_shape=leading_shape,
            horizontal_shape=horizontal_shape,
        ),
        surface_pressure_meridional_gradient=_unflatten_surface_grid(
            spm_flat,
            template=surface_pressure,
            leading_dims=leading_dims,
            horizontal_dims=horizontal_dims,
            leading_shape=leading_shape,
            horizontal_shape=horizontal_shape,
        ),
        surface_pressure_zonal_gradient=_unflatten_surface_grid(
            spl_flat,
            template=surface_pressure,
            leading_dims=leading_dims,
            horizontal_dims=horizontal_dims,
            leading_shape=leading_shape,
            horizontal_shape=horizontal_shape,
        ),
    )


def _resolve_kappa(
    kappa: xr.DataArray | None,
    *,
    specific_humidity: xr.DataArray | None,
    reference: xr.DataArray,
    hybrid_dim: str,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> xr.DataArray:
    if kappa is not None:
        if not isinstance(kappa, xr.DataArray):
            raise TypeError("kappa must be an xarray.DataArray")
        _validate_pressure_data_array(kappa, chunks=None, variable_name="kappa")
        if kappa.dims != reference.dims or kappa.shape != reference.shape:
            raise ValueError("kappa must have the same dimensions and shape as temperature")
        _validate_same_coords(reference, kappa)
        return kappa
    if specific_humidity is None:
        raise ValueError("diagnose_potential_vorticity requires either kappa=... or specific_humidity=... for native FULLPOS GPRCP")
    if not isinstance(specific_humidity, xr.DataArray):
        raise TypeError("specific_humidity must be an xarray.DataArray")
    _validate_pressure_data_array(specific_humidity, chunks=chunks, variable_name="specific_humidity")
    if specific_humidity.dims != reference.dims or specific_humidity.shape != reference.shape:
        raise ValueError("specific_humidity must have the same dimensions and shape as temperature")
    _validate_same_coords(reference, specific_humidity)
    return _compute_gprcp_kappa(specific_humidity, hybrid_dim=hybrid_dim, chunks=chunks, keep_attrs=keep_attrs)


def _compute_gprcp_kappa(
    specific_humidity: xr.DataArray,
    *,
    hybrid_dim: str,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
) -> xr.DataArray:
    add_native_runtime_dir()
    from fullpos import _vertical_native

    values = np.empty(specific_humidity.shape, dtype=np.float64)
    for selection in _leading_chunk_selections(specific_humidity, hybrid_dim=hybrid_dim, chunks=chunks):
        block = specific_humidity.isel(selection) if selection else specific_humidity
        native = _vertical_native.gprcp_kappa(_flatten_hybrid_columns(block, hybrid_dim))
        target = _output_selection(selection, specific_humidity, hybrid_dim=hybrid_dim, nlevels=specific_humidity.sizes[hybrid_dim])
        values[target] = _unflatten_output_columns(native, block, hybrid_dim, specific_humidity.sizes[hybrid_dim])
    return _wrap_prepared_like(
        specific_humidity,
        values,
        name="kappa",
        diagnostic="kappa",
        keep_attrs=keep_attrs,
        native_path="GPRCP",
    )


def _resolve_grid(obj: xr.DataArray, *, source_grid: str | GaussianGrid | None) -> GaussianGrid:
    if isinstance(source_grid, GaussianGrid):
        return source_grid
    if source_grid is not None:
        return parse_grid(str(source_grid))
    try:
        return parse_grid(infer_grid_name_from_attrs(obj.attrs))
    except ValueError:
        return parse_grid(infer_grid_name_from_shape(obj.sizes, obj.dims))


def _horizontal_dims(obj: xr.DataArray, grid: GaussianGrid) -> tuple[str, ...]:
    if grid.is_reduced:
        return (_find_packed_dim(obj, grid),)
    if "latitude" not in obj.dims or "longitude" not in obj.dims:
        raise ValueError("regular Gaussian PV diagnostics require latitude and longitude dimensions")
    if obj.sizes["latitude"] != grid.nlat or obj.sizes["longitude"] != grid.work_nlon:
        raise ValueError(
            f"regular {grid.name} expects shape ({grid.nlat}, {grid.work_nlon}) "
            "on latitude/longitude dimensions"
        )
    return ("latitude", "longitude")


def _find_packed_dim(obj: xr.DataArray, grid: GaussianGrid) -> str:
    for dim in obj.dims:
        if obj.sizes[dim] == grid.size:
            return dim
    raise ValueError(f"could not find packed reduced dimension with size {grid.size}")


def _auto_coriolis(
    reference: xr.DataArray,
    *,
    grid: GaussianGrid,
    horizontal_dims: tuple[str, ...],
    hybrid_dim: str,
) -> xr.DataArray:
    omega = 7.292115e-5
    lats = gaussian_latitudes(grid.nlat)
    if grid.is_reduced:
        assert grid.pl is not None
        packed_dim = horizontal_dims[0]
        lat_values = np.repeat(lats, np.asarray(grid.pl, dtype=np.int64))
        coriolis = xr.DataArray(
            2.0 * omega * np.sin(np.deg2rad(lat_values)),
            dims=(packed_dim,),
            coords={packed_dim: reference.coords[packed_dim]} if packed_dim in reference.coords else None,
            name="coriolis",
        )
    else:
        coriolis = xr.DataArray(
            2.0 * omega * np.sin(np.deg2rad(lats)),
            dims=("latitude",),
            coords={"latitude": reference.coords["latitude"]} if "latitude" in reference.coords else None,
            name="coriolis",
        )
    return _broadcast_surface_vector(coriolis, reference, hybrid_dim=hybrid_dim, name="coriolis")


def _normalize_ntrunc(grid: GaussianGrid, ntrunc: int | None) -> int:
    max_ntrunc = grid.n - 1
    if ntrunc is None:
        return max_ntrunc
    value = int(ntrunc)
    if value < 0 or value > max_ntrunc:
        raise ValueError(f"ntrunc must be between 0 and {max_ntrunc}, got {ntrunc}")
    return value


def _native_pl(grid: GaussianGrid) -> np.ndarray:
    if grid.pl is None:
        return np.full(grid.nlat, grid.work_nlon, dtype=np.int32)
    return np.asarray(grid.pl, dtype=np.int32)


def _validate_spectral_chunks(
    obj: xr.DataArray,
    *,
    hybrid_dim: str,
    horizontal_dims: tuple[str, ...],
    chunks: dict[str, int] | None,
) -> None:
    if chunks is None:
        return
    valid = set(obj.dims) - {hybrid_dim, *horizontal_dims}
    invalid = [dim for dim in chunks if dim not in valid]
    if invalid:
        raise ValueError(
            "PV ECTRANS diagnostics require complete vertical columns and complete horizontal grids; "
            f"chunks may only target leading dimensions {sorted(valid)}, got invalid keys {invalid}"
        )


def _spectral_chunk_selections(
    obj: xr.DataArray,
    *,
    hybrid_dim: str,
    horizontal_dims: tuple[str, ...],
    chunks: dict[str, int] | None,
):
    leading_dims = [dim for dim in obj.dims if dim != hybrid_dim and dim not in horizontal_dims]
    if not leading_dims:
        yield {}
        return
    selections = [{}]
    for dim in leading_dims:
        size = obj.sizes[dim]
        step = size if chunks is None or dim not in chunks else int(chunks[dim])
        if step <= 0:
            raise ValueError(f"chunks[{dim!r}] must be a positive integer")
        slices = [slice(start, min(start + step, size)) for start in range(0, size, step)]
        selections = [dict(sel, **{dim: slc}) for sel in selections for slc in slices]
    for selection in selections:
        yield selection


def _surface_selection(obj: xr.DataArray, selection: dict[str, slice]) -> xr.DataArray:
    surface_selection = {dim: slc for dim, slc in selection.items() if dim in obj.dims}
    return obj.isel(surface_selection) if surface_selection else obj


def _array_selection(selection: dict[str, slice], obj: xr.DataArray) -> tuple:
    return tuple(selection.get(dim, slice(None)) for dim in obj.dims)


def _flatten_hybrid_grid(
    obj: xr.DataArray,
    *,
    leading_dims: list[str],
    hybrid_dim: str,
    horizontal_dims: tuple[str, ...],
) -> np.ndarray:
    transposed = obj.transpose(*leading_dims, hybrid_dim, *horizontal_dims)
    nlev = obj.sizes[hybrid_dim]
    point_count = int(np.prod([obj.sizes[dim] for dim in horizontal_dims], dtype=np.int64))
    arr = np.asarray(transposed.values, dtype=np.float64)
    return np.ascontiguousarray(arr.reshape((-1, nlev, point_count)).reshape((-1, point_count)))


def _flatten_surface_grid(
    obj: xr.DataArray,
    *,
    leading_dims: list[str],
    horizontal_dims: tuple[str, ...],
) -> np.ndarray:
    transposed = obj.transpose(*leading_dims, *horizontal_dims)
    point_count = int(np.prod([obj.sizes[dim] for dim in horizontal_dims], dtype=np.int64))
    arr = np.asarray(transposed.values, dtype=np.float64)
    return np.ascontiguousarray(arr.reshape((-1, point_count)))


def _unflatten_hybrid_grid(
    values: np.ndarray,
    *,
    template: xr.DataArray,
    leading_dims: list[str],
    hybrid_dim: str,
    horizontal_dims: tuple[str, ...],
    leading_shape: tuple[int, ...],
    horizontal_shape: tuple[int, ...],
) -> xr.DataArray:
    nlev = template.sizes[hybrid_dim]
    arr = np.asarray(values, dtype=np.float64).reshape(leading_shape + (nlev,) + horizontal_shape)
    dims = tuple(leading_dims) + (hybrid_dim,) + tuple(horizontal_dims)
    return xr.DataArray(arr, dims=dims).transpose(*template.dims)


def _unflatten_surface_grid(
    values: np.ndarray,
    *,
    template: xr.DataArray,
    leading_dims: list[str],
    horizontal_dims: tuple[str, ...],
    leading_shape: tuple[int, ...],
    horizontal_shape: tuple[int, ...],
) -> xr.DataArray:
    arr = np.asarray(values, dtype=np.float64).reshape(leading_shape + horizontal_shape)
    dims = tuple(leading_dims) + tuple(horizontal_dims)
    return xr.DataArray(arr, dims=dims).transpose(*template.dims)


def _wrap_prepared_like(
    template: xr.DataArray,
    values: np.ndarray,
    *,
    name: str,
    diagnostic: str,
    keep_attrs: bool,
    native_path: str = "ECTRANS",
) -> xr.DataArray:
    coords = {dim: template.coords[dim].values for dim in template.dims if dim in template.coords and template.coords[dim].ndim == 1}
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["diagnostic"] = diagnostic
    attrs["vertical_backend"] = "FULLPOS"
    attrs["vertical_native_path"] = native_path
    if keep_attrs:
        attrs = append_history(attrs, f"fullpos diagnose_potential_vorticity prepare={diagnostic}")
    return xr.DataArray(values, dims=template.dims, coords=coords, attrs=attrs, name=name)


def _broadcast_surface_vector(
    obj: xr.DataArray,
    reference: xr.DataArray,
    *,
    hybrid_dim: str,
    name: str,
) -> xr.DataArray:
    if not isinstance(obj, xr.DataArray):
        raise TypeError(f"{name} must be an xarray.DataArray")
    if hybrid_dim in obj.dims:
        raise ValueError(f"{name} must not include the hybrid dimension")
    base = reference.isel({hybrid_dim: 0}, drop=True)
    try:
        return obj.broadcast_like(base)
    except Exception as exc:
        raise ValueError(f"{name} dimensions must be broadcastable to the input horizontal grid") from exc


def _wrap_diagnostic_output(
    template: xr.DataArray,
    values: np.ndarray,
    *,
    name: str,
    attrs_name: str,
    keep_attrs: bool,
) -> xr.DataArray:
    coords = {dim: template.coords[dim].values for dim in template.dims if dim in template.coords and template.coords[dim].ndim == 1}
    attrs = dict(template.attrs) if keep_attrs else {}
    attrs["diagnostic"] = attrs_name
    attrs["vertical_backend"] = "FULLPOS"
    attrs["vertical_native_path"] = "GPPVO"
    if keep_attrs:
        attrs = append_history(attrs, f"fullpos diagnose_potential_vorticity output={attrs_name}")
    return xr.DataArray(values, dims=template.dims, coords=coords, attrs=attrs, name=name)
