from __future__ import annotations

import xarray as xr

from ..grids import parse_grid
from ..metadata import append_history, format_regrid_history, infer_source_grid
from .regridder import DEFAULT_CHUNK_SIZE, Regridder


_HORIZONTAL_DIM_NAMES = {"values", "latitude", "longitude"}


def regrid(
    obj,
    *,
    source_grid=None,
    target_grid=None,
    method="linear",
    variables=None,
    missing_value=None,
    ntrunc=None,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    keep_attrs=True,
    skip_non_horizontal=True,
):
    """Regrid xarray data between supported ECMWF Gaussian grids.

    Parameters
    ----------
    obj:
        Input ``xarray.DataArray`` or ``xarray.Dataset``. Horizontal dimensions
        are inferred from GRIB metadata when possible, or from explicit
        ``source_grid``.
    source_grid, target_grid:
        Grid names such as ``"O320"``, ``"O480"``, ``"N320"``, or parsed grid
        objects. ``target_grid`` is required.
    method:
        ``"linear"``/``"spectral"`` uses the native ECTRANS spectral path.
        ``"masked"`` uses the NaN-aware surface helper for fields with bitmaps.
    variables:
        Optional dataset variable subset. Non-horizontal variables are skipped
        by default when ``variables`` is not explicit.
    ntrunc:
        Optional triangular truncation used by the spectral backend.
    chunk_size:
        Number of flattened leading fields transformed per native backend call.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        Regridded object with updated Gaussian-grid metadata and appended
        ``history``.
    """
    if target_grid is None:
        raise ValueError("target_grid is required")
    if method not in {"linear", "spectral", "masked"}:
        raise ValueError("method must be 'linear', 'spectral', or 'masked'")

    if isinstance(obj, xr.DataArray):
        regridder = Regridder(
            source_grid or _infer_source_grid(obj),
            target_grid,
            ntrunc=ntrunc,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
            method=method,
            missing_value=missing_value,
        )
        return regridder.regrid_data_array(obj, keep_attrs=keep_attrs)

    if isinstance(obj, xr.Dataset):
        selected = list(obj.data_vars) if variables is None else [str(v) for v in variables]
        missing = [name for name in selected if name not in obj.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        out = xr.Dataset(attrs=dict(obj.attrs))
        _copy_safe_coords(obj, out)
        for name in selected:
            data_array = obj[name]
            try:
                variable_source_grid = parse_grid(source_grid or infer_source_grid(data_array))
            except ValueError as exc:
                if variables is None and skip_non_horizontal:
                    continue
                raise ValueError(f"cannot infer source grid for variable {name!r}") from exc
            if not _has_horizontal_dims(data_array, variable_source_grid):
                if variables is None and skip_non_horizontal:
                    continue
                raise ValueError(
                    f"variable {name!r} does not have horizontal dimensions for "
                    f"source grid {variable_source_grid!r}"
                )
            try:
                parsed_target = parse_grid(target_grid)
                regridder = Regridder(
                    variable_source_grid,
                    parsed_target,
                    ntrunc=ntrunc,
                    chunk_size=chunk_size,
                    missing_policy=missing_policy,
                    method=method,
                    missing_value=missing_value,
                )
                out[name] = regridder.regrid_data_array(data_array, keep_attrs=keep_attrs)
            except Exception as exc:
                raise type(exc)(
                    f"failed to regrid variable {name!r} "
                    f"from {variable_source_grid!r} to {target_grid!r}: {exc}"
                ) from exc
        if variables is None and skip_non_horizontal:
            for name, data_array in obj.data_vars.items():
                if name not in out.data_vars:
                    out[name] = data_array
        out.attrs = append_history(
            out.attrs,
            format_regrid_history(
                source_grid=source_grid or "inferred",
                target_grid=target_grid,
                method=method,
                ntrunc=ntrunc,
                chunk_size=chunk_size,
                variables=selected,
            ),
        )
        return out

    raise TypeError("obj must be an xarray DataArray or Dataset")


def regrid_values(
    values,
    *,
    source_grid,
    target_grid,
    ntrunc=None,
    axis=-1,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    method="linear",
    missing_value=None,
):
    """Regrid NumPy-like values between supported ECMWF Gaussian grids.

    ``axis`` identifies the horizontal dimension(s): a single packed reduced
    dimension for O-grids, or ``(lat_axis, lon_axis)`` for regular N-grids.
    Leading dimensions are treated as independent fields and processed in
    chunks.
    """
    if method not in {"linear", "spectral", "masked"}:
        raise ValueError("method must be 'linear', 'spectral', or 'masked'")
    regridder = Regridder(
        source_grid,
        target_grid,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        method=method,
        missing_value=missing_value,
    )
    return regridder.regrid_values(values, axis=axis)


def _infer_source_grid(obj) -> str:
    return infer_source_grid(obj)


def _copy_safe_coords(source: xr.Dataset, target: xr.Dataset) -> None:
    for name, coord in source.coords.items():
        if set(coord.dims).isdisjoint(_HORIZONTAL_DIM_NAMES):
            target.coords[name] = coord


def _has_horizontal_dims(obj: xr.DataArray, source_grid) -> bool:
    grid = source_grid if hasattr(source_grid, "is_reduced") else parse_grid(source_grid)
    if grid.is_reduced:
        return any(obj.sizes[dim] == grid.size for dim in obj.dims)
    return "latitude" in obj.dims and "longitude" in obj.dims
