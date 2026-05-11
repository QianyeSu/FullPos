from __future__ import annotations

from datetime import datetime

import xarray as xr

from .grids import (
    infer_grid_name_from_attrs,
    infer_grid_name_from_shape,
)


def infer_source_grid(obj) -> str:
    """Infer the source grid name from xarray metadata or dimensions."""
    if isinstance(obj, xr.DataArray):
        try:
            return infer_grid_name_from_attrs(obj.attrs)
        except ValueError as exc:
            if has_grid_metadata(obj.attrs):
                raise
            try:
                return infer_grid_name_from_shape(obj.sizes, obj.dims)
            except ValueError:
                raise exc
    if isinstance(obj, xr.Dataset):
        for data_array in obj.data_vars.values():
            return infer_source_grid(data_array)
    raise ValueError("source_grid is required when it cannot be inferred from GRIB attrs")


def has_grid_metadata(attrs) -> bool:
    """Return whether attributes contain enough grid hints to attempt inference."""
    keys = {
        "GRIB_gridType",
        "gridType",
        "GRIB_N",
        "N",
        "GRIB_pl",
        "pl",
        "GRIB_numberOfPoints",
        "numberOfPoints",
    }
    return any(key in attrs for key in keys)


def append_history(attrs, entry: str) -> dict:
    """Prepend a CDO-style history entry while preserving existing history."""
    out = dict(attrs)
    timestamp = datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    line = f"{timestamp}: {entry}"
    previous = out.get("history")
    out["history"] = f"{line}\n{previous}" if previous else line
    return out


def format_regrid_history(
    *,
    source_grid,
    target_grid,
    method: str,
    ntrunc=None,
    chunk_size=64,
    variables=None,
) -> str:
    """Format a compact history string for regridding/filtering operations."""
    parts = [
        "fullpos regrid",
        f"source_grid={_grid_name(source_grid)}",
        f"target_grid={_grid_name(target_grid)}",
        f"method={method}",
        f"chunk_size={chunk_size}",
    ]
    if ntrunc is not None:
        parts.append(f"ntrunc={ntrunc}")
    if variables is not None:
        parts.append("variables=" + ",".join(str(v) for v in variables))
    return " ".join(parts)


def _grid_name(grid) -> str:
    return getattr(grid, "name", str(grid))
