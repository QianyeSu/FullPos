from __future__ import annotations

from datetime import datetime

import numpy as np
import xarray as xr

from .grids import (
    GaussianGrid,
    classic_reduced_grid_from_pl,
    infer_grid_from_attrs,
    infer_grid_name_from_shape,
    parse_grid,
)


def infer_source_grid(obj) -> str | GaussianGrid:
    """Infer the source grid from xarray metadata or dimensions."""
    if isinstance(obj, xr.DataArray):
        try:
            return infer_grid_from_attrs(obj.attrs)
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


def resolve_source_grid(obj, source_grid=None) -> str | GaussianGrid:
    """Resolve an explicit or inferred source grid, preserving packed metadata."""
    if isinstance(source_grid, GaussianGrid):
        return source_grid
    if source_grid is None:
        return infer_source_grid(obj)
    text = str(source_grid).strip()
    if text.upper().startswith("N"):
        try:
            return infer_source_grid(obj)
        except ValueError:
            pass
    return parse_grid(text)


def resolve_target_grid(obj, target_grid) -> str | GaussianGrid:
    """Resolve a target grid, preserving packed classic N-grid context when available."""
    if isinstance(target_grid, GaussianGrid):
        return target_grid
    text = str(target_grid).strip()
    if text.upper().startswith("N"):
        n = _resolution_from_name(text)
        attrs = getattr(obj, "attrs", {})
        grid = _classic_reduced_grid_from_current_attrs(attrs, n)
        if grid is not None:
            return grid
        grid = _classic_reduced_grid_from_preserved_source_attrs(attrs, n)
        if grid is not None:
            return grid
    return parse_grid(text)


def preserve_source_grid_attrs(attrs, source_grid: GaussianGrid) -> dict:
    """Preserve classic reduced source geometry for later ``target_grid='N*'`` calls."""
    out = dict(attrs)
    if source_grid.kind == "classic_reduced" and source_grid.pl is not None:
        out["fullpos_source_grid"] = source_grid.name
        out["fullpos_source_grid_kind"] = source_grid.kind
        out["fullpos_source_GRIB_N"] = source_grid.n
        out["fullpos_source_GRIB_gridType"] = "reduced_gg"
        out["fullpos_source_GRIB_numberOfPoints"] = source_grid.size
        out["fullpos_source_GRIB_pl"] = np.asarray(source_grid.pl, dtype=np.int64)
    return out


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


def _resolution_from_name(name: str) -> int | None:
    text = str(name).strip().upper()
    if len(text) < 2:
        return None
    try:
        return int(text[1:])
    except ValueError:
        return None


def _classic_reduced_grid_from_current_attrs(attrs, n: int | None) -> GaussianGrid | None:
    if n is None:
        return None
    try:
        grid = infer_grid_from_attrs(attrs)
    except ValueError:
        return None
    if (
        isinstance(grid, GaussianGrid)
        and grid.kind == "classic_reduced"
        and grid.n == int(n)
    ):
        return grid
    return None


def _classic_reduced_grid_from_preserved_source_attrs(
    attrs,
    n: int | None,
) -> GaussianGrid | None:
    if n is None or not attrs:
        return None
    if attrs.get("fullpos_source_grid_kind") != "classic_reduced":
        return None
    source_n = attrs.get("fullpos_source_GRIB_N")
    source_pl = attrs.get("fullpos_source_GRIB_pl")
    if source_n is None or source_pl is None:
        return None
    source_n = int(np.asarray(source_n).reshape(-1)[0])
    if source_n != int(n):
        return None
    return classic_reduced_grid_from_pl(source_n, source_pl)
