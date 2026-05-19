from __future__ import annotations

import numpy as np
import xarray as xr

from ..grids import GaussianGrid, gaussian_latitudes, regular_longitudes
from ..metadata import append_history, format_regrid_history, preserve_source_grid_attrs
from .native import spectral_regrid_chunks


def regrid_data_array(
    obj: xr.DataArray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    ntrunc: int | None = None,
    chunk_size: int | None = 64,
    missing_policy: str = "error",
    keep_attrs: bool = True,
) -> xr.DataArray:
    """Regrid an ``xarray.DataArray`` between supported Gaussian grids."""
    if source_grid.is_reduced:
        horizontal_dim = _find_packed_dim(obj, source_grid)
        leading_dims = [dim for dim in obj.dims if dim != horizontal_dim]
    else:
        horizontal_dim = None
        leading_dims = [dim for dim in obj.dims if dim not in {"latitude", "longitude"}]
        if "latitude" not in obj.dims or "longitude" not in obj.dims:
            raise ValueError("regular Gaussian input must have latitude and longitude dims")

    transposed = obj.transpose(*leading_dims, *(_horizontal_dims(obj, source_grid)))
    leading_shape = transposed.shape[: len(leading_dims)]
    flat = np.asarray(transposed.values).reshape((-1, source_grid.size))
    batched = spectral_regrid_chunks(
        flat,
        source_grid=source_grid,
        target_grid=target_grid,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
    )

    if target_grid.is_reduced:
        stacked = batched.reshape(leading_shape + (target_grid.size,))
    else:
        stacked = batched.reshape(leading_shape + (target_grid.nlat, target_grid.work_nlon))

    if target_grid.is_reduced:
        out_dims = tuple(leading_dims) + ("values",)
        coords = _leading_coords(obj, leading_dims)
        coords["values"] = np.arange(target_grid.size)
        coords["latitude"] = ("values", _packed_latitudes(target_grid))
        coords["longitude"] = ("values", _packed_longitudes(target_grid))
    else:
        out_dims = tuple(leading_dims) + ("latitude", "longitude")
        coords = _leading_coords(obj, leading_dims)
        coords["latitude"] = gaussian_latitudes(target_grid.nlat)
        coords["longitude"] = regular_longitudes(target_grid.work_nlon)

    return xr.DataArray(
        stacked,
        dims=out_dims,
        coords=coords,
        name=obj.name,
        attrs=_output_attrs(
            obj.attrs,
            source_grid,
            target_grid,
            keep_attrs=keep_attrs,
            ntrunc=ntrunc,
            chunk_size=chunk_size,
        ),
    )


def _horizontal_dims(obj: xr.DataArray, grid: GaussianGrid) -> tuple[str, ...]:
    if grid.is_reduced:
        return (_find_packed_dim(obj, grid),)
    return ("latitude", "longitude")


def _find_packed_dim(obj: xr.DataArray, grid: GaussianGrid) -> str:
    for dim in obj.dims:
        if obj.sizes[dim] == grid.size:
            return dim
    raise ValueError(f"could not find packed reduced dimension with size {grid.size}")


def _leading_coords(obj: xr.DataArray, leading_dims: list[str]) -> dict:
    coords = {}
    for name, coord in obj.coords.items():
        if set(coord.dims).issubset(set(leading_dims)):
            coords[name] = coord
    return coords


def _packed_latitudes(grid: GaussianGrid) -> np.ndarray:
    assert grid.pl is not None
    lats = gaussian_latitudes(grid.nlat)
    return np.repeat(lats, np.asarray(grid.pl, dtype=np.int64))


def _packed_longitudes(grid: GaussianGrid) -> np.ndarray:
    assert grid.pl is not None
    rows = []
    for row_nlon in grid.pl:
        rows.append(regular_longitudes(int(row_nlon)))
    return np.concatenate(rows)


def _output_attrs(
    attrs,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    *,
    keep_attrs: bool,
    ntrunc: int | None,
    chunk_size: int | None,
) -> dict:
    out = dict(attrs) if keep_attrs else {}
    out = preserve_source_grid_attrs(out, source_grid)
    out["GRIB_N"] = target_grid.n
    out["GRIB_gridType"] = "reduced_gg" if target_grid.is_reduced else "regular_gg"
    out["GRIB_numberOfPoints"] = target_grid.size
    if target_grid.is_reduced:
        out["GRIB_pl"] = np.asarray(target_grid.pl, dtype=np.int64)
    else:
        out.pop("GRIB_pl", None)
    out = append_history(
        out,
        format_regrid_history(
            source_grid=source_grid,
            target_grid=target_grid,
            method="linear",
            ntrunc=ntrunc,
            chunk_size=chunk_size,
        ),
    )
    return out
