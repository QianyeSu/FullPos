from __future__ import annotations

import numpy as np
import xarray as xr

from .grids import GaussianGrid, gaussian_latitudes, parse_grid, regular_longitudes


def masked_surface_interpolate(
    field,
    *,
    land_sea_mask: xr.DataArray,
    source_grid: str | GaussianGrid,
    target_grid: str | GaussianGrid,
    kind: str = "sea",
    method: str = "average",
    average_radius: int = 1,
    sea_threshold: float = 0.5,
    chunks: dict[str, int] | None = None,
    keep_attrs: bool = True,
):
    """Interpolate a masked surface field with native FULLPOS halo kernels.

    This is a convenience wrapper for bitmap fields such as SST. It samples a
    regular latitude/longitude land-sea mask to the source and target Gaussian
    grids, then calls ``horizontal_interpolate`` with ``source_mask`` and
    ``target_mask``. The meteorological field interpolation is performed by
    native FULLPOS ``FPAVG`` or ``FPNEAR``.

    ``True`` in the generated masks means "valid for this field". With ERA5
    land-sea mask convention, ``kind="sea"`` uses ``lsm < 0.5`` and is the
    usual choice for SST; ``kind="land"`` uses ``lsm >= 0.5``.
    """
    from .interpolation.horizontal import horizontal_interpolate

    source_mask = land_sea_mask_to_grid(
        land_sea_mask,
        target_grid=source_grid,
        sea_threshold=sea_threshold,
        kind=kind,
    )
    target_mask = land_sea_mask_to_grid(
        land_sea_mask,
        target_grid=target_grid,
        sea_threshold=sea_threshold,
        kind=kind,
    )
    out = horizontal_interpolate(
        field,
        source_grid=source_grid,
        target_grid=target_grid,
        method=method,
        average_radius=average_radius,
        source_mask=source_mask,
        target_mask=target_mask,
        chunks=chunks,
        keep_attrs=keep_attrs,
    )
    out.attrs.update(
        {
            "fullpos_surface_mask_kind": str(kind).lower(),
            "fullpos_surface_mask_threshold": float(sea_threshold),
            "fullpos_surface_mask_source_grid": _grid_name(source_grid),
            "fullpos_surface_mask_target_grid": _grid_name(target_grid),
        }
    )
    return out


def land_sea_mask_to_grid(
    lsm: xr.DataArray,
    *,
    target_grid: str | GaussianGrid,
    sea_threshold: float = 0.5,
    kind: str = "sea",
) -> xr.DataArray:
    """Sample a regular latitude/longitude land-sea mask to a Gaussian grid.

    The helper prepares boolean masks for native FULLPOS horizontal
    interpolation. It does not interpolate meteorological fields. ``True``
    marks source or target points that should participate in interpolation.

    Parameters
    ----------
    lsm:
        Land-sea mask on a regular latitude/longitude grid. ERA5 convention is
        ``0`` over sea and ``1`` over land.
    target_grid:
        Gaussian grid such as ``"O96"`` or ``"F480"``.
    sea_threshold:
        Values below this threshold are treated as sea.
    kind:
        ``"sea"`` returns sea points, and ``"land"`` returns land points.
    """
    grid = target_grid if isinstance(target_grid, GaussianGrid) else parse_grid(target_grid)
    target_lats, target_lons, dims, coords = _target_points(grid)
    values = _sample_regular_ll_nearest(lsm, target_lats, target_lons)
    normalized = str(kind).lower()
    if normalized == "sea":
        mask = values < float(sea_threshold)
    elif normalized == "land":
        mask = values >= float(sea_threshold)
    else:
        raise ValueError("kind must be 'sea' or 'land'")
    return xr.DataArray(mask, dims=dims, coords=coords, name=f"{normalized}_mask")


def _grid_name(grid: str | GaussianGrid) -> str:
    return grid.name if isinstance(grid, GaussianGrid) else str(grid)


def _target_points(grid: GaussianGrid):
    lats = gaussian_latitudes(grid.nlat)
    if grid.is_reduced:
        pl = np.asarray(grid.pl, dtype=np.int64)
        lat_values = np.repeat(lats, pl)
        lon_values = np.concatenate([regular_longitudes(int(nlon)) for nlon in pl])
        coords = {
            "values": np.arange(grid.size),
            "latitude": ("values", lat_values),
            "longitude": ("values", lon_values),
        }
        return lat_values, lon_values, ("values",), coords

    lons = regular_longitudes(grid.work_nlon)
    lat_2d, lon_2d = np.meshgrid(lats, lons, indexing="ij")
    coords = {"latitude": lats, "longitude": lons}
    return lat_2d, lon_2d, ("latitude", "longitude"), coords


def _sample_regular_ll_nearest(lsm: xr.DataArray, target_lats: np.ndarray, target_lons: np.ndarray) -> np.ndarray:
    if "latitude" not in lsm.dims or "longitude" not in lsm.dims:
        raise ValueError("lsm must have latitude and longitude dimensions")
    source = lsm.transpose("latitude", "longitude")
    src_lats = np.asarray(source["latitude"].values, dtype=np.float64)
    src_lons = np.asarray(source["longitude"].values, dtype=np.float64)
    values = np.asarray(source.values, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError("lsm must be a 2D DataArray")

    lat_idx = _nearest_lat_indices(src_lats, np.asarray(target_lats, dtype=np.float64).reshape(-1))
    lon_idx = _nearest_periodic_lon_indices(src_lons, np.asarray(target_lons, dtype=np.float64).reshape(-1))
    sampled = values[lat_idx, lon_idx]
    return sampled.reshape(np.asarray(target_lats).shape)


def _nearest_lat_indices(source_lats: np.ndarray, target_lats: np.ndarray) -> np.ndarray:
    ascending = source_lats[0] < source_lats[-1]
    work_lats = source_lats if ascending else source_lats[::-1]
    idx = np.searchsorted(work_lats, target_lats, side="left")
    idx = np.clip(idx, 1, work_lats.size - 1)
    left = idx - 1
    right = idx
    choose_right = np.abs(work_lats[right] - target_lats) < np.abs(work_lats[left] - target_lats)
    out = np.where(choose_right, right, left)
    if not ascending:
        out = source_lats.size - 1 - out
    return out.astype(np.int64)


def _nearest_periodic_lon_indices(source_lons: np.ndarray, target_lons: np.ndarray) -> np.ndarray:
    lons = np.mod(source_lons, 360.0)
    order = np.argsort(lons)
    sorted_lons = lons[order]
    targets = np.mod(target_lons, 360.0)
    idx = np.searchsorted(sorted_lons, targets, side="left")
    left = (idx - 1) % sorted_lons.size
    right = idx % sorted_lons.size
    left_dist = np.abs(((targets - sorted_lons[left] + 180.0) % 360.0) - 180.0)
    right_dist = np.abs(((targets - sorted_lons[right] + 180.0) % 360.0) - 180.0)
    chosen = np.where(right_dist < left_dist, right, left)
    return order[chosen].astype(np.int64)
