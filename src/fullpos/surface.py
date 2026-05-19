from __future__ import annotations

import numpy as np
import xarray as xr

from .grids import GaussianGrid, gaussian_latitudes, regular_longitudes
from .metadata import append_history, format_regrid_history, preserve_source_grid_attrs


def masked_regrid_values(
    values,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    axis: int | tuple[int, ...] = -1,
    chunk_size: int | None = 64,
    missing_value=None,
) -> np.ndarray:
    """Regrid NumPy-like values with simple NaN-aware horizontal interpolation.

    This helper is intended for bitmap/missing-value surface fields where a
    global spectral transform would propagate NaNs. It is not a FULLPOS
    land/sea-mask interpolation implementation.
    """
    arr = np.asarray(values)
    expected_shape = _horizontal_shape(source_grid)
    horizontal_ndim = len(expected_shape)
    moved, leading_shape, original_axes = _move_horizontal_to_end(arr, expected_shape, axis)
    flat = moved.reshape((-1, source_grid.size))
    out = masked_regrid_chunks(
        flat,
        source_grid=source_grid,
        target_grid=target_grid,
        chunk_size=chunk_size,
        missing_value=missing_value,
    )
    target_shape = _horizontal_shape(target_grid)
    result = out.reshape(leading_shape + target_shape)
    if original_axes is None:
        return result
    return _move_horizontal_from_end(result, target_shape, original_axes, horizontal_ndim)


def masked_regrid_data_array(
    obj: xr.DataArray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    chunk_size: int | None = 64,
    missing_value=None,
    keep_attrs: bool = True,
) -> xr.DataArray:
    """Regrid an xarray DataArray with NaN-aware horizontal interpolation."""
    if source_grid.is_reduced:
        horizontal_dims = (_find_packed_dim(obj, source_grid),)
        leading_dims = [dim for dim in obj.dims if dim not in horizontal_dims]
    else:
        horizontal_dims = ("latitude", "longitude")
        leading_dims = [dim for dim in obj.dims if dim not in horizontal_dims]
        if "latitude" not in obj.dims or "longitude" not in obj.dims:
            raise ValueError("regular Gaussian input must have latitude and longitude dims")

    transposed = obj.transpose(*leading_dims, *horizontal_dims)
    leading_shape = transposed.shape[: len(leading_dims)]
    flat = np.asarray(transposed.values).reshape((-1, source_grid.size))
    out = masked_regrid_chunks(
        flat,
        source_grid=source_grid,
        target_grid=target_grid,
        chunk_size=chunk_size,
        missing_value=missing_value,
    )
    target_shape = _horizontal_shape(target_grid)
    stacked = out.reshape(leading_shape + target_shape)

    if target_grid.is_reduced:
        dims = tuple(leading_dims) + ("values",)
        coords = _leading_coords(obj, leading_dims)
        coords["values"] = np.arange(target_grid.size)
        coords["latitude"] = ("values", _packed_latitudes(target_grid))
        coords["longitude"] = ("values", _packed_longitudes(target_grid))
    else:
        dims = tuple(leading_dims) + ("latitude", "longitude")
        coords = _leading_coords(obj, leading_dims)
        coords["latitude"] = gaussian_latitudes(target_grid.nlat)
        coords["longitude"] = regular_longitudes(target_grid.work_nlon)

    return xr.DataArray(
        stacked,
        dims=dims,
        coords=coords,
        name=obj.name,
        attrs=_output_attrs(
            obj.attrs,
            source_grid,
            target_grid,
            keep_attrs=keep_attrs,
            method="masked",
            chunk_size=chunk_size,
        ),
    )


def masked_regrid_chunks(
    values: np.ndarray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    chunk_size: int | None = 64,
    missing_value=None,
) -> np.ndarray:
    """Chunk a flattened field batch for masked horizontal interpolation."""
    arr = np.asarray(values)
    if arr.ndim != 2:
        raise ValueError("chunked values must have shape (nfield, npoints)")
    if arr.shape[1] != source_grid.size:
        raise ValueError(f"batch fields have {arr.shape[1]} points, expected {source_grid.size}")
    if arr.shape[0] == 0:
        raise ValueError("cannot regrid an empty field batch")

    chunks = []
    for start, stop in _chunk_slices(arr.shape[0], chunk_size):
        chunks.append(
            _masked_regrid_batch(
                np.asarray(arr[start:stop], dtype=np.float64),
                source_grid=source_grid,
                target_grid=target_grid,
                missing_value=missing_value,
            )
        )
    return np.concatenate(chunks, axis=0)


def _masked_regrid_batch(
    values: np.ndarray,
    *,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    missing_value=None,
) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if missing_value is not None:
        arr = np.where(arr == float(missing_value), np.nan, arr)

    source_lats = gaussian_latitudes(source_grid.nlat)
    target_lats = gaussian_latitudes(target_grid.nlat)
    source_offsets = _row_offsets(source_grid)
    target_rows = []
    for target_lat, target_nlon in zip(target_lats, _row_lengths(target_grid)):
        target_lons = regular_longitudes(int(target_nlon))
        row = _interpolate_latitude_row(
            arr,
            source_grid=source_grid,
            source_lats=source_lats,
            source_offsets=source_offsets,
            target_lat=float(target_lat),
            target_lons=target_lons,
        )
        target_rows.append(row.astype(np.float32, copy=False))

    if target_grid.is_reduced:
        return np.concatenate(target_rows, axis=1)
    return np.stack(target_rows, axis=1)


def _interpolate_latitude_row(
    values: np.ndarray,
    *,
    source_grid: GaussianGrid,
    source_lats: np.ndarray,
    source_offsets: np.ndarray,
    target_lat: float,
    target_lons: np.ndarray,
) -> np.ndarray:
    south_row, north_row, south_weight, north_weight = _latitude_bracket(source_lats, target_lat)
    south = _interpolate_source_row(values, source_grid, source_offsets, south_row, target_lons)
    if north_row == south_row:
        return south
    north = _interpolate_source_row(values, source_grid, source_offsets, north_row, target_lons)
    return _weighted_nanmean(
        (south, north),
        (south_weight, north_weight),
    )


def _interpolate_source_row(
    values: np.ndarray,
    grid: GaussianGrid,
    offsets: np.ndarray,
    row_index: int,
    target_lons: np.ndarray,
) -> np.ndarray:
    start = int(offsets[row_index])
    nlon = int(_row_lengths(grid)[row_index])
    row = values[:, start : start + nlon]
    return _periodic_interp_missing(row, target_lons)


def _periodic_interp_missing(row: np.ndarray, target_lons: np.ndarray) -> np.ndarray:
    nlon = row.shape[1]
    pos = (target_lons % 360.0) * (nlon / 360.0)
    left_idx = np.floor(pos).astype(np.int64) % nlon
    right_idx = (left_idx + 1) % nlon
    right_weight = pos - np.floor(pos)
    left_weight = 1.0 - right_weight

    left = row[:, left_idx]
    right = row[:, right_idx]
    return _weighted_nanmean((left, right), (left_weight, right_weight))


def _weighted_nanmean(values: tuple[np.ndarray, ...], weights: tuple[np.ndarray | float, ...]) -> np.ndarray:
    out = np.zeros_like(values[0], dtype=np.float64)
    weight_sum = np.zeros_like(values[0], dtype=np.float64)
    for value, weight in zip(values, weights):
        valid = np.isfinite(value)
        weighted = np.where(valid, value, 0.0) * weight
        out += weighted
        weight_sum += valid.astype(np.float64) * weight
    return np.divide(out, weight_sum, out=np.full_like(out, np.nan), where=weight_sum > 0.0)


def _latitude_bracket(source_lats: np.ndarray, target_lat: float) -> tuple[int, int, float, float]:
    lat_asc = source_lats[::-1]
    if target_lat <= lat_asc[0]:
        row = source_lats.size - 1
        return row, row, 1.0, 0.0
    if target_lat >= lat_asc[-1]:
        return 0, 0, 1.0, 0.0

    upper_asc = int(np.searchsorted(lat_asc, target_lat, side="left"))
    lower_asc = upper_asc - 1
    lower_lat = float(lat_asc[lower_asc])
    upper_lat = float(lat_asc[upper_asc])
    north_weight = (target_lat - lower_lat) / (upper_lat - lower_lat)
    south_weight = 1.0 - north_weight
    south_row = source_lats.size - 1 - lower_asc
    north_row = source_lats.size - 1 - upper_asc
    return south_row, north_row, south_weight, north_weight


def _row_lengths(grid: GaussianGrid) -> np.ndarray:
    if grid.pl is None:
        return np.full(grid.nlat, grid.work_nlon, dtype=np.int64)
    return np.asarray(grid.pl, dtype=np.int64)


def _row_offsets(grid: GaussianGrid) -> np.ndarray:
    row_lengths = _row_lengths(grid)
    return np.concatenate(([0], np.cumsum(row_lengths[:-1]))).astype(np.int64)


def _chunk_slices(nfield: int, chunk_size: int | None):
    if chunk_size is None:
        yield 0, nfield
        return
    chunk = int(chunk_size)
    if chunk <= 0:
        raise ValueError("chunk_size must be a positive integer or None")
    for start in range(0, nfield, chunk):
        yield start, min(start + chunk, nfield)


def _horizontal_shape(grid: GaussianGrid) -> tuple[int, ...]:
    if grid.is_reduced:
        return (grid.size,)
    return (grid.nlat, grid.work_nlon)


def _move_horizontal_to_end(
    arr: np.ndarray,
    expected_shape: tuple[int, ...],
    axis: int | tuple[int, ...],
) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...] | None]:
    horizontal_ndim = len(expected_shape)
    if horizontal_ndim == 1:
        axes = (int(axis),)
    elif isinstance(axis, tuple):
        axes = tuple(int(v) for v in axis)
    elif int(axis) == -1:
        axes = tuple(range(arr.ndim - horizontal_ndim, arr.ndim))
    else:
        raise ValueError(
            "regular Gaussian input uses two horizontal dimensions; pass axis=(lat_axis, lon_axis)"
        )

    axes = tuple(a + arr.ndim if a < 0 else a for a in axes)
    if len(axes) != horizontal_ndim or len(set(axes)) != horizontal_ndim:
        raise ValueError(f"axis must identify {horizontal_ndim} horizontal dimension(s)")
    if any(a < 0 or a >= arr.ndim for a in axes):
        raise ValueError("axis is out of bounds")
    actual_shape = tuple(arr.shape[a] for a in axes)
    if actual_shape != expected_shape:
        raise ValueError(f"horizontal shape must be {expected_shape}, got {actual_shape}")
    leading_axes = tuple(i for i in range(arr.ndim) if i not in axes)
    moved = np.transpose(arr, leading_axes + axes)
    leading_shape = moved.shape[: arr.ndim - horizontal_ndim]
    if axes == tuple(range(arr.ndim - horizontal_ndim, arr.ndim)):
        return moved, leading_shape, None
    return moved, leading_shape, axes


def _move_horizontal_from_end(
    arr: np.ndarray,
    target_horizontal_shape: tuple[int, ...],
    original_axes: tuple[int, ...],
    source_horizontal_ndim: int,
) -> np.ndarray:
    target_horizontal_ndim = len(target_horizontal_shape)
    leading_count = arr.ndim - target_horizontal_ndim
    input_ndim = leading_count + source_horizontal_ndim
    source_axes = tuple(axis if axis >= 0 else axis + input_ndim for axis in original_axes)
    leading_source_axes = [axis for axis in range(input_ndim) if axis not in source_axes]
    leading_result_axes = list(range(leading_count))
    target_result_axes = list(range(leading_count, arr.ndim))
    result_axis_by_label = dict(zip(leading_source_axes, leading_result_axes))
    insert_at = min(source_axes)
    desired_labels = []
    inserted_target = False
    for axis in range(input_ndim):
        if axis == insert_at:
            desired_labels.extend(("target", idx) for idx in range(target_horizontal_ndim))
            inserted_target = True
        if axis not in source_axes:
            desired_labels.append(("leading", axis))
    if not inserted_target:
        desired_labels.extend(("target", idx) for idx in range(target_horizontal_ndim))

    output_order = []
    for kind, value in desired_labels:
        if kind == "target":
            output_order.append(target_result_axes[value])
        else:
            output_order.append(result_axis_by_label[value])
    return np.transpose(arr, output_order)


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
    rows = [regular_longitudes(int(row_nlon)) for row_nlon in grid.pl]
    return np.concatenate(rows)


def _output_attrs(
    attrs,
    source_grid: GaussianGrid,
    target_grid: GaussianGrid,
    *,
    keep_attrs: bool,
    method: str,
    chunk_size: int | None,
) -> dict:
    out = dict(attrs) if keep_attrs else {}
    out = preserve_source_grid_attrs(out, source_grid)
    out["GRIB_N"] = target_grid.n
    out["GRIB_gridType"] = "reduced_gg" if target_grid.is_reduced else "regular_gg"
    out["GRIB_numberOfPoints"] = target_grid.size
    out["fullpos_regrid_method"] = method
    if target_grid.is_reduced:
        out["GRIB_pl"] = np.asarray(target_grid.pl, dtype=np.int64)
    else:
        out.pop("GRIB_pl", None)
    out = append_history(
        out,
        format_regrid_history(
            source_grid=source_grid,
            target_grid=target_grid,
            method=method,
            chunk_size=chunk_size,
        ),
    )
    return out
