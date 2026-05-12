from __future__ import annotations

import numpy as np
import xarray as xr

from ..errors import FullposNotImplementedError
from ..grids import GaussianGrid, gaussian_latitudes, parse_grid
from .kernels import horizontal_halo_kernel, horizontal_regular_kernel


def horizontal_interpolate(
    values,
    *,
    source_grid: str | GaussianGrid | None = None,
    source_pl=None,
    target_lats,
    target_lons,
    method: str = "bilinear",
    axis: int | tuple[int, int] | None = None,
    missing_value=None,
    source_mask=None,
    average_radius: int = 1,
    chunks: dict[str, int] | None = None,
    keep_attrs: bool = True,
    variables=None,
    skip_non_horizontal: bool = True,
) -> np.ndarray:
    """Interpolate horizontal fields with the native FULLPOS interpolation path.

    The public API is intentionally reserved for OpenIFS/FULLPOS routines
    (FPINT4, FPINT12, FPNEAR, FPAVG). It accepts the planned xarray-style
    ``chunks`` argument, for example ``{"time": 1, "hybrid": 10}``, but no
    Python fallback is exposed here.
    """
    _validate_horizontal_request(
        values,
        source_grid=source_grid,
        source_pl=source_pl,
        target_lats=target_lats,
        target_lons=target_lons,
        method=method,
        axis=axis,
        chunks=chunks,
        variables=variables,
        skip_non_horizontal=skip_non_horizontal,
    )
    _ = missing_value
    if source_mask is not None:
        raise FullposNotImplementedError(
            "native FULLPOS mask-aware horizontal interpolation is not implemented yet"
        )
    if isinstance(values, xr.Dataset):
        return _horizontal_interpolate_dataset(
            values,
            target_lats=target_lats,
            target_lons=target_lons,
            method=method,
            source_grid=source_grid,
            source_pl=source_pl,
            chunks=chunks,
            keep_attrs=keep_attrs,
            average_radius=average_radius,
            variables=variables,
            skip_non_horizontal=skip_non_horizontal,
        )
    if isinstance(values, xr.DataArray):
        return _horizontal_interpolate_data_array(
            values,
            target_lats=target_lats,
            target_lons=target_lons,
            method=method,
            source_grid=source_grid,
            source_pl=source_pl,
            chunks=chunks,
            keep_attrs=keep_attrs,
            average_radius=average_radius,
        )
    return _horizontal_interpolate_numpy(
        values,
        source_grid=source_grid,
        source_pl=source_pl,
        target_lats=target_lats,
        target_lons=target_lons,
        method=method,
        axis=axis,
        average_radius=average_radius,
    )


def nearest_interpolate(values, **kwargs) -> np.ndarray:
    """Interpolate using nearest-neighbour selection."""
    return horizontal_interpolate(values, method="nearest", **kwargs)


def bilinear_interpolate(values, **kwargs) -> np.ndarray:
    """Interpolate using 4-point bilinear weights."""
    return horizontal_interpolate(values, method="bilinear", **kwargs)


def quadratic12_interpolate(values, **kwargs) -> np.ndarray:
    """Interpolate using a 12-point high-order regular-grid stencil."""
    return horizontal_interpolate(values, method="quadratic12", **kwargs)


def average_interpolate(values, **kwargs) -> np.ndarray:
    """Interpolate using an average over a square source-grid halo."""
    return horizontal_interpolate(values, method="average", **kwargs)


def _horizontal_interpolate_numpy(
    values,
    *,
    source_grid,
    source_pl,
    target_lats,
    target_lons,
    method: str,
    axis: int | tuple[int, int],
    average_radius: int,
) -> np.ndarray:
    arr = np.asarray(values)
    nloen = _resolve_source_pl(source_pl=source_pl, source_grid=source_grid, attrs=None)
    source_shape = _resolve_numpy_source_shape(source_pl=nloen, values=arr)
    tgt_lats = np.asarray(target_lats, dtype=np.float64)
    tgt_lons = np.asarray(target_lons, dtype=np.float64)
    moved, leading_shape, original_axes = _move_source_to_end(
        arr,
        source_shape,
        axis,
    )
    flat = moved.reshape((-1, int(np.prod(source_shape))))
    out = np.empty((flat.shape[0],) + tgt_lats.shape, dtype=np.float64)
    kernel = _select_horizontal_kernel(method)
    for idx, field in enumerate(flat):
        out[idx] = kernel(
            field if nloen is not None else field.reshape(source_shape),
            source_pl=nloen,
            target_lats=tgt_lats,
            target_lons=tgt_lons,
            method=_normalize_method(method),
            **_halo_kwargs(method, average_radius),
        )
    result = out.reshape(leading_shape + tgt_lats.shape)
    if original_axes is None:
        return result
    return _move_horizontal_from_end(result, tgt_lats.shape, original_axes, len(source_shape))


def _horizontal_interpolate_data_array(
    obj: xr.DataArray,
    *,
    target_lats,
    target_lons,
    method: str,
    source_grid: str | GaussianGrid | None,
    source_pl,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    average_radius: int = 1,
) -> xr.DataArray:
    source = _resolve_xarray_source(
        obj,
        source_grid=source_grid,
        source_pl=source_pl,
    )
    leading_dims = [dim for dim in obj.dims if dim not in source["dims"]]
    transposed = obj.transpose(*leading_dims, *source["dims"])
    src_lats = source["lats"]
    src_lons = source["lons"]
    nloen = source["pl"]
    tgt_lats = np.asarray(target_lats, dtype=np.float64)
    tgt_lons = np.asarray(target_lons, dtype=np.float64)
    leading_chunks = _resolve_named_chunks(transposed, leading_dims, chunks)
    out_shape = transposed.shape[: len(leading_dims)] + tgt_lats.shape
    out_values = np.empty(out_shape, dtype=np.float64)
    if not leading_dims:
        kernel = _select_horizontal_kernel(method)
        out_values[...] = kernel(
            _native_source_field(transposed.values, nloen),
            source_lats=src_lats,
            source_lons=src_lons,
            source_pl=nloen,
            target_lats=tgt_lats,
            target_lons=tgt_lons,
            method=_normalize_method(method),
            **_halo_kwargs(method, average_radius),
        )
    else:
        kernel = _select_horizontal_kernel(method)
        chunk_sizes = leading_chunks or tuple(transposed.sizes[dim] for dim in leading_dims)
        for leading_slice in _iter_leading_slices(transposed.shape[: len(leading_dims)], chunk_sizes):
            if nloen is not None:
                block = transposed.values[leading_slice + (slice(None),)]
                flat = block.reshape((-1, int(nloen.sum())))
            else:
                block = transposed.values[leading_slice + (slice(None), slice(None))]
                flat = block.reshape((-1, src_lats.size, src_lons.size))
            block_out = np.empty((flat.shape[0],) + tgt_lats.shape, dtype=np.float64)
            for idx, field in enumerate(flat):
                block_out[idx] = kernel(
                    field,
                    source_lats=src_lats,
                    source_lons=src_lons,
                    source_pl=nloen,
                    target_lats=tgt_lats,
                    target_lons=tgt_lons,
                    method=_normalize_method(method),
                    **_halo_kwargs(method, average_radius),
                )
            out_values[leading_slice + tuple(slice(None) for _ in tgt_lats.shape)] = block_out.reshape(
                block.shape[: -len(source["dims"])] + tgt_lats.shape
            )

    target_dims, coords = _target_dims_and_coords(obj, leading_dims, tgt_lats, tgt_lons)
    attrs = dict(obj.attrs) if keep_attrs else {}
    return xr.DataArray(out_values, dims=target_dims, coords=coords, name=obj.name, attrs=attrs)


def _horizontal_interpolate_dataset(
    obj: xr.Dataset,
    *,
    target_lats,
    target_lons,
    method: str,
    source_grid: str | GaussianGrid | None,
    source_pl,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    average_radius: int,
    variables,
    skip_non_horizontal: bool,
) -> xr.Dataset:
    selected = list(obj.data_vars) if variables is None else [str(v) for v in variables]
    out = xr.Dataset(attrs=dict(obj.attrs) if keep_attrs else {})
    for name, data_array in obj.data_vars.items():
        if name not in selected:
            out[name] = data_array
            continue
        try:
            _find_xarray_source_dims(data_array, source_grid=source_grid, source_pl=source_pl)
        except ValueError:
            if variables is None and skip_non_horizontal:
                out[name] = data_array
                continue
            raise ValueError(f"variable {name!r} does not have supported horizontal dimensions")
        out[name] = _horizontal_interpolate_data_array(
            data_array,
            target_lats=target_lats,
            target_lons=target_lons,
            method=method,
            source_grid=source_grid,
            source_pl=source_pl,
            chunks=chunks,
            keep_attrs=keep_attrs,
            average_radius=average_radius,
        )
    return out


def _normalize_method(method: str) -> str:
    normalized = str(method).lower().replace("-", "_")
    aliases = {
        "4": "bilinear",
        "fpint4": "bilinear",
        "linear": "bilinear",
        "12": "quadratic12",
        "fpint12": "quadratic12",
        "quadratic": "quadratic12",
        "nearest_neighbour": "nearest",
        "nearest_neighbor": "nearest",
        "fpnear": "nearest",
        "avg": "average",
        "fpavg": "average",
        "masked_average": "average",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"nearest", "bilinear", "quadratic12", "average"}:
        raise ValueError("method must be 'nearest', 'bilinear', 'quadratic12', or 'average'")
    return normalized


def _select_horizontal_kernel(method: str):
    normalized = _normalize_method(method)
    if normalized in {"nearest", "average"}:
        return horizontal_halo_kernel
    return horizontal_regular_kernel


def _halo_kwargs(method: str, average_radius: int) -> dict:
    if _normalize_method(method) not in {"nearest", "average"}:
        return {}
    radius = int(average_radius)
    if radius <= 0:
        raise ValueError("average_radius must be positive for nearest/average interpolation")
    return {"kslwide": radius}


def _validate_horizontal_request(
    values,
    *,
    source_grid,
    source_pl,
    target_lats,
    target_lons,
    method: str,
    axis: int | tuple[int, int],
    chunks: dict[str, int] | None,
    variables,
    skip_non_horizontal: bool,
) -> None:
    _normalize_method(method)
    target_lats_arr = np.asarray(target_lats)
    target_lons_arr = np.asarray(target_lons)
    if target_lats_arr.shape != target_lons_arr.shape:
        raise ValueError("target_lats and target_lons must have matching shapes")
    if isinstance(values, xr.Dataset):
        _validate_dataset_request(
            values,
            source_grid=source_grid,
            source_pl=source_pl,
            chunks=chunks,
            variables=variables,
            skip_non_horizontal=skip_non_horizontal,
        )
        return
    if isinstance(values, xr.DataArray):
        _validate_data_array_request(values, source_grid=source_grid, source_pl=source_pl, chunks=chunks)
        return
    if chunks is not None:
        raise TypeError("chunks is only supported for xarray inputs with named dimensions")
    resolved_pl = _resolve_source_pl(source_pl=source_pl, source_grid=source_grid, attrs=None)
    if resolved_pl is not None:
        arr = np.asarray(values)
        axes = _normalize_axis(axis, arr.ndim, 1)
        if tuple(arr.shape[a] for a in axes) != (int(resolved_pl.sum()),):
            raise ValueError(
                f"packed horizontal shape must be {(int(resolved_pl.sum()),)}, "
                f"got {tuple(arr.shape[a] for a in axes)}"
            )
        return
    raise ValueError(
        "NumPy inputs are only supported for packed reduced Gaussian fields; "
        "use xarray DataArray/Dataset for regular-row horizontal interpolation"
    )


def _validate_data_array_request(
    obj: xr.DataArray,
    *,
    source_grid: str | GaussianGrid | None,
    source_pl,
    chunks: dict[str, int] | None,
) -> None:
    source_dims = _find_xarray_source_dims(obj, source_grid=source_grid, source_pl=source_pl)
    _resolve_named_chunks(obj, [dim for dim in obj.dims if dim not in source_dims], chunks)


def _validate_dataset_request(
    obj: xr.Dataset,
    *,
    source_grid: str | GaussianGrid | None,
    source_pl,
    chunks: dict[str, int] | None,
    variables,
    skip_non_horizontal: bool,
) -> None:
    selected = list(obj.data_vars) if variables is None else [str(v) for v in variables]
    missing = [name for name in selected if name not in obj.data_vars]
    if missing:
        raise KeyError(f"variables not found in dataset: {missing}")
    for name in selected:
        data_array = obj[name]
        try:
            _find_xarray_source_dims(data_array, source_grid=source_grid, source_pl=source_pl)
        except ValueError:
            if variables is None and skip_non_horizontal:
                continue
            raise ValueError(f"variable {name!r} does not have supported horizontal dimensions")
        _validate_data_array_request(
            data_array,
            source_grid=source_grid,
            source_pl=source_pl,
            chunks=chunks,
        )


def _resolve_named_chunks(
    obj: xr.DataArray,
    leading_dims: list[str],
    chunks: dict[str, int] | None,
) -> tuple[int, ...] | None:
    if chunks is None:
        return None
    unknown = set(chunks) - set(obj.dims)
    if unknown:
        raise ValueError(f"chunks contains dimensions not present in object: {sorted(unknown)}")
    resolved = []
    for dim in leading_dims:
        size = int(chunks.get(dim, obj.sizes[dim]))
        if size <= 0:
            raise ValueError("chunk sizes must be positive")
        resolved.append(size)
    return tuple(resolved)


def _iter_leading_slices(
    leading_shape: tuple[int, ...],
    leading_chunks: tuple[int, ...],
):
    ranges = [
        [slice(start, min(length, start + chunk)) for start in range(0, length, chunk)]
        for length, chunk in zip(leading_shape, leading_chunks)
    ]
    if not ranges:
        yield ()
        return
    yield from _slice_product(ranges, 0, [])


def _slice_product(ranges: list[list[slice]], idx: int, current: list[slice]):
    if idx == len(ranges):
        yield tuple(current)
        return
    for item in ranges[idx]:
        current.append(item)
        yield from _slice_product(ranges, idx + 1, current)
        current.pop()


def _target_dims_and_coords(
    obj: xr.DataArray,
    leading_dims: list[str],
    target_lats: np.ndarray,
    target_lons: np.ndarray,
) -> tuple[tuple[str, ...], dict]:
    coords = {}
    leading_set = set(leading_dims)
    for name, coord in obj.coords.items():
        if set(coord.dims).issubset(leading_set):
            coords[name] = coord
    if target_lats.ndim == 1:
        dims = tuple(leading_dims) + ("points",)
        coords["latitude"] = ("points", target_lats)
        coords["longitude"] = ("points", target_lons)
        return dims, coords
    if target_lats.ndim == 2:
        dims = tuple(leading_dims) + ("target_y", "target_x")
        coords["latitude"] = (("target_y", "target_x"), target_lats)
        coords["longitude"] = (("target_y", "target_x"), target_lons)
        return dims, coords
    target_dims = tuple(f"target_{idx}" for idx in range(target_lats.ndim))
    dims = tuple(leading_dims) + target_dims
    coords["latitude"] = (target_dims, target_lats)
    coords["longitude"] = (target_dims, target_lons)
    return dims, coords


def _resolve_numpy_source_shape(
    *,
    source_pl: np.ndarray | None,
    values: np.ndarray,
) -> tuple[int, ...]:
    if source_pl is not None:
        return (int(source_pl.sum()),)
    raise ValueError(
        "NumPy inputs are only supported for packed reduced Gaussian fields; "
        f"got regular array shape {values.shape}"
    )


def _resolve_xarray_source(
    obj: xr.DataArray,
    *,
    source_grid: str | GaussianGrid | None,
    source_pl,
) -> dict:
    dims = _find_xarray_source_dims(obj, source_grid=source_grid, source_pl=source_pl)
    pl = _resolve_source_pl(source_pl=source_pl, source_grid=source_grid, attrs=obj.attrs)
    if len(dims) == 1:
        assert pl is not None
        return {
            "dims": dims,
            "pl": pl,
            "lats": gaussian_latitudes(pl.size),
            "lons": None,
        }
    return {
        "dims": dims,
        "pl": None,
        "lats": obj["latitude"].values,
        "lons": obj["longitude"].values,
    }


def _find_xarray_source_dims(
    obj: xr.DataArray,
    *,
    source_grid: str | GaussianGrid | None,
    source_pl,
) -> tuple[str, ...]:
    pl = _resolve_source_pl(source_pl=source_pl, source_grid=source_grid, attrs=obj.attrs)
    if pl is not None:
        packed_size = int(pl.sum())
        if "values" in obj.dims and obj.sizes["values"] == packed_size:
            return ("values",)
        matches = [dim for dim in obj.dims if obj.sizes[dim] == packed_size]
        if len(matches) == 1:
            return (matches[0],)
        raise ValueError(
            "xarray packed reduced interpolation requires one dimension "
            f"with size {packed_size}"
        )
    if "latitude" in obj.dims and "longitude" in obj.dims:
        return ("latitude", "longitude")
    raise ValueError("xarray horizontal interpolation requires latitude/longitude or packed values")


def _resolve_source_pl(
    *,
    source_pl,
    source_grid: str | GaussianGrid | None,
    attrs,
) -> np.ndarray | None:
    if source_pl is not None:
        return _validate_source_pl(source_pl)
    if source_grid is not None:
        grid = source_grid if isinstance(source_grid, GaussianGrid) else parse_grid(source_grid)
        if grid.pl is not None:
            return _validate_source_pl(grid.pl)
        return None
    if attrs:
        for key in ("GRIB_pl", "pl"):
            if key in attrs:
                return _validate_source_pl(attrs[key])
    return None


def _validate_source_pl(source_pl) -> np.ndarray:
    pl = np.asarray(source_pl, dtype=np.int32).reshape(-1)
    if pl.size == 0 or np.any(pl <= 0):
        raise ValueError("source_pl must contain positive row lengths")
    return pl


def _native_source_field(values: np.ndarray, source_pl: np.ndarray | None) -> np.ndarray:
    arr = np.asarray(values)
    if source_pl is None:
        return arr
    return arr.reshape(int(source_pl.sum()))


def _normalize_axis(axis, ndim: int, horizontal_ndim: int) -> tuple[int, ...]:
    if axis is None:
        return tuple(range(ndim - horizontal_ndim, ndim))
    if horizontal_ndim == 1:
        if isinstance(axis, tuple):
            axes = tuple(int(v) for v in axis)
        else:
            axes = (int(axis),)
    else:
        if not isinstance(axis, tuple):
            if int(axis) != -1:
                raise ValueError("regular horizontal interpolation requires axis=(lat_axis, lon_axis)")
            axes = tuple(range(ndim - horizontal_ndim, ndim))
        else:
            axes = tuple(int(v) for v in axis)
    axes = tuple(a + ndim if a < 0 else a for a in axes)
    if len(axes) != horizontal_ndim or len(set(axes)) != horizontal_ndim:
        raise ValueError(f"axis must identify {horizontal_ndim} horizontal dimension(s)")
    if any(a < 0 or a >= ndim for a in axes):
        raise ValueError("axis is out of bounds")
    return axes


def _move_source_to_end(
    arr: np.ndarray,
    expected_shape: tuple[int, ...],
    axis,
) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...] | None]:
    axes = _normalize_axis(axis, arr.ndim, len(expected_shape))
    actual_shape = tuple(arr.shape[a] for a in axes)
    if actual_shape != expected_shape:
        raise ValueError(f"horizontal shape must be {expected_shape}, got {actual_shape}")
    leading_axes = tuple(i for i in range(arr.ndim) if i not in axes)
    moved = np.transpose(arr, leading_axes + axes)
    leading_shape = moved.shape[: arr.ndim - len(expected_shape)]
    if axes == tuple(range(arr.ndim - len(expected_shape), arr.ndim)):
        return moved, leading_shape, None
    return moved, leading_shape, axes


def _move_horizontal_from_end(
    arr: np.ndarray,
    target_horizontal_shape: tuple[int, ...],
    original_axes: tuple[int, int],
    source_horizontal_ndim: int,
) -> np.ndarray:
    leading_count = arr.ndim - len(target_horizontal_shape)
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
            desired_labels.extend(("target", idx) for idx in range(len(target_horizontal_shape)))
            inserted_target = True
        if axis not in source_axes:
            desired_labels.append(("leading", axis))
    if not inserted_target:
        desired_labels.extend(("target", idx) for idx in range(len(target_horizontal_shape)))
    output_order = []
    for kind, value in desired_labels:
        output_order.append(target_result_axes[value] if kind == "target" else result_axis_by_label[value])
    return np.transpose(arr, output_order)
