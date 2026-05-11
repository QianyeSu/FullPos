from __future__ import annotations

import numpy as np
import xarray as xr

from ..grids import GaussianGrid, gaussian_latitudes, parse_grid, regular_longitudes
from ..metadata import infer_source_grid
from .native import spectral_regrid_chunks
from .surface import masked_regrid_data_array, masked_regrid_values
from .xarray import regrid_data_array


DEFAULT_CHUNK_SIZE = 64
_USE_INSTANCE_CHUNK_SIZE = object()


class Regridder:
    """Reusable interpolation object for one source/target Gaussian grid pair.

    The object stores parsed grid metadata and backend options so repeated
    calls over many fields avoid rebuilding user-facing configuration.
    """

    def __init__(
        self,
        source_grid: str | GaussianGrid,
        target_grid: str | GaussianGrid,
        *,
        ntrunc: int | None = None,
        chunk_size: int | None = DEFAULT_CHUNK_SIZE,
        missing_policy: str = "error",
        method: str = "linear",
        missing_value=None,
    ) -> None:
        """Create a regridder for a fixed source and target grid."""
        self.source_grid = _ensure_grid(source_grid)
        self.target_grid = _ensure_grid(target_grid)
        self.ntrunc = ntrunc
        self.chunk_size = chunk_size
        self.missing_policy = missing_policy
        self.method = _normalize_method(method)
        self.missing_value = missing_value

    @classmethod
    def from_dataarray(
        cls,
        obj: xr.DataArray,
        *,
        target_grid: str | GaussianGrid,
        source_grid: str | GaussianGrid | None = None,
        ntrunc: int | None = None,
        chunk_size: int | None = DEFAULT_CHUNK_SIZE,
        missing_policy: str = "error",
        method: str = "linear",
        missing_value=None,
    ) -> "Regridder":
        """Build a regridder by inferring the source grid from a DataArray."""
        return cls(
            source_grid or infer_source_grid(obj),
            target_grid,
            ntrunc=ntrunc,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
            method=method,
            missing_value=missing_value,
        )

    @classmethod
    def from_dataset(
        cls,
        obj: xr.Dataset,
        *,
        target_grid: str | GaussianGrid,
        source_grid: str | GaussianGrid | None = None,
        ntrunc: int | None = None,
        chunk_size: int | None = DEFAULT_CHUNK_SIZE,
        missing_policy: str = "error",
        method: str = "linear",
        missing_value=None,
    ) -> "Regridder":
        """Build a regridder by inferring the source grid from a Dataset."""
        return cls(
            source_grid or infer_source_grid(obj),
            target_grid,
            ntrunc=ntrunc,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
            method=method,
            missing_value=missing_value,
        )

    @property
    def source_shape(self) -> tuple[int, ...]:
        """Horizontal input shape expected by :meth:`regrid_values`."""
        return _horizontal_shape(self.source_grid)

    @property
    def target_shape(self) -> tuple[int, ...]:
        """Horizontal output shape produced by this regridder."""
        return _horizontal_shape(self.target_grid)

    @property
    def source_dims(self) -> tuple[str, ...]:
        """Canonical horizontal dimension names for the source grid."""
        return _horizontal_dims(self.source_grid)

    @property
    def target_dims(self) -> tuple[str, ...]:
        """Canonical horizontal dimension names for the target grid."""
        return _horizontal_dims(self.target_grid)

    @property
    def target_coords(self) -> dict:
        """Coordinate arrays suitable for constructing target xarray objects."""
        if self.target_grid.is_reduced:
            return {
                "values": np.arange(self.target_grid.size),
                "latitude": ("values", _packed_latitudes(self.target_grid)),
                "longitude": ("values", _packed_longitudes(self.target_grid)),
            }
        return {
            "latitude": gaussian_latitudes(self.target_grid.nlat),
            "longitude": regular_longitudes(self.target_grid.work_nlon),
        }

    @property
    def info(self) -> dict:
        """Return a serializable summary of grids and backend options."""
        return {
            "source_grid": self.source_grid.name,
            "target_grid": self.target_grid.name,
            "method": self.method,
            "source_shape": self.source_shape,
            "target_shape": self.target_shape,
            "source_size": self.source_grid.size,
            "target_size": self.target_grid.size,
            "ntrunc": self.ntrunc,
            "chunk_size": self.chunk_size,
            "missing_policy": self.missing_policy,
        }

    def regrid_values(
        self,
        values,
        *,
        axis: int | tuple[int, ...] = -1,
        chunk_size: int | None | object = _USE_INSTANCE_CHUNK_SIZE,
    ) -> np.ndarray:
        """Regrid NumPy-like values using this object's configured grid pair."""
        resolved_chunk_size = self._resolve_chunk_size(chunk_size)
        if self.method == "masked":
            return masked_regrid_values(
                values,
                source_grid=self.source_grid,
                target_grid=self.target_grid,
                axis=axis,
                chunk_size=resolved_chunk_size,
                missing_value=self.missing_value,
            )

        arr = np.asarray(values)
        expected_shape = _horizontal_shape(self.source_grid)

        horizontal_ndim = len(expected_shape)
        arr, leading_shape, moved_back = _move_horizontal_to_end(arr, expected_shape, axis)
        flat = arr.reshape((-1, self.source_grid.size))
        out = spectral_regrid_chunks(
            flat,
            source_grid=self.source_grid,
            target_grid=self.target_grid,
            ntrunc=self.ntrunc,
            chunk_size=resolved_chunk_size,
            missing_policy=self.missing_policy,
        )

        target_shape = _horizontal_shape(self.target_grid)
        result = out.reshape(leading_shape + target_shape)
        if moved_back is None:
            return result
        return _move_horizontal_from_end(result, target_shape, moved_back, horizontal_ndim)

    def regrid_data_array(
        self,
        obj: xr.DataArray,
        *,
        keep_attrs: bool = True,
        chunk_size: int | None | object = _USE_INSTANCE_CHUNK_SIZE,
    ) -> xr.DataArray:
        """Regrid an ``xarray.DataArray`` using this object's grid pair."""
        resolved_chunk_size = self._resolve_chunk_size(chunk_size)
        if self.method == "masked":
            return masked_regrid_data_array(
                obj,
                source_grid=self.source_grid,
                target_grid=self.target_grid,
                chunk_size=resolved_chunk_size,
                missing_value=self.missing_value,
                keep_attrs=keep_attrs,
            )

        return regrid_data_array(
            obj,
            source_grid=self.source_grid,
            target_grid=self.target_grid,
            ntrunc=self.ntrunc,
            chunk_size=resolved_chunk_size,
            missing_policy=self.missing_policy,
            keep_attrs=keep_attrs,
        )

    def _resolve_chunk_size(self, chunk_size: int | None | object) -> int | None:
        if chunk_size is _USE_INSTANCE_CHUNK_SIZE:
            return self.chunk_size
        return chunk_size  # type: ignore[return-value]


def _ensure_grid(grid: str | GaussianGrid) -> GaussianGrid:
    if isinstance(grid, GaussianGrid):
        return grid
    return parse_grid(grid)


def _normalize_method(method: str) -> str:
    normalized = str(method).lower()
    if normalized == "spectral":
        return "linear"
    if normalized not in {"linear", "masked"}:
        raise ValueError("method must be 'linear', 'spectral', or 'masked'")
    return normalized


def _horizontal_shape(grid: GaussianGrid) -> tuple[int, ...]:
    if grid.is_reduced:
        return (grid.size,)
    return (grid.nlat, grid.work_nlon)


def _horizontal_dims(grid: GaussianGrid) -> tuple[str, ...]:
    if grid.is_reduced:
        return ("values",)
    return ("latitude", "longitude")


def _packed_latitudes(grid: GaussianGrid) -> np.ndarray:
    assert grid.pl is not None
    lats = gaussian_latitudes(grid.nlat)
    return np.repeat(lats, np.asarray(grid.pl, dtype=np.int64))


def _packed_longitudes(grid: GaussianGrid) -> np.ndarray:
    assert grid.pl is not None
    rows = [regular_longitudes(int(row_nlon)) for row_nlon in grid.pl]
    return np.concatenate(rows)


def _move_horizontal_to_end(
    arr: np.ndarray,
    expected_shape: tuple[int, ...],
    axis: int | tuple[int, ...],
) -> tuple[np.ndarray, tuple[int, ...], tuple[int, ...] | None]:
    horizontal_ndim = len(expected_shape)
    if horizontal_ndim == 1:
        axes = (int(axis),)
    else:
        if not isinstance(axis, tuple):
            if int(axis) != -1:
                raise ValueError(
                    "regular Gaussian input uses two horizontal dimensions; "
                    "pass axis=(lat_axis, lon_axis)"
                )
            axes = tuple(range(arr.ndim - horizontal_ndim, arr.ndim))
        else:
            axes = tuple(int(v) for v in axis)
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
