from __future__ import annotations

import numpy as np
import xarray as xr

from ..grids import GaussianGrid, gaussian_latitudes, parse_grid
from ..metadata import append_history, format_regrid_history, infer_source_grid
from .horizontal import horizontal_interpolate
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
    chunks: dict[str, int] | None = None,
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
        Grid names such as ``"O320"``, ``"O480"``, ``"F320"``, or parsed grid
        objects. ``N*`` regular Gaussian names are accepted as compatibility
        aliases, but ``F*`` is the preferred public spelling.
    method:
        ``"linear"``/``"spectral"`` uses the native ECTRANS spectral path.
        ``"masked"`` uses the NaN-aware surface helper for fields with bitmaps.
        For regular latitude/longitude targets, ``"linear"`` maps to native
        FULLPOS bilinear horizontal interpolation, and the native horizontal
        methods ``"nearest"``, ``"bilinear"``, ``"quadratic12"``, and
        ``"average"`` are also accepted.
    variables:
        Optional dataset variable subset. Non-horizontal variables are skipped
        by default when ``variables`` is not explicit.
    ntrunc:
        Optional triangular truncation used by the spectral backend.
    chunk_size:
        Number of flattened leading fields transformed per native backend call.
    chunks:
        Optional named leading-dimension chunk sizes for native FULLPOS
        horizontal interpolation, for example ``{"time": 1, "hybrid": 10}``.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        Regridded object with updated Gaussian-grid metadata and appended
        ``history``.
    """
    if target_grid is None:
        raise ValueError("target_grid is required")
    regular_ll = _parse_regular_latlon_target(target_grid)
    if regular_ll is not None:
        return _regrid_regular_latlon(
            obj,
            source_grid=source_grid,
            target=regular_ll,
            method=method,
            variables=variables,
            ntrunc=ntrunc,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
            chunks=chunks,
            keep_attrs=keep_attrs,
            skip_non_horizontal=skip_non_horizontal,
        )
    if method not in {"linear", "spectral", "masked"}:
        raise ValueError("method must be 'linear', 'spectral', or 'masked'")

    if isinstance(obj, xr.DataArray):
        regridder = Regridder(
            source_grid or infer_source_grid(obj),
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
    dimension for O-grids, or ``(lat_axis, lon_axis)`` for regular F-grids.
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


def spectral_interpolate(
    obj,
    *,
    source_grid=None,
    target_grid=None,
    variables=None,
    missing_value=None,
    ntrunc=None,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    keep_attrs=True,
    skip_non_horizontal=True,
):
    """Interpolate xarray data between Gaussian grids with native ECTRANS.

    This is the explicit spectral-interpolation spelling for the Gaussian
    ``O*``/``F*``/``N*`` path. It is equivalent to ``regrid(...,
    method="spectral")`` for Gaussian targets, but rejects regular
    latitude/longitude targets such as ``"LL1.0"`` because those are handled by
    native FULLPOS horizontal interpolation through ``regrid``.
    """
    if target_grid is None:
        raise ValueError("target_grid is required")
    if _parse_regular_latlon_target(target_grid) is not None:
        raise ValueError(
            "spectral_interpolate only supports Gaussian O*/F*/N* targets; "
            "use regrid(...) for regular latitude/longitude LL targets"
        )
    return regrid(
        obj,
        source_grid=source_grid,
        target_grid=target_grid,
        method="spectral",
        variables=variables,
        missing_value=missing_value,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        keep_attrs=keep_attrs,
        skip_non_horizontal=skip_non_horizontal,
    )


def spectral_regrid(
    obj,
    *,
    source_grid=None,
    target_grid=None,
    variables=None,
    missing_value=None,
    ntrunc=None,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    keep_attrs=True,
    skip_non_horizontal=True,
):
    """Regrid xarray data between Gaussian grids with native ECTRANS/FIAT.

    This is a semantic alias of ``spectral_interpolate``. It is kept separate
    from ``regrid`` so users can explicitly request the Gaussian spectral path.
    """
    return spectral_interpolate(
        obj,
        source_grid=source_grid,
        target_grid=target_grid,
        variables=variables,
        missing_value=missing_value,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        keep_attrs=keep_attrs,
        skip_non_horizontal=skip_non_horizontal,
    )


def spectral_interpolate_values(
    values,
    *,
    source_grid,
    target_grid,
    ntrunc=None,
    axis=-1,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    missing_value=None,
):
    """Interpolate NumPy-like values between Gaussian grids with ECTRANS."""
    if _parse_regular_latlon_target(target_grid) is not None:
        raise ValueError(
            "spectral_interpolate_values only supports Gaussian O*/F*/N* targets"
        )
    return regrid_values(
        values,
        source_grid=source_grid,
        target_grid=target_grid,
        ntrunc=ntrunc,
        axis=axis,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        method="spectral",
        missing_value=missing_value,
    )


def spectral_regrid_values(
    values,
    *,
    source_grid,
    target_grid,
    ntrunc=None,
    axis=-1,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    missing_value=None,
):
    """Regrid NumPy-like values between Gaussian grids with ECTRANS/FIAT."""
    return spectral_interpolate_values(
        values,
        source_grid=source_grid,
        target_grid=target_grid,
        ntrunc=ntrunc,
        axis=axis,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        missing_value=missing_value,
    )


def _regrid_regular_latlon(
    obj,
    *,
    source_grid,
    target: dict,
    method: str,
    variables,
    ntrunc,
    chunk_size,
    missing_policy: str,
    chunks: dict[str, int] | None,
    keep_attrs: bool,
    skip_non_horizontal: bool,
):
    if method == "spectral":
        raise ValueError(
            "regular latitude/longitude targets use native FULLPOS horizontal "
            "methods, not method='spectral'"
        )
    if method == "masked":
        raise ValueError(
            "regular latitude/longitude targets do not support method='masked' yet"
        )
    if not isinstance(obj, (xr.DataArray, xr.Dataset)):
        raise TypeError(
            "regular latitude/longitude regridding requires xarray DataArray "
            "or Dataset input"
        )

    latitudes = target["latitude"]
    longitudes = target["longitude"]
    target_lats, target_lons = np.meshgrid(latitudes, longitudes, indexing="ij")
    horizontal_method = "bilinear" if str(method).lower() == "linear" else method
    resolved_source_grid = source_grid
    if resolved_source_grid is None:
        try:
            resolved_source_grid = infer_source_grid(obj)
        except ValueError:
            resolved_source_grid = None
    obj, resolved_source_grid = _promote_source_if_regular_ll_needs_polar_coverage(
        obj,
        source_grid=resolved_source_grid,
        target_latitudes=latitudes,
        variables=variables,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        keep_attrs=keep_attrs,
        skip_non_horizontal=skip_non_horizontal,
    )

    out = horizontal_interpolate(
        obj,
        source_grid=resolved_source_grid,
        target_lats=target_lats,
        target_lons=target_lons,
        method=horizontal_method,
        variables=variables,
        chunks=chunks,
        keep_attrs=keep_attrs,
        skip_non_horizontal=skip_non_horizontal,
    )
    return _finalize_regular_latlon_output(
        out,
        latitudes=latitudes,
        longitudes=longitudes,
        source_grid=resolved_source_grid or source_grid or "inferred",
        target_name=target["name"],
        method=horizontal_method,
        keep_attrs=keep_attrs,
    )


def _promote_source_if_regular_ll_needs_polar_coverage(
    obj,
    *,
    source_grid,
    target_latitudes: np.ndarray,
    variables,
    ntrunc,
    chunk_size,
    missing_policy: str,
    keep_attrs: bool,
    skip_non_horizontal: bool,
):
    if source_grid is None:
        return obj, source_grid
    grid = source_grid if isinstance(source_grid, GaussianGrid) else parse_grid(source_grid)
    source_n = grid.n
    target_abs_lat = float(np.max(np.abs(target_latitudes)))
    source_abs_lat = _regular_kernel_safe_abs_lat(grid.n)
    if target_abs_lat <= source_abs_lat:
        return obj, grid
    promoted_n = source_n
    while target_abs_lat > _regular_kernel_safe_abs_lat(promoted_n):
        promoted_n *= 2
        if promoted_n > 4096:
            raise ValueError(
                "regular latitude/longitude target is too close to the poles "
                "for the available automatic Gaussian promotion limit"
            )
    promoted_grid = f"F{promoted_n}"
    promoted = regrid(
        obj,
        source_grid=grid.name,
        target_grid=promoted_grid,
        method="linear",
        variables=variables,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
        keep_attrs=keep_attrs,
        skip_non_horizontal=skip_non_horizontal,
    )
    return promoted, promoted_grid


def _regular_kernel_safe_abs_lat(n: int) -> float:
    """Maximum target latitude that keeps FULLPOS FPINT4/FPINT12 stencil valid."""
    lats = gaussian_latitudes(2 * int(n))
    return float(lats[2])


def _parse_regular_latlon_target(target_grid) -> dict | None:
    if isinstance(target_grid, str):
        text = target_grid.strip()
        upper = text.upper()
        if not upper.startswith("LL"):
            return None
        value = upper[2:].lstrip("_")
        if not value:
            raise ValueError(
                "regular latitude/longitude target must include a resolution, "
                "e.g. 'LL1.0'"
            )
        resolution = float(value)
        lat_step = lon_step = resolution
        name = f"LL{resolution:g}"
    elif isinstance(target_grid, dict):
        kind = str(target_grid.get("type", target_grid.get("grid_type", ""))).lower()
        if kind not in {"regular_ll", "latlon", "regular_latlon"}:
            return None
        resolution = target_grid.get("resolution")
        lat_step = float(target_grid.get("lat_step", resolution))
        lon_step = float(target_grid.get("lon_step", resolution))
        name = str(target_grid.get("name", f"LL{lat_step:g}x{lon_step:g}"))
    else:
        return None
    if lat_step <= 0.0 or lon_step <= 0.0:
        raise ValueError("regular latitude/longitude resolution must be positive")
    latitude = _regular_latitudes(lat_step)
    longitude = _regular_longitudes(lon_step)
    return {
        "name": name,
        "lat_step": lat_step,
        "lon_step": lon_step,
        "latitude": latitude,
        "longitude": longitude,
    }


def _regular_latitudes(step: float) -> np.ndarray:
    n_intervals = int(round(180.0 / step))
    if not np.isclose(n_intervals * step, 180.0, rtol=0.0, atol=1.0e-10):
        raise ValueError("regular latitude step must divide 180 degrees exactly")
    return 90.0 - (np.arange(n_intervals, dtype=np.float64) + 0.5) * step


def _regular_longitudes(step: float) -> np.ndarray:
    nlon = int(round(360.0 / step))
    if not np.isclose(nlon * step, 360.0, rtol=0.0, atol=1.0e-10):
        raise ValueError("regular longitude step must divide 360 degrees exactly")
    return np.arange(nlon, dtype=np.float64) * step


def _finalize_regular_latlon_output(
    out,
    *,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    source_grid,
    target_name: str,
    method: str,
    keep_attrs: bool,
):
    if isinstance(out, xr.Dataset):
        out = _assign_regular_latlon_dims(out, latitudes=latitudes, longitudes=longitudes)
        attrs = dict(out.attrs) if keep_attrs else {}
        attrs.update(_regular_ll_attrs(latitudes=latitudes, longitudes=longitudes))
        if keep_attrs:
            attrs = append_history(
                attrs,
                _regular_ll_history(
                    source_grid=source_grid,
                    target_name=target_name,
                    method=method,
                ),
            )
        out.attrs = attrs
        for name in out.data_vars:
            if {"latitude", "longitude"}.issubset(out[name].dims):
                out[name].attrs.update(
                    _regular_ll_attrs(latitudes=latitudes, longitudes=longitudes)
                )
        return out
    out = _assign_regular_latlon_dims(out, latitudes=latitudes, longitudes=longitudes)
    attrs = dict(out.attrs) if keep_attrs else {}
    attrs.update(_regular_ll_attrs(latitudes=latitudes, longitudes=longitudes))
    if keep_attrs:
        attrs = append_history(
            attrs,
            _regular_ll_history(source_grid=source_grid, target_name=target_name, method=method),
        )
    out.attrs = attrs
    return out


def _assign_regular_latlon_dims(out, *, latitudes: np.ndarray, longitudes: np.ndarray):
    rename = {}
    if "target_y" in out.dims:
        rename["target_y"] = "latitude"
    if "target_x" in out.dims:
        rename["target_x"] = "longitude"
    if rename:
        out = out.drop_vars(["latitude", "longitude"], errors="ignore").rename(rename)
    return out.assign_coords(latitude=latitudes, longitude=longitudes)


def _regular_ll_attrs(*, latitudes: np.ndarray, longitudes: np.ndarray) -> dict:
    return {
        "GRIB_gridType": "regular_ll",
        "GRIB_numberOfPoints": int(latitudes.size * longitudes.size),
        "GRIB_iDirectionIncrementInDegrees": (
            float(abs(longitudes[1] - longitudes[0]))
            if longitudes.size > 1
            else 360.0
        ),
        "GRIB_jDirectionIncrementInDegrees": (
            float(abs(latitudes[1] - latitudes[0]))
            if latitudes.size > 1
            else 180.0
        ),
    }


def _regular_ll_history(*, source_grid, target_name: str, method: str) -> str:
    return (
        f"fullpos regrid source_grid={_grid_name(source_grid)} "
        f"target_grid={target_name} method={method}"
    )

def _grid_name(grid) -> str:
    return getattr(grid, "name", str(grid))


def _copy_safe_coords(source: xr.Dataset, target: xr.Dataset) -> None:
    for name, coord in source.coords.items():
        if set(coord.dims).isdisjoint(_HORIZONTAL_DIM_NAMES):
            target.coords[name] = coord


def _has_horizontal_dims(obj: xr.DataArray, source_grid) -> bool:
    grid = source_grid if hasattr(source_grid, "is_reduced") else parse_grid(source_grid)
    if grid.is_reduced:
        return any(obj.sizes[dim] == grid.size for dim in obj.dims)
    return "latitude" in obj.dims and "longitude" in obj.dims
