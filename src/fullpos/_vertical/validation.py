from __future__ import annotations

from typing import Any

import numpy as np
import xarray as xr

from .pressure import (
    _infer_reference_pressure_from_ak,
    _midlevel_coefficients_from_half_levels,
    prepare_pressure_request,
)


def metric_block(reference: Any, candidate: Any) -> dict[str, float | int]:
    """Return finite-point error metrics between two arrays.

    NaNs and infinities are ignored pairwise. This is useful for pressure-level
    interpolation because low pressure levels can be invalid at some terrain
    points when extrapolation is disabled.
    """
    ref = np.asarray(reference, dtype=np.float64)
    cand = np.asarray(candidate, dtype=np.float64)
    if ref.shape != cand.shape:
        raise ValueError(f"metric inputs must have matching shapes: {ref.shape} vs {cand.shape}")
    mask = np.isfinite(ref) & np.isfinite(cand)
    count = int(mask.sum())
    if count == 0:
        return {
            "count": 0,
            "rmse": float("nan"),
            "mae": float("nan"),
            "max_abs": float("nan"),
            "bias": float("nan"),
        }
    diff = cand[mask] - ref[mask]
    return {
        "count": count,
        "rmse": float(np.sqrt(np.mean(diff**2))),
        "mae": float(np.mean(np.abs(diff))),
        "max_abs": float(np.max(np.abs(diff))),
        "bias": float(np.mean(diff)),
    }


def pressure_metric_summary(
    reference: xr.DataArray,
    candidate: xr.DataArray,
    *,
    level_dim: str | None = None,
) -> dict[str, Any]:
    """Compare two pressure-level fields and return overall/per-level metrics."""
    ref, cand = xr.align(reference, candidate, join="exact")
    dim = level_dim or _find_pressure_level_dim(ref)
    summary: dict[str, Any] = {
        "overall": metric_block(ref.values, cand.values),
        "per_level": {},
    }
    if dim is None:
        return summary

    for index, level in enumerate(np.asarray(ref[dim].values)):
        key = _format_level_key(level)
        summary["per_level"][key] = metric_block(
            ref.isel({dim: index}).values,
            cand.isel({dim: index}).values,
        )
    return summary


def skyborn_pressure_reference(
    data: xr.DataArray,
    *,
    levels,
    surface_pressure: xr.DataArray,
    hybrid_coefficients=None,
    method: str = "linear",
    p0: float | None = None,
) -> xr.DataArray:
    """Run Skyborn hybrid-to-pressure interpolation as a validation reference.

    This function is intentionally development-only. It imports Skyborn lazily
    so the FullPos runtime package still uses only the native FULLPOS/OpenIFS
    backend unless a validation tool explicitly asks for this reference path.
    """
    from skyborn.interp.interpolation import interp_hybrid_to_pressure

    request = prepare_pressure_request(
        data,
        levels=levels,
        surface_pressure=surface_pressure,
        hybrid_coefficients=hybrid_coefficients,
    )
    hyam, hybm = _build_midlevel_coefficients(
        data,
        hybrid_dim=request.hybrid_dim,
        ak=request.ak,
        bk=request.bk,
    )
    reference_p0 = _infer_reference_pressure_from_ak(request.ak) if p0 is None else float(p0)
    return interp_hybrid_to_pressure(
        data=data,
        ps=request.surface_pressure,
        hyam=hyam,
        hybm=hybm,
        p0=reference_p0,
        new_levels=request.levels,
        lev_dim=request.hybrid_dim,
        method=method,
    )


def skyborn_pressure_dataset_reference(
    dataset: xr.Dataset,
    *,
    levels,
    surface_pressure: xr.DataArray,
    variables=None,
    hybrid_coefficients=None,
    method: str = "linear",
    p0: float | None = None,
) -> xr.Dataset:
    """Run Skyborn pressure interpolation for selected Dataset variables."""
    names = list(dataset.data_vars) if variables is None else [str(name) for name in variables]
    missing = [name for name in names if name not in dataset.data_vars]
    if missing:
        raise KeyError(f"variables not found in dataset: {missing}")
    outputs = {
        name: skyborn_pressure_reference(
            dataset[name],
            levels=levels,
            surface_pressure=surface_pressure,
            hybrid_coefficients=hybrid_coefficients,
            method=method,
            p0=p0,
        )
        for name in names
    }
    out = xr.Dataset(outputs)
    out.attrs["reference_backend"] = "skyborn.interp.interpolation.interp_hybrid_to_pressure"
    out.attrs["reference_method"] = method
    out.attrs["reference_levels_pa"] = ",".join(f"{value:g}" for value in np.asarray(levels, dtype=np.float64))
    return out


def _build_midlevel_coefficients(
    model: xr.DataArray,
    *,
    hybrid_dim: str,
    ak: np.ndarray,
    bk: np.ndarray,
) -> tuple[xr.DataArray, xr.DataArray]:
    hyam_values, hybm_values = _midlevel_coefficients_from_half_levels(ak, bk)
    if hybrid_dim in model.coords:
        coords = {hybrid_dim: np.asarray(model.coords[hybrid_dim].values)}
    else:
        coords = {hybrid_dim: np.arange(hyam_values.size, dtype=np.int64)}
    hyam = xr.DataArray(hyam_values, dims=(hybrid_dim,), coords=coords, name="hyam")
    hybm = xr.DataArray(hybm_values, dims=(hybrid_dim,), coords=coords, name="hybm")
    return hyam, hybm


def _find_pressure_level_dim(obj: xr.DataArray) -> str | None:
    for dim in ("isobaricInhPa", "plev", "pressure", "level"):
        if dim in obj.dims:
            return dim
    return None


def _format_level_key(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:g}"
