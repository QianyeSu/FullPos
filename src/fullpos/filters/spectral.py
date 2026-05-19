from __future__ import annotations

import numpy as np
import xarray as xr

from ..metadata import append_history, format_regrid_history, resolve_source_grid
from ..interpolation import DEFAULT_CHUNK_SIZE, Regridder
from ..spectral import spectral_fit, spectral_synthesis
from .core import SpectralFilter
from .profiles import expand_profile_to_coefficients, fpfilter_profile, infer_ntrunc_from_nspec2


def spectral_filter(
    obj,
    *,
    grid=None,
    ntrunc: int,
    variables=None,
    axis=-1,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    keep_attrs=True,
    skip_non_horizontal=True,
):
    """Apply triangular spectral truncation and return grid-point values.

    This is equivalent to regridding a field from a grid back to itself with a
    lower ``ntrunc``. It uses the native ECTRANS spectral path.
    """
    if ntrunc is None:
        raise ValueError("ntrunc is required for spectral_filter")

    if isinstance(obj, xr.DataArray):
        resolved_grid = resolve_source_grid(obj, grid)
        regridder = Regridder.from_dataarray(
            obj,
            target_grid=resolved_grid,
            source_grid=resolved_grid,
            ntrunc=ntrunc,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
        )
        return regridder.regrid_data_array(obj, keep_attrs=keep_attrs)

    if isinstance(obj, xr.Dataset):
        selected = list(obj.data_vars) if variables is None else [str(v) for v in variables]
        missing = [name for name in selected if name not in obj.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        out = xr.Dataset(attrs=dict(obj.attrs))
        for name, data_array in obj.data_vars.items():
            if name not in selected:
                out[name] = data_array
                continue
            try:
                resolved_grid = resolve_source_grid(data_array, grid)
                regridder = Regridder.from_dataarray(
                    data_array,
                    target_grid=resolved_grid,
                    source_grid=resolved_grid,
                    ntrunc=ntrunc,
                    chunk_size=chunk_size,
                    missing_policy=missing_policy,
                )
                out[name] = regridder.regrid_data_array(data_array, keep_attrs=keep_attrs)
            except ValueError:
                if variables is None and skip_non_horizontal:
                    out[name] = data_array
                    continue
                raise
        out.attrs = append_history(
            out.attrs,
            format_regrid_history(
                source_grid=grid or "inferred",
                target_grid=grid or "inferred",
                method="spectral_filter",
                ntrunc=ntrunc,
                chunk_size=chunk_size,
                variables=selected,
            ),
        )
        return out

    regridder = Regridder(
        grid,
        grid,
        ntrunc=ntrunc,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
    )
    return regridder.regrid_values(obj, axis=axis)


def generic_spectral_filter(
    values,
    *,
    grid,
    filter_kind: str = "gaussian",
    ntrunc: int | None = None,
    cutoff: int | None = None,
    axis=-1,
    chunk_size=DEFAULT_CHUNK_SIZE,
    missing_policy="error",
    selectivity: float = 1.0,
    low_pass_exponent: float = 1.0,
    gaussian_exponent: float = 1.0,
) -> np.ndarray:
    """Apply a FULLPOS-style diagonal spectral filter profile to grid values.

    The operation is ``grid -> spectral coefficients -> coefficient weighting
    -> grid``. ``filter_kind`` may be ``"gaussian"``, ``"low_pass"``, or a
    prebuilt :class:`fullpos.filters.SpectralFilter`.
    """
    if isinstance(filter_kind, SpectralFilter):
        return filter_kind.apply(
            values,
            grid=grid,
            axis=axis,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
        )
    coeffs = spectral_fit(
        values,
        grid=grid,
        ntrunc=ntrunc,
        axis=axis,
        chunk_size=chunk_size,
        missing_policy=missing_policy,
    )
    filtered = filter_spectral_coefficients(
        coeffs,
        filter_kind=filter_kind,
        ntrunc=ntrunc,
        cutoff=cutoff,
        selectivity=selectivity,
        low_pass_exponent=low_pass_exponent,
        gaussian_exponent=gaussian_exponent,
    )
    return spectral_synthesis(
        filtered,
        grid=grid,
        ntrunc=ntrunc,
        axis=-1,
        chunk_size=chunk_size,
    )


def filter_spectral_coefficients(
    coefficients,
    *,
    filter_kind: str = "gaussian",
    ntrunc: int | None = None,
    cutoff: int | None = None,
    axis: int = -1,
    selectivity: float = 1.0,
    low_pass_exponent: float = 1.0,
    gaussian_exponent: float = 1.0,
) -> np.ndarray:
    """Apply a FULLPOS-style profile directly to ECTRANS spectral coefficients."""
    arr = np.asarray(coefficients)
    if arr.ndim == 0:
        raise ValueError("coefficients must have at least one dimension")
    coeff_axis = int(axis)
    if coeff_axis < 0:
        coeff_axis += arr.ndim
    if coeff_axis < 0 or coeff_axis >= arr.ndim:
        raise ValueError("axis is out of bounds")
    nspec2 = arr.shape[coeff_axis]
    resolved_ntrunc = infer_ntrunc_from_nspec2(nspec2) if ntrunc is None else int(ntrunc)
    expected = (resolved_ntrunc + 1) * (resolved_ntrunc + 2)
    if nspec2 != expected:
        raise ValueError(
            f"coefficients have {nspec2} values, but T{resolved_ntrunc} expects {expected}"
        )
    profile = fpfilter_profile(
        resolved_ntrunc,
        kind=filter_kind,
        cutoff=cutoff,
        selectivity=selectivity,
        low_pass_exponent=low_pass_exponent,
        gaussian_exponent=gaussian_exponent,
    )
    weights = expand_profile_to_coefficients(profile)
    shape = [1] * arr.ndim
    shape[coeff_axis] = weights.size
    return arr * weights.reshape(shape)
