"""Top-level package for fullpos.

Public convenience exports are resolved lazily so ``import fullpos`` stays
lightweight and does not import compiled extensions until they are needed.
"""

from __future__ import annotations

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any

try:
    __version__ = version("fullpos")
except PackageNotFoundError:  # pragma: no cover - local source tree fallback
    __version__ = "0.1.0"

_OBJECT_EXPORTS = {
    "Regridder": ("fullpos.regridder", "Regridder"),
    "SpectralFilter": ("fullpos.filters", "SpectralFilter"),
    "average_interpolate": ("fullpos.interpolation", "average_interpolate"),
    "backend_info": ("fullpos.diagnostics", "backend_info"),
    "bilinear_interpolate": ("fullpos.interpolation", "bilinear_interpolate"),
    "capabilities": ("fullpos.diagnostics", "capabilities"),
    "doctor": ("fullpos.diagnostics", "doctor"),
    "diagnose_potential_vorticity": ("fullpos.vertical", "diagnose_potential_vorticity"),
    "filter_spectral_coefficients": ("fullpos.filters", "filter_spectral_coefficients"),
    "fpfilter_profile": ("fullpos.filters", "fpfilter_profile"),
    "gaussian_filter_profile": ("fullpos.filters", "gaussian_filter_profile"),
    "generic_spectral_filter": ("fullpos.filters", "generic_spectral_filter"),
    "horizontal_interpolate": ("fullpos.interpolation", "horizontal_interpolate"),
    "load_filter_matrix": ("fullpos.filters", "load_filter_matrix"),
    "load_spectral_filter": ("fullpos.filters", "load_spectral_filter"),
    "low_pass_filter_profile": ("fullpos.filters", "low_pass_filter_profile"),
    "nearest_interpolate": ("fullpos.interpolation", "nearest_interpolate"),
    "quadratic12_interpolate": ("fullpos.interpolation", "quadratic12_interpolate"),
    "regrid": ("fullpos.api", "regrid"),
    "regrid_values": ("fullpos.api", "regrid_values"),
    "save_filter_matrix": ("fullpos.filters", "save_filter_matrix"),
    "save_spectral_filter": ("fullpos.filters", "save_spectral_filter"),
    "spectral_filter": ("fullpos.filters", "spectral_filter"),
    "spectral_fit": ("fullpos.spectral", "spectral_fit"),
    "spectral_fit_values": ("fullpos.spectral", "spectral_fit_values"),
    "spectral_synthesis": ("fullpos.spectral", "spectral_synthesis"),
    "spectral_synthesis_values": ("fullpos.spectral", "spectral_synthesis_values"),
    "vertical_capabilities": ("fullpos.vertical", "vertical_capabilities"),
    "vertical_interpolate": ("fullpos.vertical", "vertical_interpolate"),
}

__all__ = sorted((*_OBJECT_EXPORTS, "__version__"))


def __getattr__(name: str) -> Any:
    if name in _OBJECT_EXPORTS:
        module_name, object_name = _OBJECT_EXPORTS[name]
        value = getattr(import_module(module_name), object_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
