from __future__ import annotations

from ..errors import FullposNotImplementedError
from .common import SUPPORTED_VERTICAL_TARGETS, normalize_vertical_target
from .eta import interpolate_to_eta
from .model_level import interpolate_to_model_levels
from .potential_temperature import interpolate_to_potential_temperature
from .potential_vorticity import interpolate_to_potential_vorticity
from .pressure import interpolate_to_pressure
from .temperature import interpolate_to_temperature


def vertical_interpolate(*args, target: str, **kwargs):
    """Dispatch FULLPOS-style vertical interpolation targets.

    The public target names mirror the main FULLPOS post-processing level
    families controlled through ``YOMFPC``. Implementations are introduced
    target-by-target so the package can keep a stable public API while native
    wrappers are added incrementally.
    """
    normalized = normalize_vertical_target(target)
    if normalized == "pressure":
        if not args and "values" not in kwargs:
            raise TypeError("vertical interpolation target 'pressure' requires values input")
        return interpolate_to_pressure(*args, **kwargs)
    if normalized == "model_level":
        if not args and "values" not in kwargs:
            raise TypeError("vertical interpolation target 'model_level' requires values input")
        return interpolate_to_model_levels(*args, **kwargs)
    if normalized == "potential_temperature":
        if not args and "values" not in kwargs:
            raise TypeError("vertical interpolation target 'potential_temperature' requires values input")
        return interpolate_to_potential_temperature(*args, **kwargs)
    if normalized == "temperature":
        if not args and "values" not in kwargs:
            raise TypeError("vertical interpolation target 'temperature' requires values input")
        return interpolate_to_temperature(*args, **kwargs)
    if normalized == "potential_vorticity":
        if not args and "values" not in kwargs:
            raise TypeError("vertical interpolation target 'potential_vorticity' requires values input")
        return interpolate_to_potential_vorticity(*args, **kwargs)
    if normalized == "eta":
        if not args and "values" not in kwargs:
            raise TypeError("vertical interpolation target 'eta' requires values input")
        return interpolate_to_eta(*args, **kwargs)
    raise FullposNotImplementedError(
        f"vertical interpolation target {normalized!r} is planned but not implemented"
    )


def vertical_capabilities() -> dict[str, str]:
    """Return current implementation status for vertical targets."""
    capabilities = {name: "planned" for name in SUPPORTED_VERTICAL_TARGETS}
    capabilities["model_level"] = "native_pp_chain"
    capabilities["potential_temperature"] = "native_pp_chain"
    capabilities["pressure"] = "native"
    capabilities["temperature"] = "native_pp_chain"
    capabilities["potential_vorticity"] = "native_ppltp"
    capabilities["eta"] = "native_ppleta"
    return capabilities
