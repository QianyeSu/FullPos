from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

from .eta import interpolate_to_eta
from .height import (
    interpolate_to_flight_level,
    interpolate_to_height_above_orography,
    interpolate_to_height_above_sea,
)
from .model_level import interpolate_to_model_levels
from .potential_temperature import interpolate_to_potential_temperature
from .potential_vorticity import interpolate_to_potential_vorticity
from .pressure import interpolate_to_pressure
from .temperature import interpolate_to_temperature
from .common import normalize_vertical_target


@dataclass(frozen=True)
class VerticalTargetPlan:
    """Resolved native FULLPOS target and Python-side orchestration hint."""

    target: str
    backend: str
    capability: str
    handler: Callable[..., object]
    surface_handler: Callable[..., object] | None = None
    derived_handler: Callable[..., object] | None = None


_VERTICAL_TARGET_PLAN_REGISTRY: dict[str, VerticalTargetPlan] = {
    "model_level": VerticalTargetPlan("model_level", "PP_CHAIN", "native_pp_chain", interpolate_to_model_levels),
    "pressure": VerticalTargetPlan("pressure", "APACHE", "native", interpolate_to_pressure),
    "height_above_orography": VerticalTargetPlan(
        "height_above_orography",
        "GPHPRE/GPGEO/FPPS",
        "native_gpgeo_fpps",
        interpolate_to_height_above_orography,
    ),
    "height_above_sea": VerticalTargetPlan(
        "height_above_sea",
        "GPHPRE/GPGEO/FPPS",
        "native_gpgeo_fpps",
        interpolate_to_height_above_sea,
    ),
    "flight_level": VerticalTargetPlan(
        "flight_level",
        "GPHPRE/GPGEO/FPPS",
        "native_gpgeo_fpps",
        interpolate_to_flight_level,
    ),
    "potential_temperature": VerticalTargetPlan(
        "potential_temperature",
        "PP_CHAIN",
        "native_pp_chain",
        interpolate_to_potential_temperature,
    ),
    "potential_vorticity": VerticalTargetPlan(
        "potential_vorticity",
        "GPPVO/PPLTP",
        "native_ppltp",
        interpolate_to_potential_vorticity,
    ),
    "temperature": VerticalTargetPlan("temperature", "PP_CHAIN", "native_pp_chain", interpolate_to_temperature),
    "eta": VerticalTargetPlan("eta", "PPLETA", "native_ppleta", interpolate_to_eta),
}


def build_vertical_target_plan(target: str) -> VerticalTargetPlan:
    """Resolve a public target name into the native FULLPOS dispatch plan."""
    normalized = normalize_vertical_target(target)
    return _VERTICAL_TARGET_PLAN_REGISTRY[normalized]


def iter_vertical_target_plans() -> Iterable[VerticalTargetPlan]:
    """Yield the internal POS-style plan for each supported target."""
    yield from _VERTICAL_TARGET_PLAN_REGISTRY.values()


def get_vertical_target_plan(target: str) -> VerticalTargetPlan:
    """Resolve a public target name into the registered native plan."""
    return build_vertical_target_plan(target)


def vertical_target_capability(target: str) -> str:
    """Return the current capability tag for a public vertical target."""
    return get_vertical_target_plan(target).capability


def dispatch_vertical_interpolate(*args, target: str, **kwargs):
    """Dispatch a vertical request through the internal POS-style plan."""
    plan = get_vertical_target_plan(target)
    if not args and "values" not in kwargs:
        raise TypeError(f"vertical interpolation target {plan.target!r} requires values input")
    return plan.handler(*args, **kwargs)
