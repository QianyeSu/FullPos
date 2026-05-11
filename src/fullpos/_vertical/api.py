from __future__ import annotations

from .common import SUPPORTED_VERTICAL_TARGETS
from .pos import dispatch_vertical_interpolate, iter_vertical_target_plans


def vertical_interpolate(*args, target: str, **kwargs):
    """Dispatch FULLPOS-style vertical interpolation targets.

    The public target names mirror the main FULLPOS post-processing level
    families controlled through ``YOMFPC``. This function only normalizes the
    target name and dispatches into the target-specific Python wrapper. The
    numerical kernels stay in the native FULLPOS/OpenIFS Fortran path; the
    Python layer only assembles and validates the request.
    """
    return dispatch_vertical_interpolate(*args, target=target, **kwargs)


def vertical_capabilities() -> dict[str, str]:
    """Return the current Python-to-native implementation status by target.

    The values describe whether each target is still planned, handled by a
    native FULLPOS wrapper, or routed through a specific Fortran backend.
    """
    capabilities = {name: "planned" for name in SUPPORTED_VERTICAL_TARGETS}
    for plan in iter_vertical_target_plans():
        capabilities[plan.target] = plan.capability
    return capabilities
