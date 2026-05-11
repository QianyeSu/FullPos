from __future__ import annotations


SUPPORTED_VERTICAL_TARGETS = (
    "model_level",
    "pressure",
    "height_above_orography",
    "height_above_sea",
    "flight_level",
    "potential_temperature",
    "potential_vorticity",
    "temperature",
    "eta",
)


def normalize_vertical_target(target: str) -> str:
    """Normalize a public vertical target name used by the FULLPOS dispatch.

    This is a Python-side gate only; native backend handling is chosen later
    by the target-specific wrapper.
    """
    normalized = str(target).lower()
    if normalized not in SUPPORTED_VERTICAL_TARGETS:
        raise ValueError(
            "target must be one of: " + ", ".join(SUPPORTED_VERTICAL_TARGETS)
        )
    return normalized
