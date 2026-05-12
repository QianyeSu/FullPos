from __future__ import annotations

import platform
import sys

import numpy as np

from .spectral import native_backend_status
from .regridder import Regridder
from .vertical import vertical_capabilities


def backend_info() -> dict:
    """Return runtime diagnostics for the active fullpos backend."""

    info = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
    }
    info.update(native_backend_status())
    return info


def capabilities() -> dict:
    """Return a compact status summary for the current public feature set."""

    info = backend_info()
    return {
        "backend": info["backend"],
        "spectral": {
            "regrid": "native",
            "fit": "native",
            "synthesis": "native",
        },
        "filtering": {
            "spectral_filter": "native",
            "generic_spectral_filter": "native",
            "filter_matrix_io": "native",
        },
        "horizontal": {
            "regular_gaussian": "native",
            "reduced_gaussian": "native_special_cases",
            "requires_spectral_regrid_before_user_level_interpolation": True,
        },
        "vertical": vertical_capabilities(),
        "runtime": info,
    }


def doctor() -> dict:
    """Run lightweight runtime checks for the required native backend."""

    info = backend_info()
    checks = [
        _check("native prefix exists", bool(info["native_prefix_exists"]), info["native_prefix"]),
        _check("native module importable", bool(info["native_available"]), info.get("import_error")),
        _check(
            "native runtime directory exists",
            bool(info["native_runtime_dir_exists"]),
            info["native_runtime_dir"],
        ),
        _check(
            "required native libraries present",
            bool(info["required_native_libraries_present"]),
            _missing_libraries(info),
        ),
        _check(
            "external runtime libraries present",
            bool(info["external_runtime_libraries_present"]),
            _missing_external_libraries(info),
        ),
    ]
    checks.append(_native_smoke_check())
    info = backend_info()
    ok = all(check["ok"] for check in checks)
    return {"ok": ok, "backend_info": info, "checks": checks}


def _native_smoke_check() -> dict:
    try:
        values = np.ones((8, 16), dtype=np.float32)
        out = Regridder("N4", "N4", chunk_size=1).regrid_values(values)
        ok = out.shape == (8, 16) and bool(np.isfinite(out).all())
        detail = f"shape={out.shape}, finite={bool(np.isfinite(out).all())}"
    except Exception as exc:
        ok = False
        detail = f"{type(exc).__name__}: {exc}"
    return _check("native smoke regrid N4 -> N4", ok, detail)


def _missing_libraries(info: dict) -> str:
    missing = [
        name
        for name, present in info.get("required_native_libraries", {}).items()
        if not present
    ]
    if not missing:
        return "all required native libraries present"
    return "missing: " + ", ".join(missing)


def _missing_external_libraries(info: dict) -> str:
    libraries = info.get("external_runtime_libraries", {})
    if not libraries:
        return "no external runtime library list for this platform yet"
    missing = [name for name, path in libraries.items() if path is None]
    if not missing:
        return "all external runtime libraries found"
    return "missing: " + ", ".join(missing)


def _check(name: str, ok: bool, detail=None) -> dict:
    return {"name": name, "ok": bool(ok), "detail": detail}
