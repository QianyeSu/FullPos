from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


_NATIVE_LIB_STEMS = (
    "transi_dp",
    "trans_dp",
    "fiat",
    "parkind_dp",
    "parkind_sp",
)
_WINDOWS_EXTERNAL_RUNTIME_LIBRARIES = (
    "openblas.dll",
    "libgfortran-5.dll",
    "libgomp-1.dll",
    "libstdc++-6.dll",
    "libwinpthread-1.dll",
    "libgcc_s_seh-1.dll",
)


def native_runtime_dir() -> Path:
    """Return the platform-specific directory containing native runtime libraries."""
    prefix = native_prefix()
    if sys.platform == "win32":
        return prefix / "bin"
    return prefix / "lib"


def native_prefix() -> Path:
    """Return configured or bundled native ECTRANS/FIAT installation prefix."""
    override = os.environ.get("FULLPOS_NATIVE_PREFIX")
    if override:
        return Path(override)
    configured = _configured_native_prefix()
    if configured:
        return configured
    root = Path(__file__).resolve().parents[2]
    return root / "extern" / "fullpos" / "local"


def native_library_names() -> tuple[str, ...]:
    """Return native library filenames required by the package runtime."""
    if sys.platform == "win32":
        return tuple(f"lib{name}.dll" for name in _NATIVE_LIB_STEMS)
    if sys.platform == "darwin":
        return tuple(f"lib{name}.dylib" for name in _NATIVE_LIB_STEMS)
    return tuple(f"lib{name}.so" for name in _NATIVE_LIB_STEMS)


def native_library_status() -> dict[str, bool]:
    """Return whether each required native library exists in the runtime directory."""
    directory = native_runtime_dir()
    return {name: (directory / name).exists() for name in native_library_names()}


def external_runtime_library_names() -> tuple[str, ...]:
    """Return external toolchain runtime libraries checked by diagnostics."""
    if sys.platform == "win32":
        return _WINDOWS_EXTERNAL_RUNTIME_LIBRARIES
    # Linux/macOS names are toolchain-dependent; keep this explicit until those builds exist.
    return ()


def external_runtime_library_status() -> dict[str, str | None]:
    """Return discovered paths for external runtime libraries."""
    return {name: _find_runtime_library(name) for name in external_runtime_library_names()}


def add_native_runtime_dir() -> None:
    """Add the native runtime directory to Windows DLL search paths when needed."""
    if sys.platform != "win32" or not hasattr(os, "add_dll_directory"):
        return
    directory = native_runtime_dir()
    if directory.exists():
        os.add_dll_directory(str(directory))


def native_runtime_info() -> dict:
    """Return runtime path and native library availability diagnostics."""
    prefix = native_prefix()
    directory = native_runtime_dir()
    libraries = native_library_status()
    external_libraries = external_runtime_library_status()
    return {
        "native_runtime_platform": sys.platform,
        "native_prefix": str(prefix),
        "native_prefix_exists": prefix.exists(),
        "native_runtime_dir": str(directory),
        "native_runtime_dir_exists": directory.exists(),
        "required_native_libraries": libraries,
        "required_native_libraries_present": all(libraries.values()),
        "external_runtime_libraries": external_libraries,
        "external_runtime_libraries_present": all(
            path is not None for path in external_libraries.values()
        ),
    }


def _configured_native_prefix() -> Path | None:
    try:
        from ._build_config import NATIVE_PREFIX
    except Exception:
        return None
    if not NATIVE_PREFIX:
        return None
    return Path(_normalize_configured_path(NATIVE_PREFIX))


def _normalize_configured_path(path: str) -> str:
    # Meson-python currently writes Windows paths into a generated Python string.
    # If backslashes were not escaped, sequences such as \f become control chars.
    return str(path).replace("\x0c", "\\f")


def _find_runtime_library(name: str) -> str | None:
    local = native_runtime_dir() / name
    if local.exists():
        return str(local)
    found = shutil.which(name)
    if found is not None:
        return found
    for directory in _path_entries():
        candidate = directory / name
        if candidate.exists():
            return str(candidate)
    return None


def _path_entries() -> list[Path]:
    entries = [Path(entry) for entry in os.environ.get("PATH", "").split(os.pathsep) if entry]
    if sys.platform == "win32":
        entries.append(Path(sys.prefix) / "Library" / "bin")
    return entries
