from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FULLPOS_EXTERN = ROOT / "extern" / "fullpos"
FIAT_SOURCE = FULLPOS_EXTERN / "build-src" / "contrib" / "fiat"
ECTRANS_SOURCE = FULLPOS_EXTERN / "build-src" / "contrib" / "ectrans"


def main() -> None:
    args = _parse_args()
    prefix = _resolve(args.prefix)
    build_root = _resolve(args.build_root)
    ecbuild_dir = _resolve(args.ecbuild_dir)

    if args.fetch_ecbuild:
        _ensure_ecbuild(ecbuild_dir, args.ecbuild_ref)
    _require_ecbuild(ecbuild_dir)

    if args.clean:
        _safe_rmtree(build_root)
    build_root.mkdir(parents=True, exist_ok=True)
    prefix.mkdir(parents=True, exist_ok=True)

    cmake_prefix_path = _cmake_prefix_path([ecbuild_dir], os.environ.get("CMAKE_PREFIX_PATH"))

    fiat_build = build_root / "fiat"
    _configure_build_install(
        FIAT_SOURCE,
        fiat_build,
        prefix,
        cmake_prefix_path=cmake_prefix_path,
        extra_args=[
            "-DBUILD_SHARED_LIBS=ON",
            "-DENABLE_MPI=OFF",
            "-DENABLE_OMP=ON",
            "-DENABLE_TESTS=OFF",
            "-DENABLE_FORTRAN_C_INTERFACE=OFF",
            "-DENABLE_RPATHS=ON",
            "-DENABLE_RELATIVE_RPATHS=ON",
        ],
    )

    ectrans_build = build_root / "ectrans"
    ectrans_prefix_path = _cmake_prefix_path([prefix, ecbuild_dir], os.environ.get("CMAKE_PREFIX_PATH"))
    _configure_build_install(
        ECTRANS_SOURCE,
        ectrans_build,
        prefix,
        cmake_prefix_path=ectrans_prefix_path,
        extra_args=[
            "-DBUILD_SHARED_LIBS=ON",
            "-DENABLE_MPI=OFF",
            "-DENABLE_OMP=ON",
            "-DENABLE_TESTS=OFF",
            "-DENABLE_DOUBLE_PRECISION=ON",
            "-DENABLE_SINGLE_PRECISION=OFF",
            "-DENABLE_TRANSI=ON",
            "-DENABLE_FFTW=OFF",
            "-DENABLE_MKL=OFF",
            "-DENABLE_FORTRAN_C_INTERFACE=OFF",
            "-DENABLE_RPATHS=ON",
            "-DENABLE_RELATIVE_RPATHS=ON",
        ],
    )

    print(f"Native prefix ready: {prefix}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the FIAT/ECTRANS native prefix used by FullPos wheels.")
    parser.add_argument(
        "--prefix",
        default=str(FULLPOS_EXTERN / "local"),
        help="Install prefix for FIAT/ECTRANS. Defaults to extern/fullpos/local.",
    )
    parser.add_argument(
        "--build-root",
        default=str(ROOT / ".native-build"),
        help="Temporary CMake build root.",
    )
    parser.add_argument(
        "--ecbuild-dir",
        default=str(ROOT / ".ci" / "ecbuild"),
        help="Local ecbuild checkout or installation prefix.",
    )
    parser.add_argument(
        "--fetch-ecbuild",
        action="store_true",
        help="Clone ECMWF ecbuild into --ecbuild-dir when it is missing.",
    )
    parser.add_argument(
        "--ecbuild-ref",
        default="develop",
        help="Git ref used when --fetch-ecbuild clones ecbuild.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove --build-root before configuring.",
    )
    return parser.parse_args()


def _resolve(path: str) -> Path:
    return Path(path).expanduser().resolve()


def _ensure_ecbuild(ecbuild_dir: Path, ref: str) -> None:
    if _is_ecbuild_tree(ecbuild_dir):
        return
    if ecbuild_dir.exists():
        _safe_rmtree(ecbuild_dir)
    ecbuild_dir.parent.mkdir(parents=True, exist_ok=True)
    clone = [
        "git",
        "clone",
        "--depth",
        "1",
        "--branch",
        ref,
        "https://github.com/ecmwf/ecbuild.git",
        str(ecbuild_dir),
    ]
    try:
        _run(clone)
    except subprocess.CalledProcessError:
        _safe_rmtree(ecbuild_dir)
        _run(
            [
                "git",
                "clone",
                "--depth",
                "1",
                "https://github.com/ecmwf/ecbuild.git",
                str(ecbuild_dir),
            ]
        )


def _require_ecbuild(ecbuild_dir: Path) -> None:
    if not _is_ecbuild_tree(ecbuild_dir):
        raise SystemExit(
            "ecbuild was not found. Pass --ecbuild-dir or use --fetch-ecbuild."
        )


def _is_ecbuild_tree(path: Path) -> bool:
    return (
        (path / "cmake" / "ecbuild-config.cmake").exists()
        or (path / "share" / "ecbuild" / "cmake" / "ecbuild-config.cmake").exists()
    )


def _configure_build_install(
    source: Path,
    build: Path,
    prefix: Path,
    *,
    cmake_prefix_path: str,
    extra_args: list[str],
) -> None:
    configure = [
        "cmake",
        "-S",
        str(source),
        "-B",
        str(build),
        "-G",
        "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DCMAKE_INSTALL_PREFIX={prefix}",
        f"-DCMAKE_PREFIX_PATH={cmake_prefix_path}",
        *extra_args,
    ]
    _run(configure)
    _run(["cmake", "--build", str(build), "--config", "Release", "--parallel"])
    _run(["cmake", "--install", str(build), "--config", "Release"])


def _cmake_prefix_path(prefixes: list[Path], existing: str | None) -> str:
    values = [str(path) for path in prefixes]
    if existing:
        values.append(existing)
    return ";".join(values)


def _safe_rmtree(path: Path) -> None:
    if not path.exists():
        return
    root = ROOT.resolve()
    resolved = path.resolve()
    if resolved == root or root not in resolved.parents:
        raise SystemExit(f"Refusing to remove path outside the repository: {resolved}")
    shutil.rmtree(resolved)


def _run(command: list[str]) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
