from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr

from fullpos._vertical.pressure import _infer_reference_pressure_from_ak, prepare_pressure_request
from fullpos._vertical.validation import pressure_metric_summary, skyborn_pressure_reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare a pressure-level candidate output against Skyborn's "
            "compiled hybrid-to-pressure reference on the same input."
        )
    )
    parser.add_argument("model", type=Path, help="Path to model-level GRIB/NetCDF input")
    parser.add_argument("surface", type=Path, help="Path to surface-pressure GRIB/NetCDF input")
    parser.add_argument("candidate", type=Path, help="Path to pressure-level candidate NetCDF/GRIB input")
    parser.add_argument(
        "--variables",
        nargs="+",
        default=["t"],
        help="Variables to compare, for example: t u v q",
    )
    parser.add_argument(
        "--levels",
        nargs="+",
        type=float,
        required=True,
        help="Target pressure levels in Pa",
    )
    parser.add_argument(
        "--method",
        default="linear",
        choices=("linear", "log", "log-log"),
        help="Skyborn vertical interpolation method",
    )
    parser.add_argument(
        "--coeffs",
        type=Path,
        default=None,
        help="Optional external hybrid coefficient dataset",
    )
    parser.add_argument(
        "--surface-var",
        default="sp",
        help="Surface-pressure variable name",
    )
    parser.add_argument(
        "--candidate-level-dim",
        default=None,
        help="Optional pressure-level dimension name in the candidate file",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional JSON output path for the statistics",
    )
    parser.add_argument(
        "--max-rmse",
        type=float,
        default=None,
        help="Optional failure threshold for any variable overall RMSE",
    )
    parser.add_argument(
        "--max-max-abs",
        type=float,
        default=None,
        help="Optional failure threshold for any variable overall max_abs",
    )
    return parser.parse_args()


def _open_field(path: Path, short_name: str, type_of_level: str) -> xr.DataArray:
    if path.suffix.lower() in {".nc", ".nc4"}:
        return xr.open_dataset(path)[short_name]
    read_keys = ["gridType", "N", "pl", "numberOfPoints", "packingType"]
    if type_of_level == "hybrid":
        read_keys.append("pv")
    return xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": short_name, "typeOfLevel": type_of_level},
            "read_keys": read_keys,
        },
    )[short_name]


def _open_candidate(path: Path, short_name: str) -> xr.DataArray:
    if path.suffix.lower() in {".nc", ".nc4"}:
        return xr.open_dataset(path)[short_name]
    return xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": short_name, "typeOfLevel": "isobaricInhPa"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )[short_name]


def _standardize_candidate_levels(
    candidate: xr.DataArray,
    *,
    levels_pa: np.ndarray,
    level_dim: str | None,
) -> xr.DataArray:
    if "plev" in candidate.dims:
        return candidate.assign_coords(plev=("plev", levels_pa))
    dim = level_dim or _find_level_dim(candidate)
    if dim is None:
        raise ValueError("candidate output does not contain a recognizable pressure-level dimension")
    if candidate.sizes[dim] != levels_pa.size:
        raise ValueError(
            f"candidate level dimension {dim!r} has size {candidate.sizes[dim]}, expected {levels_pa.size}"
        )
    out = candidate.rename({dim: "plev"})
    return out.assign_coords(plev=("plev", levels_pa))


def _find_level_dim(obj: xr.DataArray) -> str | None:
    for dim in ("isobaricInhPa", "plev", "pressure", "level"):
        if dim in obj.dims:
            return dim
    return None


def main() -> None:
    args = parse_args()
    levels_pa = np.asarray(args.levels, dtype=np.float64)
    variable_names = [str(name) for name in args.variables]

    first = _open_field(args.model, variable_names[0], "hybrid").load()
    surface = _open_field(args.surface, args.surface_var, "surface").load()
    request = prepare_pressure_request(
        first,
        levels=levels_pa,
        surface_pressure=surface,
        hybrid_coefficients=args.coeffs,
    )
    p0 = _infer_reference_pressure_from_ak(request.ak)

    summary: dict[str, object] = {
        "reference_backend": "skyborn.interp.interpolation.interp_hybrid_to_pressure",
        "method": args.method,
        "p0": float(p0),
        "levels_pa": request.levels.tolist(),
        "variables": {},
    }

    for name in variable_names:
        model = _open_field(args.model, name, "hybrid").load()
        reference = skyborn_pressure_reference(
            model,
            levels=request.levels,
            surface_pressure=request.surface_pressure,
            hybrid_coefficients=args.coeffs,
            method=args.method,
            p0=p0,
        )
        candidate = _standardize_candidate_levels(
            _open_candidate(args.candidate, name).load(),
            levels_pa=request.levels,
            level_dim=args.candidate_level_dim,
        )
        metrics = pressure_metric_summary(reference, candidate, level_dim="plev")
        summary["variables"][name] = metrics
        overall = metrics["overall"]
        print(
            f"{name}: overall rmse={overall['rmse']:.6g} "
            f"mae={overall['mae']:.6g} max_abs={overall['max_abs']:.6g} count={overall['count']}"
        )
        for level, block in metrics["per_level"].items():
            print(
                f"{name}: {level} Pa rmse={block['rmse']:.6g} "
                f"mae={block['mae']:.6g} max_abs={block['max_abs']:.6g} count={block['count']}"
            )

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"output_json={args.output_json}")

    failures: list[str] = []
    for name, result in summary["variables"].items():
        overall = result["overall"]
        if args.max_rmse is not None and overall["rmse"] > args.max_rmse:
            failures.append(f"{name} rmse {overall['rmse']:.6g} > {args.max_rmse:.6g}")
        if args.max_max_abs is not None and overall["max_abs"] > args.max_max_abs:
            failures.append(f"{name} max_abs {overall['max_abs']:.6g} > {args.max_max_abs:.6g}")
    if failures:
        raise SystemExit("threshold check failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
