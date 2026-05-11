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
            "Compare Skyborn hybrid-to-pressure reference output against an "
            "official pressure-level ERA5/OpenIFS field on the same grid."
        )
    )
    parser.add_argument("model", type=Path, help="Path to model-level GRIB/NetCDF input")
    parser.add_argument("surface", type=Path, help="Path to surface-pressure GRIB/NetCDF input")
    parser.add_argument("truth", type=Path, help="Path to pressure-level GRIB/NetCDF truth input")
    parser.add_argument(
        "--variables",
        nargs="+",
        default=["t"],
        help="Variables to compare, for example: t u v q",
    )
    parser.add_argument(
        "--levels-hpa",
        nargs="+",
        type=float,
        required=True,
        help="Target pressure levels in hPa, for example: 300 400 500",
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


def main() -> None:
    args = parse_args()
    levels_hpa = np.asarray(args.levels_hpa, dtype=np.float64)
    levels_pa = levels_hpa * 100.0
    first_var = str(args.variables[0])

    first_model = _open_field(args.model, first_var, "hybrid").load()
    surface = _open_field(args.surface, args.surface_var, "surface").load()
    request = prepare_pressure_request(
        first_model,
        levels=levels_pa,
        surface_pressure=surface,
        hybrid_coefficients=args.coeffs,
    )
    p0 = _infer_reference_pressure_from_ak(request.ak)

    summary: dict[str, object] = {
        "method": args.method,
        "p0": float(p0),
        "levels_hpa": levels_hpa.tolist(),
        "variables": {},
    }

    for short_name in [str(name) for name in args.variables]:
        model = _open_field(args.model, short_name, "hybrid").load()
        truth = _open_field(args.truth, short_name, "isobaricInhPa").sel(isobaricInhPa=levels_hpa).load()
        reference = skyborn_pressure_reference(
            model,
            levels=request.levels,
            surface_pressure=request.surface_pressure,
            hybrid_coefficients=args.coeffs,
            method=args.method,
            p0=p0,
        )
        reference = reference.assign_coords(plev=("plev", levels_hpa)).rename(plev="isobaricInhPa")
        metrics = pressure_metric_summary(truth, reference, level_dim="isobaricInhPa")
        overall = metrics["overall"]
        per_level = metrics["per_level"]
        summary["variables"][short_name] = {
            "overall": overall,
            "per_level": per_level,
        }
        print(
            f"{short_name}: overall rmse={overall['rmse']:.6g} "
            f"mae={overall['mae']:.6g} max_abs={overall['max_abs']:.6g} count={overall['count']}"
        )
        for level in levels_hpa:
            block = per_level[f"{int(level):d}"]
            print(
                f"{short_name}: {int(level)} hPa rmse={block['rmse']:.6g} "
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
