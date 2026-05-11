from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import xarray as xr

from fullpos import vertical_interpolate
from fullpos._vertical.validation import pressure_metric_summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare native FULLPOS hybrid-to-pressure output against an "
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
        help="Target pressure levels in hPa, for example: 200 300 400 500",
    )
    parser.add_argument(
        "--surface-var",
        default="sp",
        help="Surface-pressure variable name",
    )
    parser.add_argument(
        "--chunks",
        nargs="*",
        default=["time=1", "values=10000"],
        help="Named chunks as dim=size entries, for example: time=1 values=10000",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional JSON output path for the statistics",
    )
    parser.add_argument(
        "--output-netcdf",
        type=Path,
        default=None,
        help="Optional NetCDF path for the native FULLPOS pressure output",
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


def _parse_chunks(entries: list[str]) -> dict[str, int]:
    chunks: dict[str, int] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"chunk entry must be dim=size, got {entry!r}")
        dim, value = entry.split("=", 1)
        dim = dim.strip()
        if not dim:
            raise ValueError(f"chunk entry has empty dimension: {entry!r}")
        size = int(value)
        if size <= 0:
            raise ValueError(f"chunk size must be positive: {entry!r}")
        chunks[dim] = size
    return chunks


def _open_dataset(path: Path, variables: list[str], type_of_level: str) -> xr.Dataset:
    return xr.Dataset(
        {name: _open_field(path, name, type_of_level).load() for name in variables}
    )


def _to_truth_level_coord(candidate: xr.DataArray, levels_hpa: np.ndarray) -> xr.DataArray:
    return candidate.assign_coords(pressure=("pressure", levels_hpa)).rename(pressure="isobaricInhPa")


def main() -> None:
    args = parse_args()
    variables = [str(name) for name in args.variables]
    levels_hpa = np.asarray(args.levels_hpa, dtype=np.float64)
    levels_pa = levels_hpa * 100.0
    chunks = _parse_chunks(args.chunks)

    model = _open_dataset(args.model, variables, "hybrid")
    surface = _open_field(args.surface, args.surface_var, "surface").load()
    native = vertical_interpolate(
        model,
        target="pressure",
        levels=levels_pa,
        variables=variables,
        chunks=chunks,
        surface_pressure=surface,
    )

    summary: dict[str, object] = {
        "backend": "FULLPOS",
        "levels_hpa": levels_hpa.tolist(),
        "chunks": chunks,
        "variables": {},
    }

    for short_name in variables:
        truth = _open_field(args.truth, short_name, "isobaricInhPa").sel(isobaricInhPa=levels_hpa).load()
        candidate = _to_truth_level_coord(native[short_name], levels_hpa)
        metrics = pressure_metric_summary(truth, candidate, level_dim="isobaricInhPa")
        overall = metrics["overall"]
        per_level = metrics["per_level"]
        summary["variables"][short_name] = {
            "overall": overall,
            "per_level": per_level,
        }
        print(
            f"{short_name}: overall rmse={overall['rmse']:.6g} "
            f"mae={overall['mae']:.6g} max_abs={overall['max_abs']:.6g} "
            f"bias={overall['bias']:.6g} count={overall['count']}"
        )
        for level in levels_hpa:
            key = f"{int(level):d}" if float(level).is_integer() else f"{level:g}"
            block = per_level[key]
            print(
                f"{short_name}: {level:g} hPa rmse={block['rmse']:.6g} "
                f"mae={block['mae']:.6g} max_abs={block['max_abs']:.6g} "
                f"bias={block['bias']:.6g} count={block['count']}"
            )

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(
            json.dumps(summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"output_json={args.output_json}")
    if args.output_netcdf is not None:
        args.output_netcdf.parent.mkdir(parents=True, exist_ok=True)
        native.to_netcdf(args.output_netcdf)
        print(f"output_netcdf={args.output_netcdf}")

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
