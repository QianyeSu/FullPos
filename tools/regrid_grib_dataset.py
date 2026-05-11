from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import xarray as xr

from fullpos import regrid


def main() -> None:
    args = _parse_args()
    reports = []
    for variable in args.variables:
        try:
            reports.append(_process_variable(args, variable))
        except Exception as exc:
            report = {
                "variable": variable,
                "ok": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
            reports.append(report)
            if args.strict:
                _print_report({"path": str(args.path), "results": reports}, args.json)
                raise SystemExit(1) from exc

    summary = {
        "path": str(args.path),
        "target_grid": args.target_grid,
        "method": args.method,
        "chunk_size": args.chunk_size,
        "missing_policy": args.missing_policy,
        "results": reports,
    }
    _print_report(summary, args.json)
    if args.strict and not all(row["ok"] for row in reports):
        raise SystemExit(1)


def _process_variable(args: argparse.Namespace, variable: str) -> dict:
    open_start = perf_counter()
    ds = xr.open_dataset(
        args.path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {
                "shortName": variable,
                "typeOfLevel": args.type_of_level,
            },
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    if variable not in ds:
        raise KeyError(f"{variable!r} not found after cfgrib open; data_vars={list(ds.data_vars)}")
    obj = ds[variable]
    if args.time_index is not None and "time" in obj.dims:
        obj = obj.isel(time=args.time_index)
    if args.max_levels is not None and args.level_dim in obj.dims:
        obj = obj.isel({args.level_dim: slice(0, args.max_levels)})
    open_s = perf_counter() - open_start

    regrid_start = perf_counter()
    out = regrid(
        obj,
        target_grid=args.target_grid,
        method=args.method,
        chunk_size=args.chunk_size,
        missing_policy=args.missing_policy,
    )
    regrid_s = perf_counter() - regrid_start

    input_values = obj.values
    input_finite = np.isfinite(input_values)
    values = out.values
    output_finite = np.isfinite(values)
    return {
        "variable": variable,
        "ok": True,
        "input_dims": list(obj.dims),
        "input_shape": list(obj.shape),
        "output_dims": list(out.dims),
        "output_shape": list(out.shape),
        "source_grid_type": obj.attrs.get("GRIB_gridType"),
        "source_n": obj.attrs.get("GRIB_N"),
        "target_grid_type": out.attrs.get("GRIB_gridType"),
        "target_n": out.attrs.get("GRIB_N"),
        "target_points": out.attrs.get("GRIB_numberOfPoints"),
        "input_finite_count": int(input_finite.sum()),
        "input_size": int(input_values.size),
        "input_has_missing": bool(not input_finite.all()),
        "output_finite_count": int(output_finite.sum()),
        "output_size": int(values.size),
        "finite": bool(output_finite.all()),
        "min": _nan_stat(values, np.nanmin),
        "max": _nan_stat(values, np.nanmax),
        "open_s": open_s,
        "regrid_s": regrid_s,
    }


def _parse_chunk_size(value: str) -> int | None:
    if value.lower() in {"none", "all"}:
        return None
    chunk_size = int(value)
    if chunk_size <= 0:
        raise argparse.ArgumentTypeError("chunk size must be positive or 'none'")
    return chunk_size


def _nan_stat(values, func) -> float | None:
    if not np.isfinite(values).any():
        return None
    return float(func(values))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Regrid selected variables from a GRIB file.")
    parser.add_argument("path", type=Path)
    parser.add_argument("--variables", nargs="+", required=True)
    parser.add_argument("--type-of-level", required=True)
    parser.add_argument("--target-grid", default="O160")
    parser.add_argument("--method", choices=["linear", "spectral", "masked"], default="linear")
    parser.add_argument("--chunk-size", type=_parse_chunk_size, default=64)
    parser.add_argument("--missing-policy", choices=["error", "ignore"], default="error")
    parser.add_argument("--time-index", type=int, default=0)
    parser.add_argument("--max-levels", type=int, default=None)
    parser.add_argument("--level-dim", default="hybrid")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def _print_report(report: dict, as_json: bool) -> None:
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print(f"path: {report['path']}")
    print(f"target_grid: {report['target_grid']}")
    print(f"method: {report['method']}")
    print(f"chunk_size: {report['chunk_size']}")
    print(f"missing_policy: {report['missing_policy']}")
    for row in report["results"]:
        if not row["ok"]:
            print(f"{row['variable']}: FAIL {row['error']}")
            continue
        print(
            f"{row['variable']}: "
            f"{row['input_shape']} -> {row['output_shape']} "
            f"source={row['source_grid_type']} N{row['source_n']} "
            f"target={row['target_grid_type']} N{row['target_n']} "
            f"input_finite={row['input_finite_count']}/{row['input_size']} "
            f"output_finite={row['output_finite_count']}/{row['output_size']} "
            f"open_s={row['open_s']:.3f} "
            f"regrid_s={row['regrid_s']:.3f} "
            f"min={_format_float(row['min'])} max={_format_float(row['max'])}"
        )


def _format_float(value: float | None) -> str:
    if value is None:
        return "None"
    return f"{value:.6g}"


if __name__ == "__main__":
    main()
