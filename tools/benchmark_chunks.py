from __future__ import annotations

import argparse
import json
from pathlib import Path
from statistics import mean, stdev
from time import perf_counter

import numpy as np
import xarray as xr

from fullpos import Regridder
from fullpos.grids import infer_grid_name_from_attrs


DEFAULT_PATH = Path(
    r"L:\ERA5_test\era5_reanalysis_model_level_20250102_packing_CCSDS_O320.grib2"
)


def main() -> None:
    args = _parse_args()
    open_start = perf_counter()
    field = _open_field(args.path, args.short_name, args.type_of_level, args.time_index)
    if args.max_levels is not None and "hybrid" in field.dims:
        field = field.isel(hybrid=slice(0, args.max_levels))
    field.load()
    open_s = perf_counter() - open_start

    source_grid = args.source_grid or infer_grid_name_from_attrs(field.attrs)
    results = []
    for chunk_size in args.chunk_size:
        regridder = Regridder(source_grid, args.target_grid, chunk_size=chunk_size)
        timings = []
        shape = None
        finite = None
        for _ in range(args.repeat):
            start = perf_counter()
            out = regridder.regrid_data_array(field)
            elapsed = perf_counter() - start
            timings.append(elapsed)
            shape = list(out.shape)
            finite = bool(np.isfinite(out.values).all())
        results.append(
            {
                "chunk_size": chunk_size,
                "repeat": args.repeat,
                "mean_s": mean(timings),
                "stdev_s": stdev(timings) if len(timings) > 1 else 0.0,
                "min_s": min(timings),
                "max_s": max(timings),
                "output_shape": shape,
                "finite": finite,
            }
        )

    report = {
        "path": str(args.path),
        "variable": args.short_name,
        "source_grid": source_grid,
        "target_grid": args.target_grid,
        "input_shape": list(field.shape),
        "open_and_load_s": open_s,
        "results": results,
    }
    _print_report(report, as_json=args.json)


def _open_field(path: Path, short_name: str, type_of_level: str, time_index: int) -> xr.DataArray:
    ds = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": short_name, "typeOfLevel": type_of_level},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    obj = ds[short_name]
    if "time" in obj.dims:
        obj = obj.isel(time=time_index)
    return obj


def _parse_chunk_size(value: str) -> int | None:
    if value.lower() in {"none", "all"}:
        return None
    chunk_size = int(value)
    if chunk_size <= 0:
        raise argparse.ArgumentTypeError("chunk size must be positive or 'none'")
    return chunk_size


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark native spectral regridding chunks.")
    parser.add_argument("--path", type=Path, default=DEFAULT_PATH)
    parser.add_argument("--short-name", default="t")
    parser.add_argument("--type-of-level", default="hybrid")
    parser.add_argument("--time-index", type=int, default=0)
    parser.add_argument("--source-grid", default=None)
    parser.add_argument("--target-grid", default="O480")
    parser.add_argument(
        "--chunk-size",
        type=_parse_chunk_size,
        nargs="+",
        default=[1, 8, 16, 32, 64, None],
    )
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--max-levels", type=int, default=None)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    if args.repeat <= 0:
        parser.error("--repeat must be positive")
    return args


def _print_report(report: dict, *, as_json: bool) -> None:
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return

    print(f"path: {report['path']}")
    print(f"variable: {report['variable']}")
    print(f"source_grid: {report['source_grid']}")
    print(f"target_grid: {report['target_grid']}")
    print(f"input_shape: {report['input_shape']}")
    print(f"open_and_load_s: {report['open_and_load_s']:.6g}")
    print("results:")
    for row in report["results"]:
        chunk = "None" if row["chunk_size"] is None else row["chunk_size"]
        print(
            "  "
            f"chunk_size={chunk} "
            f"mean_s={row['mean_s']:.6g} "
            f"min_s={row['min_s']:.6g} "
            f"max_s={row['max_s']:.6g} "
            f"finite={row['finite']} "
            f"output_shape={row['output_shape']}"
        )


if __name__ == "__main__":
    main()
