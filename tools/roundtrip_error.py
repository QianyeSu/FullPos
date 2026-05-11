from __future__ import annotations

import argparse
import json
from pathlib import Path
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
    t0 = perf_counter()
    field = _open_field(args.path, args.short_name, args.type_of_level, args.time_index)
    if args.max_levels is not None and "hybrid" in field.dims:
        field = field.isel(hybrid=slice(0, args.max_levels))
    t1 = perf_counter()

    source_grid = args.source_grid or infer_grid_name_from_attrs(field.attrs)
    forward = Regridder(source_grid, args.target_grid, chunk_size=args.chunk_size)
    backward = Regridder(args.target_grid, source_grid, chunk_size=args.chunk_size)

    out = forward.regrid_data_array(field)
    t2 = perf_counter()
    roundtrip = backward.regrid_data_array(out)
    t3 = perf_counter()

    original = np.asarray(field.values, dtype=np.float64)
    restored = np.asarray(roundtrip.values, dtype=np.float64)
    metrics = _error_metrics(original, restored)
    report = {
        "path": str(args.path),
        "variable": args.short_name,
        "source_grid": source_grid,
        "target_grid": args.target_grid,
        "chunk_size": args.chunk_size,
        "input_shape": list(field.shape),
        "target_shape": list(out.shape),
        "roundtrip_shape": list(roundtrip.shape),
        "open_s": t1 - t0,
        "forward_s": t2 - t1,
        "backward_s": t3 - t2,
        "total_s": t3 - t0,
        **metrics,
    }
    if args.by_level:
        report["by_level"] = _by_level_metrics(field, original, restored, args.level_dim)
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


def _error_metrics(original: np.ndarray, restored: np.ndarray) -> dict[str, float | int | bool]:
    if original.shape != restored.shape:
        raise ValueError(f"roundtrip shape mismatch: {original.shape} != {restored.shape}")
    mask = np.isfinite(original) & np.isfinite(restored)
    diff = restored[mask] - original[mask]
    if diff.size == 0:
        raise ValueError("no finite points available for error analysis")
    rmse = float(np.sqrt(np.mean(diff * diff)))
    source_std = float(np.std(original[mask]))
    return {
        "finite": bool(mask.all()),
        "count": int(diff.size),
        "bias": float(np.mean(diff)),
        "mean_abs": float(np.mean(np.abs(diff))),
        "max_abs": float(np.max(np.abs(diff))),
        "rmse": rmse,
        "source_std": source_std,
        "relative_rmse": float(rmse / source_std) if source_std else float("nan"),
    }


def _by_level_metrics(
    field: xr.DataArray,
    original: np.ndarray,
    restored: np.ndarray,
    level_dim: str,
) -> list[dict[str, float | int | bool]]:
    if level_dim not in field.dims:
        raise ValueError(f"level dimension {level_dim!r} not found in {field.dims}")
    axis = field.get_axis_num(level_dim)
    original_by_level = np.moveaxis(original, axis, 0)
    restored_by_level = np.moveaxis(restored, axis, 0)
    coord = field.coords.get(level_dim)
    rows = []
    for index, (source_level, restored_level) in enumerate(
        zip(original_by_level, restored_by_level)
    ):
        metrics = _error_metrics(
            np.asarray(source_level, dtype=np.float64),
            np.asarray(restored_level, dtype=np.float64),
        )
        level_value = int(coord.values[index]) if coord is not None else index
        rows.append({"index": index, "level": level_value, **metrics})
    return rows


def _parse_chunk_size(value: str) -> int | None:
    if value.lower() in {"none", "all"}:
        return None
    chunk_size = int(value)
    if chunk_size <= 0:
        raise argparse.ArgumentTypeError("chunk size must be positive or 'none'")
    return chunk_size


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run O/N Gaussian roundtrip error analysis.")
    parser.add_argument("--path", type=Path, default=DEFAULT_PATH)
    parser.add_argument("--short-name", default="t")
    parser.add_argument("--type-of-level", default="hybrid")
    parser.add_argument("--time-index", type=int, default=0)
    parser.add_argument("--source-grid", default=None)
    parser.add_argument("--target-grid", default="O480")
    parser.add_argument("--chunk-size", type=_parse_chunk_size, default=64)
    parser.add_argument("--max-levels", type=int, default=None)
    parser.add_argument("--by-level", action="store_true")
    parser.add_argument("--level-dim", default="hybrid")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def _print_report(report: dict, *, as_json: bool) -> None:
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    by_level = report.get("by_level")
    for key, value in report.items():
        if key == "by_level":
            continue
        if isinstance(value, float):
            print(f"{key}: {value:.6g}")
        else:
            print(f"{key}: {value}")
    if by_level:
        print("by_level:")
        for row in by_level:
            print(
                "  "
                f"index={row['index']} "
                f"level={row['level']} "
                f"rmse={row['rmse']:.6g} "
                f"max_abs={row['max_abs']:.6g} "
                f"relative_rmse={row['relative_rmse']:.6g}"
            )


if __name__ == "__main__":
    main()
