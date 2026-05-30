from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import sys
import time
from collections.abc import Callable
from pathlib import Path

import numpy as np
import xarray as xr

from fullpos import horizontal_interpolate


DEFAULT_INPUT = (
    r"L:\ERA5_test\adjustment_test"
    r"\ERA5_Reanalysis_model_level_1978123118_ml1-137_grid_av_with_lnps.grib"
)


def main() -> None:
    args = parse_args()
    path = Path(args.input)
    if not path.exists():
        raise FileNotFoundError(path)

    print_environment(args)
    q = load_q(path, max_levels=args.max_levels)
    print(
        f"input={path}\n"
        f"source_grid={args.source_grid.upper()} target_grid={args.target_grid.upper()} "
        f"shape={tuple(q.shape)} chunks={args.chunks}"
    )

    timings = []
    outputs = {}
    for chunk in args.chunks:
        label = f"q monotonic regrid chunks_hybrid={chunk}"
        timing, out = time_call(
            label,
            lambda chunk=chunk: run_q_regrid(
                q,
                source_grid=args.source_grid,
                target_grid=args.target_grid,
                hybrid_chunk=chunk,
            ),
            repeat=args.repeat,
            warmup=args.warmup,
        )
        timings.append(timing)
        outputs[chunk] = out

    if len(args.chunks) > 1:
        reference = outputs[args.chunks[0]]
        for chunk in args.chunks[1:]:
            print_comparison(
                f"chunks_hybrid={chunk} vs chunks_hybrid={args.chunks[0]}",
                outputs[chunk].values,
                reference.values,
            )
    first_output = outputs[args.chunks[0]]
    external_comparison = None
    if args.compare_output:
        baseline = np.load(args.compare_output)
        external_comparison = compare_arrays(first_output.values, baseline)
        print(
            f"compare_output={args.compare_output}: "
            f"bitwise_equal={external_comparison['bitwise_equal']} "
            f"max_abs={external_comparison['max_abs']:.9g} "
            f"rmse={external_comparison['rmse']:.9g}"
        )
    if args.save_output:
        save_path = Path(args.save_output)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(save_path, np.asarray(first_output.values))
        print(f"saved_output={save_path}")

    print_summary(timings)
    if args.json:
        payload = {
            "environment": environment_info(),
            "input": str(path),
            "source_grid": args.source_grid.upper(),
            "target_grid": args.target_grid.upper(),
            "shape": tuple(int(v) for v in q.shape),
            "chunks": [int(v) for v in args.chunks],
            "timings": [timing for timing in timings],
            "compare_output": external_comparison,
        }
        Path(args.json).write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"wrote_json={args.json}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark FULLPOS FPINT12 monotonic q regridding on real ERA5 data."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--source-grid", default="N320")
    parser.add_argument("--target-grid", default="F480")
    parser.add_argument("--chunks", type=int, nargs="+", default=[137, 64, 32, 16, 8])
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--warmup", type=int, default=0)
    parser.add_argument("--max-levels", type=int, default=0)
    parser.add_argument("--save-output", help="Optional .npy path for saving the first chunk output.")
    parser.add_argument("--compare-output", help="Optional .npy baseline path to compare against the first chunk output.")
    parser.add_argument("--json", help="Optional path for machine-readable timing output.")
    return parser.parse_args()


def print_environment(args: argparse.Namespace) -> None:
    info = environment_info()
    print(
        f"python={info['python']}\n"
        f"executable={info['executable']}\n"
        f"platform={info['platform']}\n"
        f"numpy={info['numpy']}\n"
        f"omp_num_threads={info['omp_num_threads']}\n"
        f"repeat={args.repeat} warmup={args.warmup}"
    )


def environment_info() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", ""),
        "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS", ""),
        "mkl_num_threads": os.environ.get("MKL_NUM_THREADS", ""),
    }


def load_q(path: Path, *, max_levels: int) -> xr.DataArray:
    ds = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "filter_by_keys": {"shortName": "q"},
            "indexpath": "",
            "read_keys": ["gridType", "pl", "J", "K", "M"],
        },
        decode_timedelta=False,
    )
    q = ds["q"]
    if max_levels > 0:
        q = q.isel(hybrid=slice(0, max_levels))
    return q


def run_q_regrid(
    q: xr.DataArray,
    *,
    source_grid: str,
    target_grid: str,
    hybrid_chunk: int,
) -> xr.DataArray:
    return horizontal_interpolate(
        q,
        source_grid=source_grid,
        target_grid=target_grid,
        method="quadratic12_monotonic",
        chunks={"hybrid": int(hybrid_chunk)},
    )


def time_call(
    name: str,
    func: Callable[[], xr.DataArray],
    *,
    repeat: int,
    warmup: int,
) -> tuple[dict[str, object], xr.DataArray]:
    result = None
    for _ in range(max(0, warmup)):
        result = func()
    times = []
    for index in range(repeat):
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        print(f"{name} run={index + 1} elapsed={elapsed:.3f}s")
    assert result is not None
    timing = {
        "name": name,
        "times": times,
        "median": float(statistics.median(times)),
        "mean": float(statistics.mean(times)),
        "min": float(min(times)),
        "max": float(max(times)),
    }
    print(f"{name} median={timing['median']:.3f}s mean={timing['mean']:.3f}s")
    return timing, result


def print_comparison(label: str, candidate: np.ndarray, reference: np.ndarray) -> None:
    stats = compare_arrays(candidate, reference)
    print(
        f"{label}: bitwise_equal={stats['bitwise_equal']} "
        f"max_abs={stats['max_abs']:.9g} "
        f"rmse={stats['rmse']:.9g}"
    )


def compare_arrays(candidate: np.ndarray, reference: np.ndarray) -> dict[str, object]:
    cand = np.asarray(candidate)
    ref = np.asarray(reference)
    if cand.shape != ref.shape:
        raise ValueError(f"shape mismatch: candidate={cand.shape}, reference={ref.shape}")
    diff = cand.astype(np.float64) - ref.astype(np.float64)
    return {
        "bitwise_equal": bool(np.array_equal(cand, ref)),
        "max_abs": float(np.nanmax(np.abs(diff))),
        "rmse": float(np.sqrt(np.nanmean(diff * diff))),
    }


def print_summary(timings: list[dict[str, object]]) -> None:
    print("\nsummary")
    for timing in timings:
        print(
            f"{timing['name']}: median={timing['median']:.3f}s "
            f"mean={timing['mean']:.3f}s min={timing['min']:.3f}s max={timing['max']:.3f}s"
        )
    if len(timings) > 1:
        baseline = float(timings[0]["median"])
        for timing in timings[1:]:
            print(f"{timing['name']}_vs_first_speedup={baseline / float(timing['median']):.3f}x")


if __name__ == "__main__":
    main()
