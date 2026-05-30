from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from fullpos.grids import parse_grid
from fullpos.native import add_native_runtime_dir
from fullpos.spectral import spectral_wind_synthesis
from skyborn.spharm import Spharmt


DEFAULT_INPUT = (
    r"L:\ERA5_test\adjustment_test"
    r"\ERA5_Reanalysis_model_level_1978123118_ml1-137_grid_av_with_lnps.grib"
)
DEFAULT_TARGET_GRID = "F480"
EARTH_RADIUS = 6_371_229.0


@dataclass(frozen=True)
class TimingResult:
    name: str
    times: list[float]
    median: float
    mean: float
    min: float
    max: float


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    print_environment(args)
    data = load_spectral_fields(input_path, max_levels=args.max_levels)
    target_grid = parse_grid(args.target_grid)
    if target_grid.is_reduced:
        raise ValueError("This benchmark targets regular Gaussian output, for example F480")
    target_pl = np.full(target_grid.nlat, target_grid.work_nlon, dtype=np.int32)
    ntrunc = infer_ntrunc(data["t"].shape[-1])
    print(
        f"input={input_path}\n"
        f"target_grid={args.target_grid.upper()} nlat={target_grid.nlat} "
        f"nlon={target_grid.work_nlon} npoints={target_grid.size}\n"
        f"levels={data['t'].shape[0]} ntrunc=T{ntrunc}"
    )

    add_native_runtime_dir()
    from fullpos import _ectrans

    flat_tvd = np.ascontiguousarray(
        np.concatenate([data["t"], data["vo"], data["d"]], axis=0),
        dtype=np.float64,
    )
    flat_vd = np.ascontiguousarray(
        np.concatenate([data["vo"], data["d"]], axis=0),
        dtype=np.float64,
    )

    # Warm the import path and the native transform plan.
    _ectrans.synthesis(flat_tvd[:1], target_pl, int(ntrunc))

    stages = set(args.stages)
    run_synthesis = "all" in stages or "synthesis" in stages
    run_wind = "all" in stages or "wind" in stages
    run_vordiv = "all" in stages or "vordiv" in stages

    baseline_synth = None
    synth = None
    if args.compare_old_synthesis and run_synthesis:
        with temporary_env(
            {
                "FULLPOS_ECTRANS_DISABLE_CACHE": "1",
                "FULLPOS_ECTRANS_DISABLE_GLOBAL_FASTPATH": "1",
            }
        ):
            _old_timing, baseline_synth = time_call(
                "spectral t/vo/d synthesis old C wrapper path",
                lambda: np.asarray(
                    _ectrans.synthesis(flat_tvd, target_pl, int(ntrunc)),
                    dtype=np.float64,
                ),
                repeat=args.old_repeat,
                warmup=0,
            )
    else:
        _old_timing = None

    if run_synthesis or run_wind:
        current_timing, synth = time_call(
            "spectral t/vo/d synthesis current",
            lambda: np.asarray(
                _ectrans.synthesis(flat_tvd, target_pl, int(ntrunc)),
                dtype=np.float64,
            ),
            repeat=args.repeat,
            warmup=args.warmup if run_synthesis else 0,
        )
        if baseline_synth is not None:
            print_comparison("current vs old synthesis", synth, baseline_synth)
    else:
        current_timing = None

    timings = []
    if _old_timing is not None:
        timings.append(_old_timing)
    if current_timing is not None and run_synthesis:
        timings.append(current_timing)

    baseline_wind = None
    if run_wind:
        assert synth is not None
        sph = Spharmt(
            target_grid.work_nlon,
            target_grid.nlat,
            rsphere=EARTH_RADIUS,
            gridtype="gaussian",
            legfunc=args.legfunc,
            precision=args.precision,
        )
        vo_grid, div_grid = split_vordiv_grid(synth, data["vo"].shape[0], target_grid)

        baseline_wind_timing, baseline_wind = time_call(
            "wind synthesis current script path",
            lambda: wind_current_script_path(sph, vo_grid, div_grid, ntrunc),
            repeat=args.repeat,
            warmup=args.warmup,
        )
        batched_wind_timing, batched_wind = time_call(
            "wind synthesis batched grdtospec candidate",
            lambda: wind_batched_grdtospec_path(sph, vo_grid, div_grid, ntrunc),
            repeat=args.repeat,
            warmup=args.warmup,
        )
        print_pair_comparison("batched grdtospec vs current wind", batched_wind, baseline_wind)
        timings.extend([baseline_wind_timing, batched_wind_timing])

    if args.include_experimental_vordiv and run_vordiv:
        experimental_timing, experimental_wind = time_call(
            "wind synthesis ECTRANS vordiv experimental",
            lambda: wind_ectrans_vordiv_path(data["vo"], data["d"], ntrunc, target_grid),
            repeat=args.repeat,
            warmup=args.warmup,
        )
        if baseline_wind is not None:
            print_pair_comparison("ECTRANS vordiv vs current wind", experimental_wind, baseline_wind)
        timings.append(experimental_timing)

    print_summary(timings)
    if args.json:
        payload = {
            "environment": environment_info(),
            "input": str(input_path),
            "target_grid": args.target_grid.upper(),
            "levels": int(data["t"].shape[0]),
            "ntrunc": int(ntrunc),
            "timings": [timing.__dict__ for timing in timings],
        }
        Path(args.json).write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"wrote_json={args.json}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark the ERA5 mass-adjustment spectral t/vo/d synthesis and "
            "wind synthesis hotspots without changing FULLPOS/ECTRANS Fortran code."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--target-grid", default=DEFAULT_TARGET_GRID)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--old-repeat", type=int, default=1)
    parser.add_argument("--max-levels", type=int, default=0)
    parser.add_argument("--precision", choices=("single", "double", "auto"), default="single")
    parser.add_argument("--legfunc", choices=("stored", "computed"), default="stored")
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=("all", "synthesis", "wind", "vordiv"),
        default=["all"],
        help=(
            "Benchmark only selected stages. Use this for full 137-level runs "
            "when the full combination would exceed local runtime limits."
        ),
    )
    parser.add_argument(
        "--compare-old-synthesis",
        action="store_true",
        help="Also time the pre-cache/pre-fast-path C wrapper behavior using env escape hatches.",
    )
    parser.add_argument(
        "--include-experimental-vordiv",
        action="store_true",
        help=(
            "Also time native ECTRANS vorticity/divergence to U/V synthesis. "
            "This is fast but is reported separately because it may not be bitwise-identical "
            "to the current skyborn wind path."
        ),
    )
    parser.add_argument("--json", help="Optional path for machine-readable timing output.")
    return parser.parse_args()


def print_environment(args: argparse.Namespace) -> None:
    info = environment_info()
    print(
        f"python={info['python']}\n"
        f"executable={info['executable']}\n"
        f"platform={info['platform']}\n"
        f"numpy={info['numpy']}\n"
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


def load_spectral_fields(path: Path, *, max_levels: int) -> dict[str, np.ndarray]:
    dataset = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "filter_by_keys": {"shortName": ["t", "vo", "d"]},
            "indexpath": "",
            "read_keys": ["gridType", "J", "K", "M"],
        },
        decode_timedelta=False,
    )
    out: dict[str, np.ndarray] = {}
    for name in ("t", "vo", "d"):
        raw = np.asarray(dataset[name].transpose("hybrid", "values").values, dtype=np.float64)
        if max_levels > 0:
            raw = raw[:max_levels]
        out[name] = np.ascontiguousarray(raw, dtype=np.float64)
    return out


def infer_ntrunc(nspec2: int) -> int:
    ntrunc = int((np.sqrt(1 + 4 * int(nspec2)) - 3) / 2)
    if (ntrunc + 1) * (ntrunc + 2) != int(nspec2):
        raise ValueError(f"spectral coefficient count is not triangular: {nspec2}")
    return ntrunc


def time_call(
    name: str,
    func: Callable[[], object],
    *,
    repeat: int,
    warmup: int,
) -> tuple[TimingResult, object]:
    if repeat <= 0:
        raise ValueError("repeat must be positive")
    result = None
    for _ in range(max(0, warmup)):
        result = func()
    times: list[float] = []
    for index in range(repeat):
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        print(f"{name} run={index + 1} elapsed={elapsed:.3f}s")
    assert result is not None
    timing = TimingResult(
        name=name,
        times=times,
        median=float(statistics.median(times)),
        mean=float(statistics.mean(times)),
        min=float(min(times)),
        max=float(max(times)),
    )
    print(f"{name} median={timing.median:.3f}s mean={timing.mean:.3f}s")
    return timing, result


def split_vordiv_grid(
    synth: np.ndarray,
    nlev: int,
    target_grid,
) -> tuple[np.ndarray, np.ndarray]:
    arr = np.asarray(synth, dtype=np.float64)
    vo = arr[nlev : 2 * nlev].astype(np.float32, copy=False)
    div = arr[2 * nlev : 3 * nlev].astype(np.float32, copy=False)
    vo_grid = vo.reshape((nlev, target_grid.nlat, target_grid.work_nlon)).transpose(1, 2, 0)
    div_grid = div.reshape((nlev, target_grid.nlat, target_grid.work_nlon)).transpose(1, 2, 0)
    return (
        np.ascontiguousarray(vo_grid[:, :, None, :], dtype=np.float32),
        np.ascontiguousarray(div_grid[:, :, None, :], dtype=np.float32),
    )


def wind_current_script_path(
    sph: Spharmt,
    vo_grid: np.ndarray,
    div_grid: np.ndarray,
    ntrunc: int,
) -> tuple[np.ndarray, np.ndarray]:
    vort_spec = sph.grdtospec(np.ascontiguousarray(vo_grid, dtype=np.float32), ntrunc=ntrunc)
    div_spec = sph.grdtospec(np.ascontiguousarray(div_grid, dtype=np.float32), ntrunc=ntrunc)
    u, v = sph.getuv(vort_spec, div_spec)
    return np.ascontiguousarray(u, dtype=np.float32), np.ascontiguousarray(v, dtype=np.float32)


def wind_batched_grdtospec_path(
    sph: Spharmt,
    vo_grid: np.ndarray,
    div_grid: np.ndarray,
    ntrunc: int,
) -> tuple[np.ndarray, np.ndarray]:
    nlev = vo_grid.shape[-1]
    combined = np.ascontiguousarray(np.concatenate([vo_grid, div_grid], axis=-1), dtype=np.float32)
    spec = sph.grdtospec(combined, ntrunc=ntrunc)
    u, v = sph.getuv(spec[..., :nlev], spec[..., nlev:])
    return np.ascontiguousarray(u, dtype=np.float32), np.ascontiguousarray(v, dtype=np.float32)


def wind_ectrans_vordiv_path(
    raw_vo: np.ndarray,
    raw_d: np.ndarray,
    ntrunc: int,
    target_grid,
) -> tuple[np.ndarray, np.ndarray]:
    nlev = raw_vo.shape[0]
    u, v = spectral_wind_synthesis(
        np.ascontiguousarray(raw_vo, dtype=np.float64),
        np.ascontiguousarray(raw_d, dtype=np.float64),
        grid=target_grid,
        ntrunc=int(ntrunc),
    )
    u = np.asarray(u, dtype=np.float64).reshape((nlev, target_grid.nlat, target_grid.work_nlon)).transpose(1, 2, 0)
    v = np.asarray(v, dtype=np.float64).reshape((nlev, target_grid.nlat, target_grid.work_nlon)).transpose(1, 2, 0)
    return (
        np.ascontiguousarray(u[:, :, None, :], dtype=np.float32),
        np.ascontiguousarray(v[:, :, None, :], dtype=np.float32),
    )


def print_pair_comparison(
    label: str,
    candidate: tuple[np.ndarray, np.ndarray],
    baseline: tuple[np.ndarray, np.ndarray],
) -> None:
    print_comparison(label + " u", candidate[0], baseline[0])
    print_comparison(label + " v", candidate[1], baseline[1])


def print_comparison(label: str, candidate: np.ndarray, baseline: np.ndarray) -> None:
    cand = np.asarray(candidate)
    base = np.asarray(baseline)
    diff = cand.astype(np.float64) - base.astype(np.float64)
    max_abs = float(np.nanmax(np.abs(diff)))
    rmse = float(np.sqrt(np.nanmean(diff * diff)))
    mean_abs = float(np.nanmean(np.abs(diff)))
    print(
        f"{label}: bitwise_equal={bool(np.array_equal(cand, base))} "
        f"max_abs={max_abs:.9g} rmse={rmse:.9g} mean_abs={mean_abs:.9g}"
    )


def print_summary(timings: list[TimingResult]) -> None:
    print("\nsummary")
    for timing in timings:
        print(
            f"{timing.name}: median={timing.median:.3f}s "
            f"mean={timing.mean:.3f}s min={timing.min:.3f}s max={timing.max:.3f}s"
        )
    lookup = {timing.name: timing for timing in timings}
    old = lookup.get("spectral t/vo/d synthesis old C wrapper path")
    current = lookup.get("spectral t/vo/d synthesis current")
    if old is not None and current is not None:
        print(f"spectral_current_vs_old_speedup={old.median / current.median:.3f}x")
    baseline = lookup.get("wind synthesis current script path")
    batched = lookup.get("wind synthesis batched grdtospec candidate")
    if baseline is not None and batched is not None:
        print(f"wind_batched_vs_current_speedup={baseline.median / batched.median:.3f}x")
    experimental = lookup.get("wind synthesis ECTRANS vordiv experimental")
    if baseline is not None and experimental is not None:
        print(f"wind_ectrans_vordiv_vs_current_speedup={baseline.median / experimental.median:.3f}x")


class temporary_env:
    def __init__(self, values: dict[str, str]) -> None:
        self.values = values
        self.old: dict[str, str | None] = {}

    def __enter__(self) -> None:
        for key, value in self.values.items():
            self.old[key] = os.environ.get(key)
            os.environ[key] = value

    def __exit__(self, exc_type, exc, tb) -> None:
        for key, value in self.old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


if __name__ == "__main__":
    main()
