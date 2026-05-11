from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

from fullpos._vertical.pressure import _infer_reference_pressure_from_ak, prepare_pressure_request
from fullpos._vertical.validation import skyborn_pressure_reference


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Skyborn's compiled hybrid-to-pressure interpolation as a "
            "development reference for FullPos vertical work."
        )
    )
    parser.add_argument("model", type=Path, help="Path to model-level GRIB/NetCDF input")
    parser.add_argument("surface", type=Path, help="Path to surface-pressure GRIB/NetCDF input")
    parser.add_argument(
        "--coeffs",
        type=Path,
        default=None,
        help=(
            "Optional path to an external hybrid coefficient dataset. "
            "If omitted, GRIB pv metadata from the model field is used."
        ),
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        default=["t"],
        help="Model-level variables to interpolate, for example: t u v q",
    )
    parser.add_argument(
        "--surface-var",
        default="sp",
        help="Surface-pressure variable to open from the surface file",
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
        help="Skyborn interpolation method",
    )
    parser.add_argument(
        "--p0",
        type=float,
        default=None,
        help=(
            "Optional hybrid reference pressure passed to Skyborn. "
            "If omitted, the tool infers 1 for pressure-unit A coefficients "
            "such as ERA5 ap/hyam, otherwise 100000."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional NetCDF output path for the interpolated reference dataset",
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
    variable_names = [str(name) for name in args.variables]
    first = _open_field(args.model, variable_names[0], "hybrid").load()
    surface = _open_field(args.surface, args.surface_var, "surface").load()

    request = prepare_pressure_request(
        first,
        levels=args.levels,
        surface_pressure=surface,
        hybrid_coefficients=args.coeffs,
    )
    p0 = args.p0
    if p0 is None:
        p0 = _infer_reference_pressure_from_ak(request.ak)

    outputs: dict[str, xr.DataArray] = {}
    for name in variable_names:
        field = _open_field(args.model, name, "hybrid").load()
        out = skyborn_pressure_reference(
            field,
            levels=request.levels,
            surface_pressure=request.surface_pressure,
            hybrid_coefficients=args.coeffs,
            method=args.method,
            p0=p0,
        )
        outputs[name] = out
        finite = int(np.isfinite(out.values).sum())
        total = int(out.size)
        print(f"{name}: dims={out.dims} shape={out.shape} finite={finite}/{total}")
        if "plev" in out.dims:
            for index, level in enumerate(np.asarray(out["plev"].values, dtype=np.float64)):
                level_total = int(out.isel(plev=index).size)
                level_finite = int(np.isfinite(out.isel(plev=index).values).sum())
                print(f"{name}: plev={level:g} finite={level_finite}/{level_total}")

    dataset = xr.Dataset(outputs)
    dataset.attrs["reference_backend"] = "skyborn.interp.interpolation.interp_hybrid_to_pressure"
    dataset.attrs["reference_method"] = args.method
    dataset.attrs["reference_levels_pa"] = ",".join(f"{value:g}" for value in request.levels)

    print(f"levels={request.levels.tolist()}")
    print(f"hybrid_dim={request.hybrid_dim}")
    print(f"p0={p0}")
    print(f"surface_pressure_shape={tuple(request.surface_pressure.shape)}")

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        dataset.to_netcdf(args.output)
        print(f"output={args.output}")


if __name__ == "__main__":
    main()
