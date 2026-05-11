from __future__ import annotations

import argparse
from pathlib import Path

import xarray as xr

from fullpos._vertical.pressure import prepare_pressure_request


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate ERA5/OpenIFS pressure-level vertical interpolation inputs."
    )
    parser.add_argument("model", type=Path, help="Path to model-level GRIB/NetCDF input")
    parser.add_argument("surface", type=Path, help="Path to surface-pressure GRIB/NetCDF input")
    parser.add_argument(
        "--coeffs",
        type=Path,
        default=None,
        help="Optional path to an external hybrid A/B coefficient dataset",
    )
    parser.add_argument(
        "--model-var",
        default="t",
        help="Model-level variable to open from the model file",
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
    return parser.parse_args()


def _open_field(path: Path, short_name: str, type_of_level: str):
    if path.suffix.lower() in {".nc", ".nc4"}:
        return xr.open_dataset(path)[short_name]
    return xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": short_name, "typeOfLevel": type_of_level},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "pv", "packingType"],
        },
    )[short_name]


def main() -> None:
    args = parse_args()
    model = _open_field(args.model, args.model_var, "hybrid").load()
    surface = _open_field(args.surface, args.surface_var, "surface").load()

    request = prepare_pressure_request(
        model,
        levels=args.levels,
        surface_pressure=surface,
        hybrid_coefficients=args.coeffs,
    )

    print("request_ok=True")
    print(f"hybrid_dim={request.hybrid_dim}")
    print(f"levels={request.levels.tolist()}")
    print(f"ak_len={request.ak.size}")
    print(f"bk_len={request.bk.size}")
    print(f"surface_pressure_shape={tuple(request.surface_pressure.shape)}")
    if "time" in request.surface_pressure.coords:
        print(f"time_start={request.surface_pressure.time.values[0]}")
        print(f"time_end={request.surface_pressure.time.values[-1]}")


if __name__ == "__main__":
    main()
