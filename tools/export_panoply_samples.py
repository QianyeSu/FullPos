from __future__ import annotations

from pathlib import Path
from time import perf_counter

import numpy as np
import xarray as xr

from fullpos import regrid


MODEL_PATH = Path(
    r"L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2"
)
SURFACE_PATH = Path(
    r"L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19950710_hourly_O96.grib"
)
OUT_DIR = Path(r"L:\ERA5_test")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    export_model_sample()
    export_surface_sample()


def export_model_sample() -> None:
    output = OUT_DIR / "fullpos_o96_to_n160_model_tuvq_time0_ml1-10.nc"
    variables = ["t", "u", "v", "q"]
    arrays = []
    start = perf_counter()
    for variable in variables:
        da = _open_variable(MODEL_PATH, variable, "hybrid")
        da = da.isel(time=0, hybrid=slice(0, 10))
        out = regrid(da, target_grid="N160", chunk_size=16)
        out.name = variable
        arrays.append(out)
        print(
            f"model {variable}: {da.shape} -> {out.shape}, "
            f"finite={bool(np.isfinite(out.values).all())}"
        )
    ds = xr.merge(arrays, compat="override")
    ds.attrs.update(
        {
            "title": "FullPos O96 to N160 model-level sample",
            "source_file": str(MODEL_PATH),
            "target_grid": "N160",
            "note": "First time step, first 10 hybrid model levels.",
        }
    )
    ds.to_netcdf(output)
    print(f"wrote {output} in {perf_counter() - start:.3f}s")


def export_surface_sample() -> None:
    output = OUT_DIR / "fullpos_o96_to_n160_surface_sst_sp_tcwv_msl_time0.nc"
    variables = ["sst", "sp", "tcwv", "msl"]
    arrays = []
    original_arrays = []
    start = perf_counter()
    for variable in variables:
        da = _open_variable(SURFACE_PATH, variable, "surface")
        if "time" in da.dims:
            da = da.isel(time=0)
        original = da.copy(deep=False)
        original.name = f"{variable}_o96_packed"
        original_arrays.append(original)
        input_finite = int(np.isfinite(da.values).sum())
        if input_finite != da.size:
            out = regrid(da, target_grid="N160", method="masked", chunk_size=16)
            out.name = variable
            arrays.append(out)
            print(
                f"surface {variable}: {da.shape} -> {out.shape}, method=masked, "
                f"input_finite={input_finite}/{da.size}, "
                f"output_finite={int(np.isfinite(out.values).sum())}/{out.size}"
            )
            continue
        out = regrid(da, target_grid="N160", chunk_size=16)
        out.name = variable
        arrays.append(out)
        print(
            f"surface {variable}: {da.shape} -> {out.shape}, "
            f"input_finite={input_finite}/{da.size}, "
            f"output_finite={int(np.isfinite(out.values).sum())}/{out.size}"
        )
    ds = xr.merge(arrays + original_arrays, compat="override")
    ds.attrs.update(
        {
            "title": "FullPos O96 to N160 surface sample",
            "source_file": str(SURFACE_PATH),
            "target_grid": "N160",
            "note": (
                "First time step. Variables named *_o96_packed are original packed O96 fields. "
                "Fields with missing values, for example SST with a GRIB bitmap, are "
                "regridded with missing-aware masked interpolation; finite fields use spectral interpolation."
            ),
        }
    )
    ds.to_netcdf(output)
    print(f"wrote {output} in {perf_counter() - start:.3f}s")


def _open_variable(path: Path, short_name: str, type_of_level: str) -> xr.DataArray:
    ds = xr.open_dataset(
        path,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": short_name, "typeOfLevel": type_of_level},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    if short_name not in ds:
        raise KeyError(f"{short_name!r} not found in {path}")
    return ds[short_name]


if __name__ == "__main__":
    main()
