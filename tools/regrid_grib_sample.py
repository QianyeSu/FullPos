from __future__ import annotations

from pathlib import Path

import xarray as xr

from fullpos import regrid


PATH = Path(r"L:\ERA5_test\era5_reanalysis_model_level_20250102_packing_CCSDS_O320.grib2")


def main() -> None:
    ds = xr.open_dataset(
        PATH,
        engine="cfgrib",
        backend_kwargs={
            "indexpath": "",
            "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
            "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
        },
    )
    field = ds["t"].isel(time=0, hybrid=0)
    out = regrid(field, source_grid="O320", target_grid="O480")
    print(out)
    print("finite output points:", int(out.count()), "/", out.size)
    print("input min/max:", float(field.min()), float(field.max()))
    print("output min/max:", float(out.min()), float(out.max()))


if __name__ == "__main__":
    main()
