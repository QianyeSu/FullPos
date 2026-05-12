# fullpos

Independent Python package for ECMWF/OpenIFS FULLPOS-style regridding and
interpolation.

The interpolation backend is the native OpenIFS/ECTRANS/FULLPOS path. The
package does not silently fall back to Skyborn or another spectral
implementation.

## v0.1.0a0 Preview Boundary

v0.1.0a0 previews native spectral `regrid`, `fit`, `synthesis`, and filtering
paths, plus the selected native vertical pressure and APACHE-related routes
that are already wired in this tree. It is not the complete FULLPOS `POS`
model-level branch, and it does not include the full reduced-grid horizontal
workflow.

Install into the active environment with the project’s normal package
installation flow, then validate with `pytest tests/test_diagnostics.py` and a
quick `fullpos.capabilities()` check. The diagnostics test suite currently
includes a native smoke check, so backend availability must be working for the
runtime report to pass.

## Spectral regridding

```python
import xarray as xr
from fullpos import Regridder, regrid

ds = xr.open_dataset(
    "path/to/model-level-o320.grib2",
    engine="cfgrib",
    backend_kwargs={
        "indexpath": "",
        "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
        "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
    },
)

# source_grid can be inferred from GRIB_gridType, GRIB_N, and GRIB_pl.
out = regrid(ds["t"].isel(time=0), target_grid="O480", chunk_size=64)

regridder = Regridder.from_dataarray(ds["t"], target_grid="O480", chunk_size=64)
out2 = regridder.regrid_data_array(ds["t"].isel(time=0))
```

`chunk_size` controls how many fields are sent to the native ECTRANS batch call
at once. Use a smaller value to reduce peak memory, or `chunk_size=None` to
process the whole batch in one call.

Spectral regridding requires finite global input. Fields decoded from GRIB
bitmaps, for example SST over land, contain missing values and are rejected by
default with `missing_policy="error"`. Full native FULLPOS land/sea-mask
weighted SST interpolation is not implemented yet.

`regrid(..., target_grid="LL1.0")` and `regrid(..., target_grid="LL0.25")`
are supported for xarray inputs. These regular lat/lon targets stay on the
native FULLPOS horizontal path, return cell-center latitude coordinates, and
may first promote a Gaussian source through native spectral regridding when
the stencil needs extra polar coverage. This is not a Python fallback.
Current stable `regrid` support is Gaussian/O/F/N -> Gaussian/O/F/N and
Gaussian/O/F/N -> regular lat/lon (`LL1.0`, `LL0.25`). `LL` is output-only at
this stage and is not supported as `source_grid`.

Dataset inputs can be regridded directly. Variables without horizontal grid
dimensions are skipped by default unless explicitly requested with
`variables=`.

```python
out = regrid(
    ds.isel(time=0),
    target_grid="O480",
    variables=["t", "u", "v"],
    chunk_size=64,
)
```

## Horizontal interpolation

Use `horizontal_interpolate()` for native FULLPOS horizontal interpolation. The
supported methods are:

- `method="bilinear"`: `SUHOW1/SUHOW2` plus `FPINT4`.
- `method="quadratic12"`: `SUHOW1/SUHOW2` plus `FPINT12`.
- `method="nearest"`: `SUHOX1` plus `FPNEAR`.
- `method="average"`: `SUHOX1` plus `FPAVG`.

```python
from fullpos import horizontal_interpolate

out = horizontal_interpolate(
    ds[["t", "u", "v"]],
    target_lats=target_lats,
    target_lons=target_lons,
    method="bilinear",
    variables=["t", "u", "v"],
    chunks={"time": 1, "hybrid": 10},
)
```

Packed reduced Gaussian input is supported when `source_grid="O96"`,
`source_pl=...`, or xarray GRIB metadata contains `GRIB_pl`.

Mask-aware horizontal interpolation with `source_mask` is intentionally still
reported as not implemented until the corresponding native FULLPOS land/sea
mask workflow is wired in.
Mask-aware SST interpolation is still not implemented.

## Vertical pressure interpolation

Use `vertical_interpolate(..., target="pressure")` for native FULLPOS
hybrid-to-pressure interpolation. The current backend calls the OpenIFS/FULLPOS
PP-chain routines `PPINIT`, `PPFLEV`, `PPQ`, `PPUV`, `PPT`, and `PPSTA`.

```python
from fullpos import vertical_interpolate

out = vertical_interpolate(
    ds,
    target="pressure",
    levels=[20000.0, 30000.0, 50000.0],
    variables=["t", "u", "v", "q"],
    chunks={"time": 1, "values": 10000},
    surface_pressure=surface_ds["sp"],
)
```

Temperature uses `PPT`, paired `u`/`v` variables use `PPUV`, and other scalar
variables use `PPQ`. This is not yet the complete `APACHE`/`LESCALE` workflow.

`target="model_level"` is also available as a native PP-chain path when target
hybrid coefficients are supplied with `target_hybrid_coefficients=...`. It uses
Fortran-computed per-column target pressures and the same `PPQ`/`PPUV`/`PPT`
kernels; it is not the complete FULLPOS `POS` model-level copy branch.

`target="potential_temperature"` is available for theta surfaces:

```python
out = vertical_interpolate(
    ds,
    target="potential_temperature",
    levels=[300.0, 320.0],
    variables=["t", "u", "v"],
    chunks={"time": 1, "values": 10000},
    surface_pressure=surface_ds["sp"],
)
```

The target pressure calculation uses native FULLPOS `GPTET` and `PPLTETA`, then
interpolates with the same native `PPQ`/`PPUV`/`PPT` kernels. The current theta
path uses a dry-air kappa constant and does not yet include the full moist
`GPRCP` context.

`target="temperature"` is available for isothermal surfaces:

```python
out = vertical_interpolate(
    ds,
    target="temperature",
    levels=[250.0, 270.0],
    variables=["t", "u", "v", "q"],
    chunks={"time": 1, "values": 10000},
    surface_pressure=surface_ds["sp"],
)
```

The current vendored source set does not include a usable `PPLTEMP`; the native
temperature path uses available FULLPOS `PPLTW`/`FPPS` to compute target
pressures before `PPQ`/`PPUV`/`PPT` interpolation.

## Spectral filtering

```python
from fullpos import spectral_filter

smooth = spectral_filter(ds["t"].isel(time=0), grid="O320", ntrunc=159, chunk_size=64)
```

The current filter API applies filter profiles to native ECTRANS spectral
coefficients. Gaussian, low-pass, generic diagonal coefficient filters, and
filter matrix read/write helpers are available.

## GRIB output

Use `to_grib()` to write xarray output back to GRIB with an existing GRIB
template:

```python
from fullpos import regrid, to_grib

out = regrid(ds["t"].isel(time=0), target_grid="O480")
to_grib(out, "t_o480.grib2", template="template_o480_t.grib2")
```

The writer uses ecCodes to clone matching template messages and replace only
the field values. This preserves template grid, level, time, packing, bitmap,
and product metadata. It also means the template must already match the desired
output layout; `to_grib()` does not convert an `O320` template into `O480`
metadata and does not infer complete GRIB messages from xarray metadata.

Install the optional GRIB dependencies with `pip install -e .[grib]
--no-build-isolation`.

## Validation tools

Run the test suite:

```powershell
python -m pytest -q
```

Run a roundtrip error check:

```powershell
python tools/roundtrip_error.py --max-levels 4 --chunk-size 4
```

Benchmark batch sizes after loading the GRIB field once:

```powershell
python tools/benchmark_chunks.py --chunk-size 16 64 none
```

Check the saved validation report against regression thresholds:

```powershell
python tools/check_regression.py
```

Inspect runtime backend status:

```python
from fullpos import backend_info

print(backend_info())
```

Run a native backend health check:

```powershell
python tools/doctor.py
```

## Native prefix

By default, builds use:

```text
extern/fullpos/local
```

Use a different FIAT/ECTRANS installation prefix at build time with Meson:

```powershell
python -m pip install -e . --no-build-isolation --no-deps --config-settings=setup-args="-Dfullpos_native_prefix=<path-to-native-prefix>"
```

For runtime diagnostics or local testing, `FULLPOS_NATIVE_PREFIX` can override
the configured prefix.

Wheel builds install the FIAT/ECTRANS shared libraries from the selected prefix
into `fullpos/_native` by default. This keeps installed wheels from depending on
the source-tree `extern/fullpos/local` path. The Meson option
`-Dfullpos_bundle_native_runtime=false` disables that copy for local external
dependency builds.

## Native libraries

On Windows, `.dll` files are dynamic libraries. The
`_ectrans.cp312-win_amd64.pyd` file is also a dynamic library, but it is a
Python extension module that Python imports directly. The `.pyd` depends on
FIAT/ECTRANS `.dll` files, and those DLLs can further depend on external
runtime DLLs such as OpenBLAS and gfortran.

On Linux, the equivalent native libraries normally use `.so`. On macOS, they
normally use `.dylib`.

For wheel builds, FIAT/ECTRANS are bundled under `fullpos/_native`, while
platform repair tools (`auditwheel`, `delocate`, `delvewheel`) are expected to
vendor external runtime dependencies such as OpenBLAS and gfortran. Source-tree
editable installs can still use external dependency mode through
`FULLPOS_NATIVE_PREFIX` and the current shell `PATH`/loader search path.
