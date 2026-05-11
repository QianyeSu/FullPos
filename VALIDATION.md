# FullPos Validation Notes

This file records local validation results for the required native ECTRANS spectral
interpolation path. Commands assume the `skyborn_dev` environment.

## Environment

```powershell
$env:PATH='F:\PythonPackage\FullPos\extern\fullpos\local\bin;F:\Anaconda3\envs\skyborn_dev;F:\Anaconda3\envs\skyborn_dev\Library\mingw-w64\bin;F:\Anaconda3\envs\skyborn_dev\Library\usr\bin;F:\Anaconda3\envs\skyborn_dev\Library\bin;F:\Anaconda3\envs\skyborn_dev\Scripts;'+$env:PATH
F:\Anaconda3\envs\skyborn_dev\python.exe -m pytest -q
```

Latest result:

```text
35 passed
```

Runtime backend diagnostics:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe -c "from fullpos import backend_info; print(backend_info())"
F:\Anaconda3\envs\skyborn_dev\python.exe tools\doctor.py
```

Expected local status:

```text
backend: native
native_available: True
native_module: fullpos._ectrans
native_runtime_platform: win32
native_runtime_dir_exists: True
required_native_dlls_present: True
native smoke regrid N4 -> N4: PASS
```

There is intentionally no Skyborn or regular-grid fallback in the interpolation path. If the native
FULLPOS/ECTRANS backend is unavailable or fails, the package should raise an error.

The configured native prefix is generated at build time from Meson option
`fullpos_native_prefix`. If the option is not set, it defaults to:

```text
extern/fullpos/local
```

Runtime diagnostics can temporarily override it with:

```powershell
$env:FULLPOS_NATIVE_PREFIX='F:\PythonPackage\FullPos\extern\fullpos\local'
```

External dependency mode:

```text
FullPos does not bundle OpenBLAS, gfortran, OpenMP, C++ runtime, or pthread runtime yet.
doctor checks whether these libraries can be found from the current environment.
```

Current Windows external runtime libraries:

```text
openblas.dll
libgfortran-5.dll
libgomp-1.dll
libstdc++-6.dll
libwinpthread-1.dll
libgcc_s_seh-1.dll
```

## Real Data

Input file:

```text
L:\ERA5_test\era5_reanalysis_model_level_20250102_packing_CCSDS_O320.grib2
```

Variable:

```text
shortName=t, typeOfLevel=hybrid, time=0
```

Detected source grid:

```text
O320, reduced_gg, 137 hybrid levels, 421120 packed points per field
```

## O320 To O480 Benchmark

Command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\benchmark_chunks.py --chunk-size 16 64 none
```

Result:

```text
input_shape: [137, 421120]
open_and_load_s: 7.07125
chunk_size=16: 13.83 s
chunk_size=64: 12.5448 s
chunk_size=None: 11.1412 s
output_shape: [137, 938880]
finite: True
```

Native-only rerun after removing fallback:

```text
input_shape: [137, 421120]
open_and_load_s: 13.5122
chunk_size=64: 17.9093 s
output_shape: [137, 938880]
finite: True
```

Interpretation:

- `chunk_size=None` is fastest on this machine but has the highest peak memory pressure.
- `chunk_size=64` is the current default because it keeps memory bounded with only a modest speed cost.
- Local I/O and runtime have noticeable run-to-run variance; compare multiple repeats before treating small speed changes as meaningful.

## Dataset Path

Command shape:

```python
out = regrid(ds.isel(time=0), target_grid="O480", chunk_size=64)
```

Latest local result:

```text
data_vars: ['t']
dims: {'hybrid': 137, 'values': 938880}
t dims: ('hybrid', 'values')
t shape: (137, 938880)
finite: True
attrs: GRIB_gridType=reduced_gg, GRIB_N=480, GRIB_numberOfPoints=938880
elapsed_s: 13.896
```

## O96 Multi-Variable GRIB Checks

Model-level file:

```text
L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2
```

Command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\regrid_grib_dataset.py "L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2" --variables t u v q --type-of-level hybrid --target-grid O160 --chunk-size 16 --max-levels 10 --strict
```

Result summary:

```text
t: [10, 40320] -> [10, 108160], finite=True
u: [10, 40320] -> [10, 108160], finite=True
v: [10, 40320] -> [10, 108160], finite=True
q: [10, 40320] -> [10, 108160], finite=True
```

Surface file:

```text
L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19950710_hourly_O96.grib
```

Command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\regrid_grib_dataset.py "L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19950710_hourly_O96.grib" --variables sst sp tcwv msl --type-of-level surface --target-grid O160 --chunk-size 16 --strict
```

Spectral result summary:

```text
sst: FAIL ValueError: spectral regridding requires finite input; found 11571 non-finite value(s) out of 40320. Fields with GRIB bitmap or missing values, for example SST over land, need masked surface interpolation instead of global spectral interpolation.
sp: input_finite=40320/40320, output_finite=108160/108160
tcwv: input_finite=40320/40320, output_finite=108160/108160
msl: input_finite=40320/40320, output_finite=108160/108160
```

Interpretation:

- O96 model-level `t/u/v/q` works for finite fields.
- O96 surface `sp/tcwv/msl` works.
- `sst` has `bitmapPresent=1` and `numberOfMissing=11571`; it is now rejected before ECTRANS instead of silently producing all-NaN output.
- SST-like fields need masked surface interpolation, not global spectral interpolation.

Masked SST command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\regrid_grib_dataset.py "L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19950710_hourly_O96.grib" --variables sst --type-of-level surface --target-grid N160 --method masked --chunk-size 16
```

Masked SST result:

```text
sst: [40320] -> [320, 640] source=reduced_gg N96 target=regular_gg N160 input_finite=28749/40320 output_finite=144868/204800 open_s=2.439 regrid_s=0.121 min=271.455 max=305.844
```

Panoply export:

```text
L:\ERA5_test\fullpos_o96_to_n160_surface_sst_sp_tcwv_msl_time0.nc
```

The exported `sst` variable uses `method="masked"` and keeps NaN where no valid source neighbors exist.

OpenIFS reference points:

## O96 Skyborn Vertical Reference Against Official Pressure Levels

Model-level file:

```text
L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2
```

Surface-pressure file:

```text
L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19781201_hourly_O96.grib
```

Pressure-level truth file:

```text
L:\ERA5_Complete\Reanalysis\pressure_level\ERA5_Reanalysis_19781201_6hourly_O96.grib2
```

Command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\compare_pressure_reference.py "L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2" "L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19781201_hourly_O96.grib" "L:\ERA5_Complete\Reanalysis\pressure_level\ERA5_Reanalysis_19781201_6hourly_O96.grib2" --variables t u v q --levels-hpa 300 400 500 --method log
```

Key result:

```text
reference backend: skyborn.interp.interpolation.interp_hybrid_to_pressure
p0: 1.0
grid: O96 reduced_gg packed values
times: 4
levels: 300 / 400 / 500 hPa
```

`log` interpolation summary:

```text
t overall: rmse=0.216124 mae=0.0779591 max_abs=7.1803
t 300 hPa: rmse=0.132792 mae=0.0478971 max_abs=3.69124
t 400 hPa: rmse=0.22049  mae=0.0829634 max_abs=5.21048
t 500 hPa: rmse=0.271806 mae=0.103017  max_abs=7.1803

u overall: rmse=0.153601 mae=0.0543891 max_abs=7.68178
u 300 hPa: rmse=0.0937532 mae=0.037921  max_abs=3.20418
u 400 hPa: rmse=0.147285  mae=0.0551373 max_abs=7.39463
u 500 hPa: rmse=0.200741  mae=0.0701092 max_abs=7.68178

v overall: rmse=0.113571 mae=0.0459642 max_abs=3.90626
v 300 hPa: rmse=0.0812846 mae=0.0344817 max_abs=2.43432
v 400 hPa: rmse=0.11154   mae=0.0462958 max_abs=3.30546
v 500 hPa: rmse=0.140168  mae=0.0571151 max_abs=3.90626

q overall: rmse=4.63937e-06 mae=5.50950e-07 max_abs=5.85705e-04
q 300 hPa: rmse=4.88502e-07 mae=8.05731e-08 max_abs=5.08147e-05
q 400 hPa: rmse=2.43887e-06 mae=3.64340e-07 max_abs=2.25896e-04
q 500 hPa: rmse=7.64097e-06 mae=1.20794e-06 max_abs=5.85705e-04
```

Comparison note:

- For this case, `log` is slightly better than `linear` for `u` and `v`, almost identical for `t`, and mixed but still very close for `q`.
- Errors increase downward from `300 hPa` to `500 hPa`, which is plausible because gradients are usually stronger lower in the column.
- This comparison is against the official ERA5 pressure-level product, not against FULLPOS native vertical output.

## O96 Native FULLPOS Pressure Levels Against Official Pressure Levels

Command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\compare_native_pressure_reference.py "L:\ERA5_Complete\Reanalysis\model_level\ERA5_Reanalysis_19781201_6hourly_ml1-137_O96.grib2" "L:\ERA5_Complete\Reanalysis\surface\ERA5_Reanalysis_surface_19781201_hourly_O96.grib" "L:\ERA5_Complete\Reanalysis\pressure_level\ERA5_Reanalysis_19781201_6hourly_O96.grib2" --variables t u v q --levels-hpa 200 300 400 500 --chunks time=1 values=10000 --output-json "L:\ERA5_test\fullpos_native_pressure_reference_metrics_19781201_o96_200_300_400_500.json"
```

Output files:

```text
L:\ERA5_test\fullpos_native_pressure_reference_metrics_19781201_o96_200_300_400_500.json
L:\ERA5_test\fullpos_native_pressure_o96_tuvq_200_300_400_500.nc
```

Key result:

```text
backend: FULLPOS PP-chain
grid: O96 reduced_gg packed values
times: 4
levels: 200 / 300 / 400 / 500 hPa
t overall: rmse=0.188270 mae=0.0630102 max_abs=7.17754
u overall: rmse=0.134141 mae=0.0447220 max_abs=7.68178
v overall: rmse=0.0999003 mae=0.0385203 max_abs=3.90626
q overall: rmse=4.00771e-06 mae=4.39510e-07 max_abs=5.84819e-04
```

Per-level highlights:

```text
t 200 hPa: rmse=0.0405560 max_abs=1.35688
u 200 hPa: rmse=0.0345799 max_abs=1.36310
v 200 hPa: rmse=0.0349981 max_abs=1.08247
q 200 hPa: rmse=5.35056e-08 max_abs=2.61334e-06
```

Implementation note:

- This native path calls the original OpenIFS/FULLPOS `PPINIT`, `PPFLEV`, `PPQ`, `PPUV`, `PPT`, and `PPSTA` routines.
- It is a native PP-chain wrapper, not the complete `APACHE`/`LESCALE` workflow.

## O96 Native FULLPOS Model-Level PP-Chain Smoke

Output file:

```text
L:\ERA5_test\fullpos_native_model_level_o96_tuvq_time0_ml1-10.nc
```

Key result:

```text
target: model_level
backend: FULLPOS PP-chain column-pressure wrapper
shape: t/u/v/q(model_level=10, values=40320)
max_abs(out - input):
  q: 0.0
  t: 0.5664007958922639
  u: 4.017446673275224
  v: 2.701768687905453
```

Interpretation:

- The column-pressure backend runs on real O96 data and writes a Panoply-readable NetCDF sample.
- This is not the complete FULLPOS `POS` `CDCONF='M'` model-level copy branch. `PPT` and `PPUV` apply variable-specific FULLPOS logic, so identical source/target hybrid coefficients are not guaranteed to be identity for `t/u/v`.

- `F:\openifs-main\ifs-source\arpifs\module\yomfpc.F90`: `NFPMASK` selects no mask, land-sea mask, or separate land/sea masks for surface-field interpolation.
- `F:\openifs-main\ifs-source\arpifs\module\yomwfpb.F90`: `WLAN*` and `WSEA*` store land/sea masked interpolation weights for surface physical fields.
- `F:\openifs-main\ifs-source\arpifs\interpol\fpint4x.F90`: `FPINT4X` is explicitly for fields with missing values.
- `F:\openifs-main\ifs-source\arpifs\climate\updclie.F90`: SST update code tests missing SST values before use.

## Roundtrip Error

Small sample command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\roundtrip_error.py --max-levels 4 --chunk-size 4
```

Small sample result:

```text
input_shape: [4, 421120]
target_shape: [4, 938880]
roundtrip_shape: [4, 421120]
forward_s: 1.0493
backward_s: 0.9180
rmse: 0.00273946
max_abs: 0.336487
relative_rmse: 0.000222724
finite: True
```

Full 137-level command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\roundtrip_error.py --chunk-size 64
```

Full 137-level result:

```text
input_shape: [137, 421120]
target_shape: [137, 938880]
roundtrip_shape: [137, 421120]
open_s: 6.46239
forward_s: 13.7963
backward_s: 12.9420
total_s: 33.2007
finite: True
count: 57693440
bias: -4.22452e-07
mean_abs: 0.0441958
max_abs: 8.52429
rmse: 0.14457
source_std: 31.3297
relative_rmse: 0.00461448
```

Interpretation:

- The full-column roundtrip is numerically stable and finite.
- The full-column error is much larger than the first four levels, so the next validation task should compute per-level error statistics before treating this stage as final.

Per-level report command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\roundtrip_error.py --chunk-size 64 --by-level --json > validation_roundtrip_o320_o480_o320_t_137.json
```

Largest per-level relative RMSE values:

```text
level=137 rmse=0.357148 max_abs=8.12909 relative_rmse=0.0236437
level=136 rmse=0.345497 max_abs=8.32877 relative_rmse=0.0233122
level=135 rmse=0.331991 max_abs=8.38217 relative_rmse=0.0226848
level=134 rmse=0.322113 max_abs=8.43134 relative_rmse=0.0222438
level=133 rmse=0.314871 max_abs=8.48563 relative_rmse=0.0219596
level=132 rmse=0.308727 max_abs=8.52429 relative_rmse=0.0217374
level=131 rmse=0.303587 max_abs=8.48856 relative_rmse=0.0215810
level=130 rmse=0.298995 max_abs=8.45291 relative_rmse=0.0214632
level=129 rmse=0.294415 max_abs=8.41910 relative_rmse=0.0213450
level=128 rmse=0.289807 max_abs=8.38626 relative_rmse=0.0212242
```

Interpretation:

- The largest roundtrip errors are concentrated in the lowest hybrid levels.
- Early upper levels have much smaller error; for example levels 1-4 have relative RMSE around `2e-4` to `3e-4`.
- Future regression thresholds should be per-variable and per-level-aware. A single global threshold hides where the error is coming from.

## Gate Before Next Feature

The spectral interpolation stage is ready to support the next feature only when:

- Unit tests pass in `skyborn_dev`.
- Real O320 to O480 interpolation works for one level and all 137 levels.
- Output xarray attrs describe the target grid, not the source grid.
- `chunk_size` works for small chunks, `64`, and `None`.
- Roundtrip error is recorded for the full 137-level field.
- Per-level roundtrip error is available, so large errors can be attributed to specific levels instead of hidden in one global metric.
- A regression threshold is chosen for future tests. Suggested initial thresholds for temperature roundtrip are global `relative_rmse < 0.01`, per-level `relative_rmse < 0.03`, finite output, and no shape/metadata mismatch.

The per-level report now exists, so the next feature can start after the regression thresholds above are accepted or adjusted.

Regression check command:

```powershell
F:\Anaconda3\envs\skyborn_dev\python.exe tools\check_regression.py
```

Latest result:

```text
Regression check passed:
  global relative_rmse=0.00461448
  worst level level=137 relative_rmse=0.0236437 rmse=0.357148 max_abs=8.12909
```
