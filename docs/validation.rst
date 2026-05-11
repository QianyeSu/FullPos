Validation
==========

Unit Tests
----------

Run the test suite from the project root:

.. code-block:: powershell

   python -m pytest -q

The latest recorded full local run before this documentation page was written:

.. code-block:: text

   105 passed in 30.76s

Native Kernel Checks
--------------------

Low-level native horizontal kernels were checked against the OpenIFS/FULLPOS
formulas:

.. code-block:: text

   FPINT4  max_abs_error 0.0
   FPINT12 max_abs_error 0.0
   FPAVG   max_abs_error 0.0
   FPNEAR  max_abs_error 0.0

Real Data Checks
----------------

The native regular-grid horizontal path was checked on real ERA5 F96 data. The
exact local file path is machine-specific; the recorded sample used a decoded
F96 ``regular_gg`` model-level temperature field.

.. code-block:: text

   input:  <path-to-f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: nearest
   RMSE:   0.0
   max:    0.0

.. code-block:: text

   input:  <path-to-f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: bilinear
   RMSE:   4.584093342540204e-15
   max:    1.1368683772161603e-13

.. code-block:: text

   input:  <path-to-f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: quadratic12
   RMSE:   6.8512788032415326e-15
   max:    1.1368683772161603e-13

.. code-block:: text

   input:  <path-to-f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: average
   RMSE:   0.0
   max:    0.0

Native FULLPOS pressure-level interpolation was checked on real ERA5 O96
reduced-grid inputs. A model-level ``(time, hybrid, values)`` field plus
matching surface ``sp(time, values)`` were accepted, the GRIB ``pv`` vector was
split into 138 half-level ``A/B`` coefficients, and the hourly surface pressure
was aligned to the 6-hourly model times.

The native output was compared against the official ERA5 pressure-level product
on the same O96 grid for ``t/u/v/q`` at 200, 300, 400, and 500 hPa:

.. code-block:: text

   backend: FULLPOS PP-chain
   output:  t/u/v/q(time=4, pressure=4, values=40320)
   t overall RMSE: 0.188270 K,     max_abs: 7.17754 K
   u overall RMSE: 0.134141 m s-1, max_abs: 7.68178 m s-1
   v overall RMSE: 0.099900 m s-1, max_abs: 3.90626 m s-1
   q overall RMSE: 4.00771e-06,    max_abs: 5.84819e-04
   200 hPa t/u/v RMSE: 0.04056 / 0.03458 / 0.03500

The full JSON metric record from this run is
``fullpos_native_pressure_reference_metrics_19781201_o96_200_300_400_500.json``.
A Panoply-readable NetCDF sample was also written as
``fullpos_native_pressure_o96_tuvq_200_300_400_500.nc``.

The native hybrid-target wrapper was also smoke-tested on the same O96
model-level data for ``target="model_level"``. Target full-level pressures are
computed from target hybrid half-level coefficients in Fortran:

.. code-block:: text

   output: fullpos_native_model_level_o96_tuvq_time0_ml1-10.nc
   shape:  t/u/v/q(model_level=10, values=40320)
   max_abs(out - input): q=0, t=0.5664, u=4.01745, v=2.70177

The non-zero ``t/u/v`` differences are expected for this PP-chain path because
``PPT`` and ``PPUV`` apply FULLPOS variable-specific interpolation logic. This
check is not a validation of the complete ``POS`` ``CDCONF='M'`` copy branch.

Native potential-temperature target-pressure calculation is covered by a
low-level synthetic check. The test constructs a temperature profile from a
known theta profile, calls native ``GPTET``/``PPLTETA`` through
``theta_pressures``, and verifies the resulting pressure of internal theta
surfaces. Dataset-level tests also verify that ``target="potential_temperature"``
uses the dataset ``t`` variable to locate theta surfaces and then interpolates
``u``/``v`` through the native wind pair path.

Native temperature-surface calculation is covered by a low-level check that
calls ``temperature_pressures`` and verifies finite positive target pressures
from vendored ``PPLTW``/``FPPS``. Dataset-level tests verify
``target="temperature"`` uses dataset ``t`` to locate target surfaces and then
interpolates ``u``/``v``/``q`` through native PP-chain kernels.

Skyborn's compiled ``interp_hybrid_to_pressure`` path was then used as a
development reference on the same O96 packed ERA5 input. For ERA5 the helper
uses ``p0=1`` because the hybrid ``A`` coefficients are pressure-unit
``ap/hyam`` values:

.. code-block:: text

   input:  <path-to-o96-model-level-grib>
   field:  t(time, hybrid, values)
   levels: [100000, 85000, 50000] Pa
   p0:     1.0
   output: t(time=4, plev=3, values=40320)
   finite: 418732 / 483840
   100000 Pa finite: 103890 / 161280
   85000 Pa finite:  153562 / 161280
   50000 Pa finite:  161280 / 161280

The NaNs at 1000 hPa and part of 850 hPa are expected when
``extrapolate=False`` because some reduced-grid points sit above those target
pressure surfaces.

Validation Tools
----------------

Roundtrip validation:

.. code-block:: powershell

   python tools/roundtrip_error.py --max-levels 4 --chunk-size 4

Chunk-size benchmark:

.. code-block:: powershell

   python tools/benchmark_chunks.py --chunk-size 16 64 none

Regression threshold check:

.. code-block:: powershell

   python tools/check_regression.py

Skyborn vertical reference:

.. code-block:: powershell

   python tools/skyborn_pressure_reference.py model.grib2 surface.grib --variables t --levels 100000 85000 50000

Compare Skyborn reference against an official pressure-level product:

.. code-block:: powershell

   python tools/compare_pressure_reference.py model.grib2 surface.grib truth.grib2 --variables t u v q --levels-hpa 300 400 500 --method log --max-rmse 1

Compare native FULLPOS pressure output against an official pressure-level
product:

.. code-block:: powershell

   python tools/compare_native_pressure_reference.py model.grib2 surface.grib truth.grib2 --variables t u v q --levels-hpa 200 300 400 500 --output-json native_metrics.json

Compare a native FULLPOS vertical NetCDF output against the Skyborn development
reference:

.. code-block:: powershell

   python tools/compare_pressure_to_skyborn.py model.grib2 surface.grib candidate.nc --variables t u v q --levels 30000 40000 50000 --method log
