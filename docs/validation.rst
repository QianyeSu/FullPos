Validation
==========

Current validation scope
------------------------

The current stable path to validate is narrow:

* native spectral regridding,
* native spectral filtering, and
* native pressure-level and first-level vertical interpolation through the
  OpenIFS/FULLPOS-backed paths documented in ``vertical.rst``.

Unit Tests
----------

Run the test suite from the project root:

.. code-block:: powershell

   python -m pytest -q

The latest recorded full local run before this documentation page was written:

.. code-block:: text

   120 passed in 39.45s

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

The native regular-grid horizontal path was checked on real ERA5 F96 data.
Use a locally available decoded F96 ``regular_gg`` model-level temperature
field when reproducing the check.

.. code-block:: text

   input:  <f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: nearest
   RMSE:   0.0
   max:    0.0

.. code-block:: text

   input:  <f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: bilinear
   RMSE:   4.584093342540204e-15
   max:    1.1368683772161603e-13

.. code-block:: text

   input:  <f96-model-level-grib>
   field:  t(time=0, hybrid=0)
   method: quadratic12
   RMSE:   6.8512788032415326e-15
   max:    1.1368683772161603e-13

.. code-block:: text

   input:  <f96-model-level-grib>
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

The full JSON metric record from this run was written as a project-local
metrics artifact, and a Panoply-readable NetCDF sample was also produced for
the same pressure-grid comparison.

For reduced Gaussian inputs, validate the native spectral regrid first, then
apply user-level horizontal interpolation only after the field has been
converted to a regular-row grid. Packed reduced-Gaussian direct horizontal
cases remain the narrower native special cases documented in
``horizontal.rst``.

The real-data pytest smoke for this stage uses the same native pressure path
but does not require a fixed local file location. It reads the model-level and
surface inputs from ``FULLPOS_ERA5_MODEL_FILE`` and
``FULLPOS_ERA5_SURFACE_FILE`` when those environment variables are set; if
either file is missing the test is skipped. The check targets 200, 300, and
500 hPa for ``t/u/v/q``, verifies the native ``FULLPOS``/``APACHE`` attrs, and
uses an internal self-consistency comparison between the per-variable pressure
outputs.

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

The current native FULLPOS/OpenIFS-backed vertical validation path is the
pressure-level route: confirm that the Dataset input includes ``t``, ``u``,
``v``, and ``q``, provide surface pressure, and check the native ``APACHE``
attrs plus the pressure-level output against a known reference or the
built-in smoke tests. The same dataset path is used for the implemented
first-level targets, which remain wrappers around native pressure lookup and
interpolation rather than a full ``POS`` workflow.

Native temperature-surface calculation is covered by a low-level check that
calls ``temperature_pressures`` and verifies finite positive target pressures
from vendored ``PPLTW``/``FPPS``. Dataset-level tests verify
``target="temperature"`` uses dataset ``t`` to locate target surfaces and then
interpolates ``u``/``v``/``q`` through native PP-chain kernels.

Native height-above-orography target-pressure calculation is covered by a
low-level check that calls ``height_above_orography_pressures`` and verifies
finite pressures, exact agreement between the 0 m target and surface pressure,
and pressure decrease with increasing height. Dataset-level tests verify
``target="height_above_orography"`` uses dataset ``t`` and optional ``q``,
accepts ECMWF surface geopotential ``z``, broadcasts static orography over
time, and interpolates ``u``/``v``/``q`` through native PP-chain kernels.

The same path was checked on real ERA5 O96 reduced-grid data. A model-level
``(time=1, hybrid=137, values=40320)`` sample, matching surface pressure, and
O96 surface geopotential were processed at 0, 100, 500, 1000, and 3000 m above
orography:

.. code-block:: text

   backend: FULLPOS GPHPRE/GPGEO/FPPS + PP-chain
   output:  t/u/v/q(time=1, height_above_orography=5, values=40320)
   target pressure finite fraction: 1.0
   height 0 m surface-pressure max_abs_error: 0.0 Pa
   pressure strictly decreases with height fraction: 1.0
   target pressure range: 33238.64 .. 104378.45 Pa

The Panoply-readable sample output was recorded as a project-local NetCDF
artifact with a matching metrics JSON file.

Native height-above-sea target-pressure calculation is covered by the same
low-level and Dataset-level checks. The low-level check verifies finite
pressures, exact agreement between the 0 m sea-level target and surface
pressure for a zero-orography column, greater-than-surface pressure for a
positive-orography column at 0 m, and pressure decrease with increasing
absolute height.

The same ERA5 O96 sample was processed at 0, 100, 500, 1000, and 3000 m above
sea level:

.. code-block:: text

   backend: FULLPOS GPHPRE/GPGEO/FPPS + PP-chain
   output:  t/u/v/q(time=1, height_above_sea=5, values=40320)
   target pressure finite fraction: 1.0
   pressure strictly decreases with height fraction: 1.0
   target pressure range: 64941.76 .. 105808.47 Pa
   surface height range: -100.54 .. 5535.45 m

The Panoply-readable sample output was recorded as a project-local NetCDF
artifact with a matching metrics JSON file.

Native flight-level interpolation is covered by a Dataset-level check that
compares ``target="flight_level"`` against the equivalent
``target="height_above_sea"`` request after converting flight-level numbers to
metres. For example, ``FL10`` and ``FL20`` are checked against 304.8 m and
609.6 m above mean sea level. The numerical pressure lookup and field
interpolation remain on the native FULLPOS ``GPHPRE``/``GPGEO``/``FPPS`` plus
PP-chain path.

The same ERA5 O96 sample was processed at ``FL100``, ``FL200``, ``FL300``, and
``FL350`` for ``t/u/v/q``. The output was compared with
``target="height_above_sea"`` at 3048, 6096, 9144, and 10668 m:

.. code-block:: text

   backend: FULLPOS GPHPRE/GPGEO/FPPS + PP-chain
   output:  t/u/v/q(time=1, flight_level=4, values=40320)
   finite fraction: 1.0 for all variables
   equivalence RMSE: 0.0 for t/u/v/q
   equivalence max_abs: 0.0 for t/u/v/q

The Panoply-readable sample output was recorded as a project-local NetCDF
artifact with a matching metrics JSON file.

Skyborn's compiled ``interp_hybrid_to_pressure`` path was then used as a
development reference on the same O96 packed ERA5 input. For ERA5 the helper
uses ``p0=1`` because the hybrid ``A`` coefficients are pressure-unit
``ap/hyam`` values:

.. code-block:: text

   input:  <o96-model-level-grib>
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
