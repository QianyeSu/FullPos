Vertical Interpolation
======================

Status
------

The public vertical API stays at ``fullpos.vertical``. The implementation has
been moved behind an internal ``fullpos._vertical`` package so new target
families can be added without growing one large module.

Current structure:

* ``fullpos.vertical``: stable public shim exported by the top-level package.
* ``fullpos._vertical.api``: public dispatch and capability reporting logic.
* ``fullpos._vertical.pos``: internal POS-style plan registry and dispatch
  layer for the implemented targets.
* ``fullpos._vertical.common``: shared target naming and validation.
* ``fullpos._vertical.pressure``: first implementation target for native
  hybrid-to-pressure interpolation.

Current implementation status:

* Public target names are defined and validated.
* ``target="pressure"`` uses the native OpenIFS/FULLPOS pressure-level
  backend. Dataset requests containing ``t``, ``u``, ``v``, and ``q`` use the
  native ``APACHE`` core for those variables; other pressure requests use the
  native PP-chain kernels directly.
* ``target="model_level"`` uses a native PP-chain column-pressure backend. It
  accepts target hybrid half-level coefficients and interpolates to the
  resulting per-column target full-level pressures. The target full-level
  pressures are computed inside the native Fortran wrapper.
* ``target="potential_temperature"`` uses native ``GPTET``/``PPLTETA`` to
  compute each target theta surface pressure, then interpolates fields through
  the same native ``PPQ``/``PPUV``/``PPT`` kernels.
* ``target="temperature"`` uses vendored native ``PPLTW`` and ``FPPS`` to
  compute each target temperature surface pressure, then interpolates fields
  through the same native ``PPQ``/``PPUV``/``PPT`` kernels.
* ``target="height_above_orography"`` uses native ``GPHPRE``/``GPGEO`` and
  ``FPPS`` to convert heights above the supplied surface geopotential into
  per-column target pressures, then interpolates fields through the same
  native ``PPQ``/``PPUV``/``PPT`` kernels.
* ``target="height_above_sea"`` uses the same native ``GPHPRE``/``GPGEO`` and
  ``FPPS`` path with target geopotential ``g * height_m``.
* ``target="flight_level"`` uses the same native absolute-height path as
  ``height_above_sea``. Public levels are standard flight-level numbers, so
  ``levels=[350]`` means ``FL350`` or 35,000 ft.
* ``target="potential_vorticity"`` uses native ``PPLTP`` on a provided PV
  field and Coriolis input, then interpolates fields through the same native
  ``PPQ``/``PPUV``/``PPT`` kernels.
* ``target="eta"`` uses native ``PPLETA``/``GPHPRE`` to compute target
  pressures for integer FULLPOS eta/model-level indexes, then interpolates
  fields through the same native ``PPQ``/``PPUV``/``PPT`` kernels.
* ``diagnose_potential_vorticity`` exposes the native FULLPOS ``GPPVO`` model
  level diagnostic. It can now auto-prepare relative vorticity, horizontal
  gradients, ``kappa`` from specific humidity, and Coriolis from Gaussian
  latitude metadata when those inputs are not supplied explicitly.
* The pressure-level request shape now validates ``levels``, dataset
  ``variables``, xarray ``chunks``, and hybrid/model-level dimensions.
* Pressure requests currently require surface pressure input
  (``surface_pressure=sp`` or ``surface_pressure=lnsp``).
* Hybrid A/B coefficients are taken from GRIB ``pv`` metadata when available,
  or from an external coefficients dataset passed through
  ``hybrid_coefficients=...``.
* For ERA5-style ``ap + b * ps`` coefficients, the native path uses the
  pressure-unit ``A`` coefficients directly.
* The compiled native pressure wrapper calls original OpenIFS/FULLPOS
  routines ``PPINIT``, ``PPFLEV``, ``PPQ``, ``PPUV``, ``PPT``, and ``PPSTA``.
  The potential-temperature target path also calls ``GPTET`` and ``PPLTETA``.
  The temperature target path calls ``PPLTW``, ``FPPS``, and ``PPPMER``.
  The height-above-orography target path calls ``GPHPRE``, ``GPGEO``,
  ``GPRCP`` when ``q`` is available, ``FPPS``, and ``PPPMER``.
  The height-above-sea target path uses the same native routines with absolute
  target geopotential.
  The flight-level target path converts ``FL`` numbers to metres above mean
  sea level and then uses the same native ``GPHPRE``/``GPGEO``/``FPPS`` path.
  The eta target path calls ``PPLETA`` and ``GPHPRE``.

Current boundary: the Python layer now owns a small POS-style plan registry and
dispatch layer for the implemented targets. Pressure Dataset interpolation for
``t``/``u``/``v``/``q`` is connected to the native ``APACHE`` core, and
pressure requests can opt into the APACHE ``LESCALE`` branch with
``lescale=True`` and ``target_surface_pressure=...``. This is still a Python
orchestration layer around native kernels, not a full ``pos.F90`` binding.

Native FULLPOS source boundary
------------------------------

The pressure-level path in OpenIFS/FULLPOS is not a single small routine.
The relevant chain is:

* ``VPOS`` / ``POS`` organize the grid-point vertical post-processing request.
* ``PPLETA`` computes pressure-level target pressure arrays for pressure-like
  vertical coordinates.
* ``APACHE`` performs the basic-field vertical interpolation.
* ``LESCALE`` adds the later post-processing and surface/derived-field
  handling that surrounds the core interpolation.
* ``GPHPRE``, ``GPGEO``, ``GPRCP``, and ``GPRH`` prepare hydrostatic pressure,
  geopotential, moist thermodynamic constants, and relative humidity context.

Completed
---------

* Pressure Dataset interpolation for ``t``/``u``/``v``/``q`` now reaches the
  native ``APACHE`` core.
* Pressure Dataset interpolation can enable the native APACHE ``LESCALE``
  branch with ``lescale=True`` and an optional ``target_surface_pressure``.
* Native pressure lookups for the implemented targets use the Fortran PP-chain
  path rather than a Python fallback.
* ``diagnose_potential_vorticity`` can auto-prepare the required diagnostic
  inputs when the Dataset provides the standard model-level fields.

Partial
-------

* The current backend is still a smaller native wrapper around APACHE and the
  PP-chain, not the full ``POS`` workflow.
* Several preparation routines are now wrapped directly for specific targets,
  including ``GPHPRE``, ``GPGEO``, ``GPRCP``, and ``FPPS`` for height targets.
* The full path still needs additional OpenIFS geometry, moist thermodynamic,
  relative-humidity, and ESCALE state such as ``FPVIEW``, ``GPRH``, and
  ``PPGEOP``.

Not complete
------------

* The complete ``POS`` model-level copy branch.
* Full ``POS``-driven surface and derived-field post-processing around the
  pressure requests.
* Surface-field special processing and the remaining geometry-dependent
  branches.
* Any claim that the current implementation is the exact full native vertical
  post-processing path.

Next milestones
---------------

The active vertical milestone is pressure-level interpolation. The immediate
goal is to keep ``target="pressure"`` stable for ERA5/OpenIFS model-level
Datasets, including the native ``APACHE`` path for ``t``/``u``/``v``/``q``.

The POS-style registry should remain a small internal extension point. Full
``POS`` orchestration, surface special cases, height/flight derived workflows,
and additional target families should be developed only when they become the
selected next feature.

Usage
-----

Use ``vertical_interpolate`` for native FULLPOS pressure-level interpolation of
model-level ERA5/OpenIFS fields:

.. code-block:: python

   from fullpos import vertical_interpolate

   out = vertical_interpolate(
       ds,
       target="pressure",
       levels=[100000, 85000, 50000],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
   )

Temperature uses ``PPT``, paired ``u``/``v`` variables use ``PPUV``, and other
scalar variables use ``PPQ``. The Python layer only validates metadata, chunks
leading dimensions, aligns surface pressure, and restores xarray metadata.

When the Dataset contains ``t``, ``u``, ``v``, and ``q``, those variables are
sent through native ``APACHE`` together. To enable APACHE ``LESCALE`` for a
pressure request, pass a target surface pressure field:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="pressure",
       levels=[30000, 50000, 85000],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
       target_surface_pressure=target_surface_ds["sp"],
       lescale=True,
   )

``target_surface_pressure`` must be an xarray ``DataArray`` aligned to the same
non-hybrid dimensions as ``surface_pressure``. If it is omitted, APACHE uses the
source surface pressure. ``lescale=True`` is only valid for Dataset inputs that
contain all four APACHE variables ``t``, ``u``, ``v``, and ``q``.

Model-level Target
------------------

Use ``target="model_level"`` when the target levels are represented by hybrid
half-level coefficients:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="model_level",
       variables=["t", "u", "v", "q"],
       target_hybrid_coefficients=target_coeffs,
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
   )

The current implementation passes target hybrid half-level coefficients to the
native wrapper, which computes per-column target full-level pressures in
Fortran and then calls the same native ``PPQ``/``PPUV``/``PPT`` kernels. It is
not the complete FULLPOS ``POS`` ``CDCONF='M'`` model-level copy branch. If the
target coefficients are the same as the source coefficients, ``PPQ`` scalar
fields can be nearly identity, but ``PPT`` and ``PPUV`` still apply their native
variable-specific logic.

Potential-temperature Target
----------------------------

Use ``target="potential_temperature"`` for theta surfaces. The target
``levels`` are potential temperature values in K:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="potential_temperature",
       levels=[300.0, 320.0, 340.0],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
   )

For ``Dataset`` input, variable ``t`` is used automatically to locate the
theta surfaces. For ``DataArray`` input, pass ``temperature=...`` unless the
array itself is the temperature field.

The target pressure calculation is native Fortran: ``GPTET`` computes the
full-level potential temperature profile and ``PPLTETA`` computes the pressure
of each requested theta surface. The current wrapper uses a dry-air kappa
constant, matching the minimal ``GPTET`` path. It does not yet wire the full
moist ``GPRCP`` thermodynamic context from ``POS/APACHE``.

Temperature Target
------------------

Use ``target="temperature"`` for isothermal surfaces. The target ``levels`` are
temperature values in K:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="temperature",
       levels=[250.0, 270.0],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
   )

For ``Dataset`` input, variable ``t`` is used automatically to locate the
temperature surfaces. For ``DataArray`` input, pass ``temperature=...`` unless
the array itself is the temperature field.

The current vendored source tree does not contain a usable ``PPLTEMP`` source;
both the local OpenIFS tree and the vendored tree only expose a dummy
``PPLTEMP`` that aborts. The implemented native path therefore uses the
available FULLPOS ``PPLTW`` routine to locate the target-temperature
geopotential and ``FPPS`` to convert that surface geopotential to pressure
before interpolation. This keeps the numerical path in Fortran/FULLPOS, but it
is not yet the exact missing ``PPLTEMP`` branch from full ``POS``.

Height Above Orography Target
-----------------------------

Use ``target="height_above_orography"`` for surfaces defined by geometric
height in metres above local orography:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="height_above_orography",
       levels=[0.0, 100.0, 500.0, 1000.0],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
       orography_geopotential=orography_ds["z"].isel(time=0, drop=True),
   )

The target ``levels`` are metres above orography. ``orography_geopotential``
must be the ECMWF surface geopotential in ``m2 s-2``. For ``Dataset`` input,
variable ``t`` is used automatically for the hydrostatic geopotential
calculation, and variable ``q`` is used automatically by native ``GPRCP`` to
compute moist ``R`` when it is present. If no ``q`` field is supplied, the
native wrapper uses dry-air ``R``.

The target pressure calculation is native Fortran: ``GPHPRE`` computes
half-level pressure and pressure metrics, ``GPGEO`` reconstructs half-level
geopotential from temperature and ``R``, and ``FPPS`` converts the requested
height surfaces to per-column pressures. The field interpolation then reuses
the same native ``PPQ``/``PPUV``/``PPT`` kernels as pressure-level output.

Height Above Sea Target
-----------------------

Use ``target="height_above_sea"`` for surfaces defined by absolute geometric
height in metres above mean sea level:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="height_above_sea",
       levels=[0.0, 100.0, 500.0, 1000.0],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
       orography_geopotential=orography_ds["z"].isel(time=0, drop=True),
   )

The input requirements match ``target="height_above_orography"``: temperature
is required for ``GPGEO`` and ECMWF surface geopotential ``z`` is required for
the source column geometry. The difference is the target surface definition:
``height_above_orography`` uses ``surface_geopotential + g * level`` while
``height_above_sea`` uses ``g * level``. Therefore the 0 m sea-level target is
not generally equal to the model surface pressure over land or below-sea
terrain.

Flight Level Target
-------------------

Use ``target="flight_level"`` for aircraft-style flight levels:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="flight_level",
       levels=[100.0, 200.0, 300.0, 350.0],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
       orography_geopotential=orography_ds["z"].isel(time=0, drop=True),
   )

The target ``levels`` are flight-level numbers in hundreds of feet. For
example, ``350`` means ``FL350`` or 35,000 ft, which is converted internally to
10,668 m above mean sea level. This conversion is only unit handling in the
Python layer; the target-pressure lookup and field interpolation still use the
native FULLPOS/OpenIFS Fortran path.

Internally, ``flight_level`` uses the same absolute-height pressure lookup as
``height_above_sea`` after converting ``FL`` values to metres. The wrapper then
interpolates fields with the native ``PPQ``/``PPUV``/``PPT`` kernels. The
output dimension is named ``flight_level`` and stores the original ``FL``
numbers.

Potential-vorticity Target
--------------------------

When opening ERA5/OpenIFS model-level GRIB files with ``cfgrib``, request the
hybrid A/B coefficient array explicitly:

.. code-block:: python

   ds = xr.open_dataset(
       "model-level.grib2",
       engine="cfgrib",
       backend_kwargs={"read_keys": ["pv"]},
   )

The resulting data variables should contain ``GRIB_pv`` in their attrs. FULLPOS
needs those coefficients for every vertical target.

Use ``target="potential_vorticity"`` for iso-PV surfaces. This currently
accepts either a supplied full-level PV field or a Dataset that contains the
native diagnostic inputs ``u``, ``v``, ``t``/``temperature``, and
``q``/``specific_humidity``. In the second case the wrapper first diagnoses
PV with ECTRANS + FULLPOS ``GPRCP``/``GPPVO`` and then locates the iso-PV
surface with ``PPLTP``:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="potential_vorticity",
       levels=[2.0e-6, 4.0e-6, 6.0e-6],
       variables=["t", "u", "v", "q"],
       surface_pressure=surface_ds["sp"],
       source_grid="O320",
       chunks={"time": 1},
   )

When ``potential_vorticity=...`` is supplied explicitly, the wrapper skips the
diagnostic step and uses that field directly. A Coriolis field may be supplied
with ``coriolis=...``; otherwise it is inferred from latitude coordinates or
from ``source_grid`` for supported Gaussian grids. After ``PPLTP`` locates
each iso-PV surface pressure, the wrapper reuses the same native
``PPQ``/``PPUV``/``PPT`` kernels for interpolation.

Potential-vorticity Diagnostic
------------------------------

Use ``diagnose_potential_vorticity`` when the goal is to compute model-level
PV and theta with native FULLPOS ``GPPVO`` before using
``target="potential_vorticity"``:

.. code-block:: python

   from fullpos import diagnose_potential_vorticity

   diag = diagnose_potential_vorticity(
       u=ds["u"],
       v=ds["v"],
       temperature=ds["t"],
       surface_pressure=surface_ds["sp"],
       specific_humidity=ds["q"],
       source_grid="N4",
   )

The numerical diagnostic path is native Fortran ``GPPVO``. When you only pass
``u``, ``v``, ``temperature``, ``surface_pressure``, and ``specific_humidity``,
the wrapper now prepares the missing native inputs through ECTRANS and
FULLPOS ``GPRCP``. Explicit ``relative_vorticity`` / gradient / ``kappa``
inputs still override the automatic preparation path, and ``coriolis`` is
derived automatically from Gaussian latitude metadata when it is not supplied.

Eta Target
----------

Use ``target="eta"`` for FULLPOS eta/model-level-index output. In the
OpenIFS/FULLPOS ``PPLETA`` interface, ``levels`` are 1-based integer indexes of
the output eta system, not arbitrary continuous sigma values:

.. code-block:: python

   out = vertical_interpolate(
       ds,
       target="eta",
       levels=[1, 10, 50, 100, 137],
       variables=["t", "u", "v", "q"],
       chunks={"time": 1, "values": 10000},
       surface_pressure=surface_ds["sp"],
   )

The wrapper computes target pressures with native ``PPLETA``/``GPHPRE`` and
then reuses native ``PPQ``/``PPUV``/``PPT`` for field interpolation. The output
dimension is named ``eta`` and stores the requested integer indexes.

Input validation helper
-----------------------

Use the local helper script to check that a model-level file, a surface
pressure file, and optional external hybrid coefficients form a valid pressure
request:

.. code-block:: powershell

   python tools/check_vertical_inputs.py model.grib2 surface.grib --levels 100000 85000 50000

Skyborn reference helper
------------------------

For development-only validation, use Skyborn's compiled
``interp_hybrid_to_pressure`` implementation as a reference on the same packed
ERA5/OpenIFS inputs. This does not change the FullPos runtime backend; it is
only a comparison tool.

For ERA5 coefficients the helper auto-detects ``p0=1`` because the hybrid
``A`` values are pressure-unit ``ap`` coefficients, not normalized ``a``
coefficients:

.. code-block:: powershell

   python tools/skyborn_pressure_reference.py model.grib2 surface.grib --variables t u v q --levels 100000 85000 50000

The helper can also write a NetCDF reference file:

.. code-block:: powershell

   python tools/skyborn_pressure_reference.py model.grib2 surface.grib --variables t --levels 30000 40000 50000 --method log --output skyborn_reference.nc

Candidate comparison helper
---------------------------

Compare native FULLPOS pressure output against an official pressure-level
product with:

.. code-block:: powershell

   python tools/compare_native_pressure_reference.py model.grib2 surface.grib truth.grib2 --variables t u v q --levels-hpa 200 300 400 500

Compare a native NetCDF candidate against the Skyborn development reference
with:

.. code-block:: powershell

   python tools/compare_pressure_to_skyborn.py model.grib2 surface.grib candidate.nc --variables t u v q --levels 30000 40000 50000 --method log

The comparison reports overall and per-level ``RMSE``, ``MAE``, ``max_abs``,
``bias``, and finite-point counts. Optional thresholds make it suitable for
local regression checks:

.. code-block:: powershell

   python tools/compare_pressure_to_skyborn.py model.grib2 surface.grib candidate.nc --variables t --levels 30000 40000 50000 --method log --max-rmse 0.5 --max-max-abs 10
