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
* ``fullpos._vertical.common``: shared target naming and validation.
* ``fullpos._vertical.pressure``: first implementation target for native
  hybrid-to-pressure interpolation.

Current implementation status:

* Public target names are defined and validated.
* ``target="pressure"`` uses the native OpenIFS/FULLPOS pressure-level
  PP-chain backend.
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
  The eta target path calls ``PPLETA`` and ``GPHPRE``.

Native FULLPOS source boundary
------------------------------

The pressure-level path in OpenIFS/FULLPOS is not a single small routine.
The relevant chain is:

* ``VPOS`` / ``POS`` organize the grid-point vertical post-processing request.
* ``PPLETA`` computes pressure-level target pressure arrays for pressure-like
  vertical coordinates.
* ``APACHE`` performs the basic-field vertical interpolation.
* ``GPHPRE``, ``GPGEO``, ``GPRCP``, and ``GPRH`` prepare hydrostatic pressure,
  geopotential, moist thermodynamic constants, and relative humidity context.

The current backend is intentionally a smaller native PP-chain wrapper, not the
complete ``APACHE``/``LESCALE`` workflow. ``APACHE`` also pulls in OpenIFS
geometry, moist thermodynamic, geopotential, and ESCALE state such as
``FPVIEW``, ``GPHPRE``, ``GPRCP``, ``GPRH``, and ``PPGEOP``. Those dependencies
are still outside the Python extension boundary.

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
