Horizontal Interpolation
========================

Backend Policy
--------------

The public horizontal interpolation API is reserved for native OpenIFS/FULLPOS
components. There is no public Python regular-lat/lon fallback.

Current native coverage:

* ``bilinear`` uses ``SUHOW1/SUHOW2`` address and weight generation plus
  ``FPINT4``.
* ``quadratic12`` uses ``SUHOW1/SUHOW2`` address and weight generation plus
  ``FPINT12``.
* ``nearest`` uses ``SUHOX1`` halo address generation plus ``FPNEAR``.
* ``average`` uses ``SUHOX1`` halo address generation plus ``FPAVG``.

For ``quadratic12``, the user-level API also exposes the native
``FPINT12`` monotonic clamp through ``shape_preserving=True``.
The explicit alias ``method="quadratic12_monotonic"`` is equivalent.
This is the preferred horizontal path for bounded positive scalar fields such
as specific humidity ``q`` and total column water vapour ``tcwv`` when spectral
ringing must be avoided.

Supported User-Level Input
--------------------------

The current user-level path supports:

* Global regular-row fields where every latitude row has the same longitude
  count. In practice this covers regular Gaussian ``N`` grids and regular
  Gaussian fields decoded as ``(latitude, longitude)``.
* Regular latitude/longitude ``LL`` targets for xarray inputs, including
  ``LL1.0`` and ``LL0.25``. The output uses cell-center latitude coordinates,
  and the wrapper stays on the native FULLPOS path.
* Packed reduced Gaussian fields when ``source_pl`` or ``source_grid="O<N>"``
  is supplied, or when xarray attributes contain ``GRIB_pl``.

Reduced Gaussian input must be brought onto a regular-row grid with native
spectral regridding before the user-level horizontal interpolation path is
used. The packed reduced-Gaussian cases documented below are the direct native
special cases.

``LL`` is output-only in the current stable surface. It is not accepted as a
``source_grid`` for ``regrid``.

Example:

.. code-block:: python

   import numpy as np
   import xarray as xr
   from fullpos import horizontal_interpolate

   lats = np.linspace(89.0, -89.0, 8)
   lons = np.linspace(0.0, 337.5, 16)
   da = xr.DataArray(
       np.ones((8, 16)),
       dims=("latitude", "longitude"),
       coords={"latitude": lats, "longitude": lons},
   )
   target_lats_1d = np.linspace(89.0, -89.0, 181)
   target_lons_1d = np.linspace(0.0, 359.0, 360)
   target_lons, target_lats = np.meshgrid(target_lons_1d, target_lats_1d)

   out = horizontal_interpolate(
       da,
       target_lats=target_lats,
       target_lons=target_lons,
       method="bilinear",
       chunks={"time": 1, "hybrid": 10},
   )

If the input is plain NumPy data, wrap it in an ``xarray.DataArray`` with
``latitude`` and ``longitude`` coordinates before calling the public API.

xarray DataArray Usage
----------------------

For xarray inputs, the API detects ``latitude`` and ``longitude`` dimensions.
You can either pass explicit target coordinates or a regular Gaussian target
grid such as ``F480``:

.. code-block:: python

   out = horizontal_interpolate(
       da,
       target_grid="F480",
       method="quadratic12_monotonic",
       chunks={"time": 1, "hybrid": 10},
   )

``shape_preserving=True`` or ``method="quadratic12_monotonic"`` is a direct wrapper around the original
OpenIFS/FULLPOS ``LDMONO`` branch in ``FPINT12``. It is only valid for
``method="quadratic12"``. It does not apply to spectral regridding.

``chunks`` uses named xarray dimensions. It controls how leading dimensions are
split before calling the native kernel. Each chunk is submitted as a native
FULLPOS ``KFIELDS`` batch, so ``chunks={"hybrid": 137}`` processes all model
levels in one geometry/weight setup when memory allows. Smaller chunks reduce
peak memory at the cost of more native calls. ``chunks`` is not used for NumPy
inputs.

Packed Reduced Gaussian Usage
-----------------------------

Packed reduced input can be passed directly to the native horizontal path:

.. code-block:: python

   out = horizontal_interpolate(
       packed_values,
       source_grid="O96",
       target_grid="F480",
       method="quadratic12_monotonic",
   )

For xarray GRIB data, ``GRIB_pl`` is detected automatically:

.. code-block:: python

   out = horizontal_interpolate(
       da,
       target_grid="F480",
       method="quadratic12_monotonic",
       chunks={"time": 1, "hybrid": 10},
   )

For ERA5 water-vapour workflows this keeps the calculation on the native
FULLPOS horizontal path:

.. code-block:: python

   q_f480 = horizontal_interpolate(
       model_ds,
       source_grid="N320",
       target_grid="F480",
       method="quadratic12_monotonic",
       variables=["q"],
       chunks={"time": 1, "hybrid": 137},
   )

   tcwv_f480 = horizontal_interpolate(
       surface_ds,
       source_grid="O320",
       target_grid="F480",
       method="quadratic12_monotonic",
       variables=["tcwv"],
       chunks={"time": 1},
   )

The output has regular Gaussian ``latitude`` and ``longitude`` dimensions and
updated ``GRIB_N``/``GRIB_gridType`` metadata. This path does not use spectral
interpolation.

Dataset Usage
-------------

Datasets can be processed variable by variable:

.. code-block:: python

   out = horizontal_interpolate(
       ds,
       target_lats=target_lats,
       target_lons=target_lons,
       method="bilinear",
       variables=["t", "u", "v"],
       chunks={"time": 1, "hybrid": 10},
   )

Variables without ``latitude`` and ``longitude`` dimensions are skipped by
default when ``variables`` is not explicit.

Reduced Gaussian Limits
-----------------------

Packed reduced input is supported for unmasked ``bilinear`` and
``quadratic12`` interpolation through ``SUHOW1/SUHOW2`` plus
``FPINT4``/``FPINT12``. It is also supported for unmasked ``nearest`` and
``average`` interpolation through ``SUHOX1`` plus ``FPNEAR``/``FPAVG``.
Mask-aware ``nearest`` and ``average`` are available through the native
``FPNEAR``/``FPAVG`` halo path by passing ``source_mask``. ``bilinear`` and
``quadratic12`` still use the unmasked ``FPINT4``/``FPINT12`` path.

Nearest, Average, and Masks
---------------------------

The native ``FPNEAR`` and ``FPAVG`` kernels are available through
``method="nearest"`` and ``method="average"``:

.. code-block:: python

   nearest = horizontal_interpolate(
       da,
       target_lats=target_lats,
       target_lons=target_lons,
       method="nearest",
   )

   averaged = horizontal_interpolate(
       da,
       target_lats=target_lats,
       target_lons=target_lons,
       method="average",
       average_radius=1,
   )

``average_radius`` maps to FULLPOS ``KSLWIDE``. The default ``1`` uses a 2x2
halo.

For bitmap surface fields, pass a boolean ``source_mask``. ``True`` marks
source points that can participate in interpolation. Points outside the mask
and non-finite source values are converted to the native FULLPOS undefined
value before ``FPAVG``/``FPNEAR`` is called; undefined native output is returned
to Python as ``NaN``. A boolean ``target_mask`` can also be passed to force
invalid target points, for example SST over land, to ``NaN``:

.. code-block:: python

   from fullpos import horizontal_interpolate, land_sea_mask_to_grid

   source_sea = land_sea_mask_to_grid(lsm, target_grid="O96", kind="sea")
   target_sea = land_sea_mask_to_grid(lsm, target_grid="F160", kind="sea")

   sst_f160 = horizontal_interpolate(
       sst,
       source_grid="O96",
       target_grid="F160",
       method="average",
       source_mask=source_sea,
       target_mask=target_sea,
       chunks={"time": 1},
   )

``land_sea_mask_to_grid`` samples a regular latitude/longitude land-sea mask
to an ``O``/``F``/``N`` Gaussian grid for mask preparation. ERA5 land-sea mask
convention is ``0`` over sea and ``1`` over land; ``kind="sea"`` uses
``lsm < 0.5``.

The convenience wrapper prepares both source and target masks, then dispatches
the field itself through native FULLPOS:

.. code-block:: python

   from fullpos import masked_surface_interpolate

   sst_f480 = masked_surface_interpolate(
       sst,
       land_sea_mask=lsm,
       source_grid="O96",
       target_grid="F480",
       kind="sea",
       method="average",
       chunks={"time": 1},
   )

Land points in the returned SST field are ``NaN``. This wrapper currently uses
the native ``FPAVG``/``FPNEAR`` mask-aware halo path. The separate native
FULLPOS physical 4/12-point land/sea weighting path remains future work.
