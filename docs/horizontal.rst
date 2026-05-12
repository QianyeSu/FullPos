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

For xarray inputs, the API detects ``latitude`` and ``longitude`` dimensions:

.. code-block:: python

   out = horizontal_interpolate(
       da,
       target_lats=target_lats,
       target_lons=target_lons,
       method="quadratic12",
       chunks={"time": 1, "hybrid": 10},
   )

``chunks`` uses named xarray dimensions. It controls how leading dimensions are
split before calling the native kernel. It is not used for NumPy inputs.

Packed Reduced Gaussian Usage
-----------------------------

Packed reduced input can be passed directly to the native horizontal path:

.. code-block:: python

   out = horizontal_interpolate(
       packed_values,
       source_grid="O96",
       target_lats=target_lats,
       target_lons=target_lons,
       method="bilinear",
   )

For xarray GRIB data, ``GRIB_pl`` is detected automatically:

.. code-block:: python

   out = horizontal_interpolate(
       da,
       target_lats=target_lats,
       target_lons=target_lons,
       method="quadratic12",
       chunks={"time": 1, "hybrid": 10},
   )

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
Mask-aware reduced-grid interpolation still needs the corresponding FULLPOS
mask-generation layer.

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
halo. Mask-aware interpolation with ``source_mask`` is still not implemented.
