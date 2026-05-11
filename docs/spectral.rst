Spectral Regridding
===================

Backend
-------

Spectral regridding uses the native ECTRANS/FIAT path exposed by the
``fullpos._ectrans`` extension. It does not fall back to Skyborn or another
Python spectral transform implementation.

Supported Gaussian grids:

* ``O<N>`` octahedral reduced Gaussian grids, stored as packed 1D fields.
* ``N<N>`` regular Gaussian grids, stored as ``(latitude, longitude)`` fields.

Basic xarray Usage
------------------

.. code-block:: python

   import xarray as xr
   from fullpos import regrid

   ds = xr.open_dataset(
       "<path-to-o320-model-level-grib>",
       engine="cfgrib",
       backend_kwargs={
           "indexpath": "",
           "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
           "read_keys": ["gridType", "N", "pl", "numberOfPoints", "packingType"],
       },
   )

   out = regrid(ds["t"].isel(time=0), target_grid="O480", chunk_size=64)

When GRIB metadata contains ``GRIB_gridType``, ``GRIB_N``, and ``GRIB_pl``,
``source_grid`` can be inferred. Otherwise pass it explicitly:

.. code-block:: python

   out = regrid(ds["t"], source_grid="O320", target_grid="O480", chunk_size=16)

NumPy Usage
-----------

Use ``regrid_values`` when the grid is already known:

.. code-block:: python

   from fullpos import regrid_values

   result = regrid_values(
       values,
       source_grid="O320",
       target_grid="O480",
       axis=-1,
       chunk_size=64,
   )

For reduced grids, ``axis`` identifies the packed horizontal dimension. For
regular Gaussian grids, pass ``axis=(lat_axis, lon_axis)`` when the horizontal
dimensions are not the last two axes.

Missing Values
--------------

The spectral path requires finite global input by default. Fields decoded from
GRIB bitmaps, for example SST over land, contain missing values and are rejected
with ``missing_policy="error"``.

Use ``missing_policy="ignore"`` only when debugging native NaN propagation. It
does not make spectral interpolation physically mask-aware.

Repeated Calls
--------------

For repeated operations on the same grid pair, reuse ``Regridder``:

.. code-block:: python

   from fullpos import Regridder

   regridder = Regridder("O320", "O480", chunk_size=64)
   out = regridder.regrid_data_array(ds["t"].isel(time=0))

Spectral Fit and Synthesis
--------------------------

``spectral_fit`` converts grid-point values to the native ECTRANS global
real/imaginary spectral coefficient layout. ``spectral_synthesis`` converts
those coefficients back to grid-point values.

.. code-block:: python

   from fullpos import spectral_fit, spectral_synthesis

   coeffs = spectral_fit(values, grid="O320", ntrunc=319, axis=-1)
   restored = spectral_synthesis(coeffs, grid="O320", ntrunc=319, axis=-1)

The last coefficient dimension has length ``(ntrunc + 1) * (ntrunc + 2)``.
