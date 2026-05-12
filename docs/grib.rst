GRIB Output
===========

``fullpos.to_grib`` writes xarray data to GRIB by cloning messages from an
existing GRIB template and replacing only the data values.

This is intentionally different from ``xarray.to_netcdf``. NetCDF stores a
self-describing dataset. GRIB is a stream of messages, and each message carries
grid, level, time, product definition, packing, bitmap, and discipline keys.
FullPos therefore requires a template file instead of trying to infer a safe
GRIB message from xarray metadata alone.

Basic Usage
-----------

.. code-block:: python

   import xarray as xr
   from fullpos import regrid, to_grib

   ds = xr.open_dataset(
       "input_o320.grib2",
       engine="cfgrib",
       backend_kwargs={
           "indexpath": "",
           "filter_by_keys": {"shortName": "t", "typeOfLevel": "hybrid"},
           "read_keys": ["gridType", "N", "pl", "numberOfPoints", "shortName"],
       },
   )

   out = regrid(ds["t"].isel(time=0), target_grid="O480")
   to_grib(out, "t_o480.grib2", template="template_o480_t.grib2")

The template must already describe the output grid and metadata that you want
in the file. For example, an ``O480`` output should use an ``O480`` template;
the writer will not convert an ``O320`` template into ``O480`` metadata.

Dataset Output
--------------

Dataset variables are matched to template messages by GRIB ``shortName``. A
variable's ``GRIB_shortName`` or ``shortName`` attribute takes precedence; the
xarray variable name is used otherwise.

.. code-block:: python

   to_grib(
       ds_out,
       "tuv_o480.grib2",
       template="template_o480_tuv.grib2",
       variables=["t", "u", "v"],
   )

By default, ``strict=True`` requires every written variable to have exactly the
same number of fields as the matching template messages. This prevents silently
writing a 10-level subset into a 137-level template with misleading level
metadata. Use ``strict=False`` only when you intentionally want to write the
leading subset of matching template messages.

Supported Shapes
----------------

``to_grib`` accepts:

* packed reduced-grid fields with a ``values`` dimension
* regular latitude/longitude or regular Gaussian fields with
  ``latitude``/``longitude`` or ``lat``/``lon`` dimensions
* leading dimensions such as ``time`` or ``hybrid``; these are flattened in
  xarray order and written as one GRIB message per 2D or packed field

Current Limitations
-------------------

``to_grib`` is an ecCodes-based I/O helper, not a FULLPOS interpolation kernel.
It does not change grids, levels, or packing metadata. Run ``regrid`` or
``vertical_interpolate`` first, then write using a template that already matches
the output geometry and message layout.

Install the optional GRIB dependencies before using this writer:

.. code-block:: powershell

   python -m pip install -e .[grib] --no-build-isolation
