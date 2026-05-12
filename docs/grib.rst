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

Packing and Compression
-----------------------

By default, output messages keep the template ``packingType``. If the template
uses ``grid_ccsds``, the output stays CCSDS-packed. If the template uses
``grid_simple``, the output stays simple-packed.

You can request an explicit ecCodes packing override:

.. code-block:: python

   to_grib(
       out,
       "t_o480_ccsds.grib2",
       template="template_o480_t.grib2",
       packing_type="ccsds",
       bits_per_value=16,
   )

``packing_type="simple"`` maps to ``grid_simple``. ``packing_type="ccsds"``
or ``"aec"`` maps to ``grid_ccsds``. You may also pass an ecCodes ``grid_*``
packing name directly.

CCSDS/AEC packing is a GRIB edition 2 feature and requires an ecCodes build
with libaec support. If the active ecCodes library cannot encode the requested
packing, ``to_grib`` will raise the underlying ecCodes error.

Metadata Overrides
------------------

By default, ``to_grib`` preserves the template metadata. Optional overrides are
available for common GRIB keys:

.. code-block:: python

   to_grib(
       out,
       "t_o480_20250102.grib2",
       template="template_o480_t.grib2",
       centre="ecmf",
       sub_centre=0,
       generating_process_identifier=255,
       data_date=20250102,
       data_time=600,
       step_type="instant",
       forecast_time=0,
       type_of_level="isobaricInhPa",
       level=500,
   )

The keyword names use Python style and map to ecCodes keys:

* ``centre`` -> ``centre``
* ``sub_centre`` -> ``subCentre``
* ``generating_process_identifier`` -> ``generatingProcessIdentifier``
* ``data_date`` -> ``dataDate``
* ``data_time`` -> ``dataTime``
* ``step_type`` -> ``stepType``
* ``step_range`` -> ``stepRange``
* ``forecast_time`` -> ``forecastTime``
* ``param_id`` -> ``paramId``
* ``short_name`` -> ``shortName``
* ``type_of_level`` -> ``typeOfLevel``
* ``level`` -> ``level``
* ``edition`` -> ``edition``

``param_id`` and ``short_name`` are mutually exclusive. ``step_range`` and
``forecast_time`` are also mutually exclusive. ``param_id`` and ``short_name``
are also used to select matching template messages, so they can only be used
with a single DataArray or one selected Dataset variable.

For advanced ecCodes keys, use ``key_overrides``:

.. code-block:: python

   to_grib(
       out,
       "t_o480_custom.grib2",
       template="template_o480_t.grib2",
       key_overrides={"localDefinitionNumber": 1},
   )

``key_overrides`` is applied after the named arguments, so it can intentionally
override them. Avoid changing geometry keys such as ``numberOfPoints`` or
``pl``; use a matching template instead.

GRIB1 and GRIB2
---------------

``to_grib`` preserves the template edition by default. A GRIB2 template writes
GRIB2 output; a GRIB1 template writes GRIB1 output. You can pass
``edition=2`` or ``edition=1`` to ask ecCodes to convert the cloned message,
but a matching template in the desired edition is safer.

GRIB1 is the older format and has fewer product-definition options. GRIB2 has
more explicit templates for products, ensembles, probability fields, chemical
and wave products, and modern packing such as CCSDS/AEC. For new FullPos
outputs, use GRIB2 templates unless you need compatibility with an older GRIB1
workflow.

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
It does not change grids or levels. Run ``regrid`` or ``vertical_interpolate``
first, then write using a template that already matches the output geometry and
message layout. Packing can be preserved from the template or changed through
``packing_type``.

Install the optional GRIB dependencies before using this writer:

.. code-block:: powershell

   python -m pip install -e .[grib] --no-build-isolation
