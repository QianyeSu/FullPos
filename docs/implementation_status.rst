Implementation Status
=====================

This page records the current development boundary. It should be updated when
new native FULLPOS functionality is wired into the Python API.

Completed
---------

* Native ECTRANS/FIAT spectral regridding through ``regrid`` and
  ``regrid_values``.
* xarray ``DataArray`` and ``Dataset`` spectral regridding with metadata and
  history preservation.
* Native spectral fitting through ``spectral_fit`` and ``spectral_fit_values``.
* Native spectral synthesis through ``spectral_synthesis`` and
  ``spectral_synthesis_values``.
* Spectral truncation filtering through ``spectral_filter``.
* Generic diagonal spectral coefficient filtering with Gaussian and low-pass
  profiles.
* Filter read/write helpers.
* Native low-level FULLPOS horizontal kernels: ``FPINT4``, ``FPINT12``,
  ``FPAVG``, and ``FPNEAR``.
* User-level regular-row native horizontal interpolation for ``bilinear`` and
  ``quadratic12`` via ``SUHOW1/SUHOW2`` plus ``FPINT4``/``FPINT12``.
* User-level packed reduced Gaussian horizontal interpolation for unmasked
  ``bilinear`` and ``quadratic12`` via ``SUHOW1/SUHOW2`` plus
  ``FPINT4``/``FPINT12``.
* User-level regular-row and packed reduced Gaussian horizontal interpolation
  for unmasked ``nearest`` and ``average`` via ``SUHOX1`` plus
  ``FPNEAR``/``FPAVG``.
* Named xarray-style ``chunks={...}`` handling for user-level horizontal
  interpolation.
* Native pressure-level vertical interpolation for xarray ERA5/OpenIFS-style
  model-level inputs through the OpenIFS/FULLPOS PP-chain:
  ``PPINIT``/``PPFLEV`` plus ``PPQ`` for scalar fields, paired ``PPUV`` for
  winds, and ``PPT``/``PPSTA`` for temperature.
* Dataset-level vertical pressure interpolation with ``variables=[...]``,
  surface pressure alignment, GRIB ``pv`` hybrid coefficient extraction,
  ``lnsp`` normalization, metadata/history preservation, and xarray-style
  ``chunks={...}``.
* Native PP-chain column-pressure backend for per-column target pressures.
* Native hybrid-target wrappers for ``target="model_level"``. Target hybrid
  half-level coefficients are converted to per-column full-level pressures in
  Fortran before calling the same ``PPQ``/``PPUV``/``PPT`` kernels.
* Native potential-temperature target pressures through OpenIFS/FULLPOS
  ``GPTET`` and ``PPLTETA``, followed by the same native
  ``PPQ``/``PPUV``/``PPT`` interpolation kernels.
* Native temperature target pressures through vendored FULLPOS ``PPLTW`` and
  ``FPPS``, followed by the same native ``PPQ``/``PPUV``/``PPT`` interpolation
  kernels.
* Native iso-PV target pressure lookup through FULLPOS ``PPLTP`` on a provided
  PV field and Coriolis input, followed by the same native
  ``PPQ``/``PPUV``/``PPT`` interpolation kernels.
* Native dependency diagnostics through ``backend_info`` and ``doctor``.

Partially Complete
------------------

* Mask-aware SST/surface-field handling is not a completed native FULLPOS
  user-level workflow.
* Pressure-level vertical interpolation is a native PP-chain wrapper, not yet
  the complete ``APACHE``/``LESCALE`` path. The complete path still needs
  OpenIFS geometry, moist thermodynamic, geopotential, and ESCALE dependencies
  such as ``FPVIEW``, ``GPHPRE``, ``GPRCP``, ``GPRH``, and ``PPGEOP``.
* ``target="model_level"`` currently means native PP-chain interpolation to
  target hybrid full-level pressures computed in the Fortran wrapper. It is not
  the complete FULLPOS ``POS`` ``CDCONF='M'`` model-level copy branch, so
  variable-specific kernels such as ``PPT`` and ``PPUV`` can change values even
  when source and target hybrid coefficients are the same.
* ``target="potential_temperature"`` currently uses the dry-air ``GPTET`` path
  to compute theta. It does not yet include the full moist ``GPRCP`` context
  used by complete ``POS/APACHE`` workflows.
* ``target="temperature"`` currently uses ``PPLTW``/``FPPS`` because this
  vendored FULLPOS/OpenIFS source set does not include a usable ``PPLTEMP``
  implementation. This is native Fortran, but not the exact missing
  ``PPLTEMP`` branch from full ``POS``.
* ``target="potential_vorticity"`` currently uses native ``PPLTP`` on a
  provided PV field and Coriolis input. The full automatic ``GPPVO`` diagnostic
  from ``u``/``v``/``t`` and geometry inputs is still pending.
* Development-only Skyborn compiled reference checks remain available for
  comparison, including correct ``p0=1`` handling for ERA5 ``ap``
  coefficients.
* Filter profiles are applied in Python to native ECTRANS coefficients. This is
  not yet a complete wrapper of every ``fpfilter.F90`` branch.

Not Complete
------------

* FULLPOS land/sea-mask weighted SST interpolation.
* Complete ``POS/APACHE``/``LESCALE`` vertical workflow, including the
  ``CDCONF='M'`` model-level branch.
* Full automatic ``GPPVO``-based potential-vorticity diagnostic from
  ``u``/``v``/``t`` and geometry inputs.
* Vertical interpolation to height, eta, and related levels.
* Surface-field special FULLPOS processing.
* Stretched-geometry filtering.
* Post-processing level-specific filtering.
* Wheel bundling of external native runtimes such as OpenBLAS and Fortran
  runtime libraries.

Current Grid Handling
---------------------

The spectral path supports ``O<N>`` octahedral reduced Gaussian grids and
``N<N>`` regular Gaussian grids.

The user-level horizontal path supports regular-row input and packed reduced
Gaussian input for unmasked ``bilinear``/``quadratic12``/``nearest``/``average``
interpolation. Real decoded examples:

* F96 ``regular_gg``: rectangular ``(192, 384)`` and usable by the current
  horizontal path.
* N96 ``reduced_gg``: packed classic reduced Gaussian field. Direct support
  requires providing compatible ``source_pl`` row lengths.
* O96 ``reduced_gg``: packed octahedral reduced Gaussian field, usable through
  ``source_grid="O96"`` or ``GRIB_pl``.

Next Practical Milestone
------------------------

The next major feature should either:

* add FULLPOS land/sea-mask weighted surface interpolation, or
* extend vertical interpolation from the current PP-chain wrapper toward the
  complete vendored OpenIFS/FULLPOS ``POS/APACHE`` chain.
