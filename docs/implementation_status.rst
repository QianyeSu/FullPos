Implementation Status
=====================

This page records the current development boundary. It should be updated when
new native FULLPOS functionality is wired into the Python API.

Current stable scope
--------------------

The documented stable path is intentionally narrow:

* native spectral regridding,
* native spectral filtering, and
* native pressure-level vertical interpolation through the current
  OpenIFS/FULLPOS-backed dataset paths.
* native Gaussian/O/F/N -> Gaussian/O/F/N and Gaussian/O/F/N -> ``LL1.0`` /
  ``LL0.25`` regridding. ``LL`` is output-only and is not supported as
  ``source_grid``.

Completed
---------

* Native ECTRANS/FIAT spectral regridding through ``regrid`` and
  ``regrid_values``.
* xarray ``DataArray`` and ``Dataset`` spectral regridding with metadata and
  history preservation.
* Native regular latitude/longitude ``LL1.0`` and ``LL0.25`` output through
  ``regrid(...)`` on the FULLPOS horizontal path, with cell-center latitude
  coordinates and no Python fallback.
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
* Python-side POS-style vertical orchestration through a small internal plan
  registry for the implemented targets.
* Dataset-level vertical pressure interpolation with ``variables=[...]``,
  surface pressure alignment, GRIB ``pv`` hybrid coefficient extraction,
  ``lnsp`` normalization, metadata/history preservation, and xarray-style
  ``chunks={...}``.
* Dataset-level ``t``/``u``/``v``/``q`` pressure interpolation through native
  ``APACHE``, with explicit APACHE ``LESCALE`` opt-in via ``lescale=True`` and
  ``target_surface_pressure=...``.
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
* Native height-above-orography target pressure lookup through FULLPOS
  ``GPHPRE``/``GPGEO``/``FPPS``. The wrapper accepts heights in metres above
  ECMWF surface geopotential, uses native ``GPRCP`` for moist ``R`` when
  ``q`` is available, and then reuses the same native
  ``PPQ``/``PPUV``/``PPT`` interpolation kernels.
* Native height-above-sea target pressure lookup through the same FULLPOS
  ``GPHPRE``/``GPGEO``/``FPPS`` path with target geopotential
  ``g * height_m``.
* Native flight-level target pressure lookup through the same absolute-height
  FULLPOS ``GPHPRE``/``GPGEO``/``FPPS`` path. Public levels are standard
  ``FL`` numbers in hundreds of feet, converted internally to metres above
  mean sea level before calling the native pressure lookup.
* Native iso-PV target pressure lookup through FULLPOS ``PPLTP`` on either a
  provided PV field or an automatically diagnosed native ``GPPVO`` PV field,
  followed by the same native ``PPQ``/``PPUV``/``PPT`` interpolation kernels.
* Native eta/model-level-index target pressure lookup through FULLPOS
  ``PPLETA``/``GPHPRE``, followed by the same native
  ``PPQ``/``PPUV``/``PPT`` interpolation kernels.
* Native model-level potential-vorticity and potential-temperature diagnostic
  through FULLPOS ``GPPVO``. Missing relative vorticity, temperature
  gradients, surface-pressure gradients, and ``kappa`` (``R/Cp``) are now
  prepared natively through ECTRANS and FULLPOS ``GPRCP`` when the wrapper is
  given ``u``/``v``/``t``/``q`` plus grid metadata.
* Native dependency diagnostics through ``backend_info`` and ``doctor``.

Partial
-------

* Mask-aware SST/surface-field handling is not a completed native FULLPOS
  user-level workflow.
* Pressure-level Dataset interpolation for ``t``/``u``/``v``/``q`` is wired
  into the native ``APACHE`` core, including explicit APACHE ``LESCALE`` opt-in
  for pressure requests. The Python side now routes targets through a small
  internal POS-style registry/plan layer, but it is still not the full
  ``pos.F90`` workflow. Some preparation routines are still wrapped directly
  for specific targets, including ``GPHPRE``/``GPGEO`` and ``GPRCP`` for
  ``height_above_orography``. The remaining path still needs additional
  OpenIFS geometry and surface-processing state to expose full ``POS``.
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
* ``target="potential_vorticity"`` uses native ``PPLTP`` and can now
  auto-diagnose the required PV field through native ECTRANS + FULLPOS
  ``GPRCP``/``GPPVO`` when the input Dataset contains ``u``/``v``/``t``/``q``
  and Gaussian grid metadata.
* ``target="eta"`` currently follows the native ``PPLETA`` convention where
  requested levels are integer eta/model-level indexes. It is not a continuous
  0..1 sigma-coordinate interpolation API.
* ``target="height_above_orography"``, ``target="height_above_sea"``, and
  ``target="flight_level"`` currently implement the native FULLPOS
  ``GPHPRE``/``GPGEO``/``FPPS`` pressure lookups plus the native PP-chain
  interpolation kernels. They are not yet the complete ``APACHE`` output-field
  workflow for all derived fields.
* Reduced Gaussian fields still need spectral regridding before the regular-row
  user-level horizontal interpolation path is used. The packed reduced-Gaussian
  native routes documented above remain the narrower special cases.
* Development-only Skyborn compiled reference checks remain available for
  comparison, including correct ``p0=1`` handling for ERA5 ``ap``
  coefficients.
* Filter profiles are applied in Python to native ECTRANS coefficients. This is
  not yet a complete wrapper of every ``fpfilter.F90`` branch.

Not complete
------------

* FULLPOS land/sea-mask weighted SST interpolation.
* Complete ``POS`` vertical workflow, including the ``CDCONF='M'`` model-level
  branch and the remaining geometry and surface processing hooks.
* Surface-field special FULLPOS processing.
* Stretched-geometry filtering.
* Post-processing level-specific filtering.
* Fully proven release publishing on all wheel targets. The wheel build
  workflow now prepares a native FIAT/ECTRANS prefix and bundles those shared
  libraries, but the repaired wheels still need to be validated on GitHub
  Actions before a public release is cut.

Current Grid Handling
---------------------

The spectral path supports ``O<N>`` octahedral reduced Gaussian grids and
``F<N>`` full regular Gaussian grids. ``N<N>`` remains accepted as a
compatibility alias for the regular Gaussian layout.

The user-level horizontal path supports regular-row input and packed reduced
Gaussian input for unmasked ``bilinear``/``quadratic12``/``nearest``/``average``
interpolation. Real decoded examples:

* F96 ``regular_gg``: rectangular ``(192, 384)`` and usable by the current
  horizontal path.
* N96 ``reduced_gg``: packed classic reduced Gaussian field. Direct support
  requires providing compatible ``source_pl`` row lengths.
* O96 ``reduced_gg``: packed octahedral reduced Gaussian field, usable through
  ``source_grid="O96"`` or ``GRIB_pl``.
Next milestones
---------------

The current development milestone is intentionally narrow:

* keep the native spectral regridding API stable,
* keep the native spectral filtering API stable, and
* stabilize the native pressure-level vertical path, especially the
  Dataset-level ``APACHE`` route for ``t``/``u``/``v``/``q``.

The internal POS-style registry is only an extension point for later work.
Surface special processing, full ``POS`` orchestration, height/flight derived
post-processing, and additional vertical targets should stay out of the active
milestone unless they are explicitly selected as the next feature.
