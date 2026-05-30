API Reference
=============

Top-Level API
-------------

The top-level ``fullpos`` package re-exports the main user-facing functions and
classes from the modules below. Prefer the module sections for complete
docstrings and object indexes.

Common imports:

.. code-block:: python

   from fullpos import (
       Regridder,
       SpectralFilter,
       regrid,
       regrid_values,
       spectral_interpolate,
       spectral_interpolate_values,
       spectral_regrid,
       spectral_regrid_values,
       spectral_fit,
       spectral_synthesis,
       spectral_filter,
       generic_spectral_filter,
       horizontal_interpolate,
       land_sea_mask_to_grid,
       masked_surface_interpolate,
       capabilities,
       diagnose_potential_vorticity,
       to_grib,
       vertical_interpolate,
   )

Use ``capabilities()`` to inspect the current native feature boundary from
Python. It reports the native spectral, filtering, horizontal, and vertical
paths that are wired today, including the reduced-Gaussian note that user-level
horizontal interpolation should follow native spectral regridding.

For the explicit native ECTRANS/FIAT spectral path, prefer
``spectral_regrid(...)`` and ``spectral_regrid_values(...)``.
``spectral_interpolate(...)`` remains available as a readable alias.
``regrid(...)`` remains available as the broader compatibility entry point.

``regrid(...)`` also accepts regular latitude/longitude targets such as
``target_grid="LL1.0"`` and ``target_grid="LL0.25"`` for xarray inputs. Those
targets are dispatched through the native FULLPOS horizontal path and may
promote a Gaussian source through native spectral regridding first when the
interpolation stencil needs extra polar coverage. The ``LL`` outputs use
cell-center latitude coordinates and are not a Python fallback.
Current stable ``regrid`` support is Gaussian/O/F/N -> Gaussian/O/F/N and
Gaussian/O/F/N -> ``LL1.0`` / ``LL0.25``. ``LL`` is not supported as
``source_grid``.

For Gaussian-grid naming, ``F<N>`` is the canonical full regular Gaussian
name, while ``N<N>`` remains accepted as a compatibility alias for the same
regular Gaussian layout.

Spectral Module
---------------

.. automodule:: fullpos.spectral
   :members:

Filters
-------

.. automodule:: fullpos.filters
   :members:

Horizontal Interpolation
------------------------

.. automodule:: fullpos.interpolation.horizontal
   :members:

Surface Masks
-------------

.. automodule:: fullpos.masks
   :members:

Vertical Interpolation
----------------------

.. automodule:: fullpos.vertical
   :members:

GRIB I/O
--------

.. automodule:: fullpos.grib
   :members:

Low-Level Native Kernels
------------------------

.. automodule:: fullpos.interpolation.kernels
   :members:
