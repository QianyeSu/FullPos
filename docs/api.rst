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
       spectral_fit,
       spectral_synthesis,
       spectral_filter,
       generic_spectral_filter,
       horizontal_interpolate,
       diagnose_potential_vorticity,
       vertical_interpolate,
   )

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

Vertical Interpolation
----------------------

.. automodule:: fullpos.vertical
   :members:

Low-Level Native Kernels
------------------------

.. automodule:: fullpos.interpolation.kernels
   :members:
