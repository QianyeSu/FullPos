Spectral Filtering
==================

Backend
-------

Spectral filtering is implemented as:

.. code-block:: text

   grid values -> native ECTRANS spectral coefficients -> coefficient weights -> grid values

The transform steps use the native ECTRANS/FIAT backend. The filter profiles
are currently built in Python and applied to native ECTRANS coefficient arrays.

Triangular Truncation Filter
----------------------------

``spectral_filter`` applies a lower triangular truncation and returns values on
the same grid:

.. code-block:: python

   from fullpos import spectral_filter

   smooth = spectral_filter(
       ds["t"].isel(time=0),
       grid="O320",
       ntrunc=159,
       chunk_size=64,
   )

Generic Filter Profiles
-----------------------

``generic_spectral_filter`` supports FULLPOS-style diagonal coefficient
profiles:

.. code-block:: python

   from fullpos import generic_spectral_filter

   filtered = generic_spectral_filter(
       values,
       grid="O320",
       filter_kind="gaussian",
       ntrunc=319,
       cutoff=160,
       axis=-1,
   )

Supported profile kinds:

* ``"gaussian"``
* ``"low_pass"``

Direct Coefficient Filtering
----------------------------

When coefficients are already available, apply the profile directly:

.. code-block:: python

   from fullpos import filter_spectral_coefficients, spectral_fit

   coeffs = spectral_fit(values, grid="O320", ntrunc=319)
   filtered_coeffs = filter_spectral_coefficients(
       coeffs,
       filter_kind="low_pass",
       ntrunc=319,
       cutoff=160,
   )

Filter Matrix I/O
-----------------

Filter objects and matrix-like profiles can be saved and loaded:

.. code-block:: python

   from fullpos import load_spectral_filter, save_spectral_filter

   save_spectral_filter("filter.npz", spectral_filter_object)
   restored = load_spectral_filter("filter.npz")

Current Scope
-------------

The implemented filtering path covers generic spectral coefficient weighting,
Gaussian profiles, low-pass profiles, and read/write helpers. FULLPOS
post-processing level-specific filtering and stretched-geometry filtering are
not complete user-level APIs yet.
