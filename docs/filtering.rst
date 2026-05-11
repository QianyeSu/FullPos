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

Completed
---------

* Spectral filtering through ``spectral_filter`` and
  ``generic_spectral_filter``.
* Direct coefficient filtering through ``filter_spectral_coefficients``.
* Filter read/write helpers for stored filter objects.

Partial
-------

* The current Python implementation covers generic spectral coefficient
  weighting and the common Gaussian and low-pass profiles.
* The code operates on native ECTRANS coefficient arrays, but the filter
  profiles are still assembled in Python rather than by a full wrapper around
  every native ``fpfilter.F90`` branch.

Not complete
------------

* ``stretched geometry filter`` means the FULLPOS filtering branch that is
  aware of stretched-grid geometry and the associated geometry-specific
  interpolation rules. That branch is not exposed as a user-level API yet.
* ``post-processing level-specific filtering`` means the FULLPOS branch where
  the applied filter depends on the post-processing level or output level
  family. That branch is also not complete in the current Python wrapper.
* FULLPOS parity for every ``fpfilter.F90`` branch.

Next milestones
---------------

* Add a native wrapper for the remaining ``fpfilter.F90`` branches.
* Define the stretched-geometry path explicitly at the user API boundary.
* Add level-specific filter selection once the corresponding FULLPOS semantics
  are exposed cleanly.
