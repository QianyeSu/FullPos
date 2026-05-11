Metadata and History
====================

Attribute Preservation
----------------------

Most xarray-facing APIs preserve variable attributes by default with
``keep_attrs=True``.

Dataset-level APIs copy safe non-horizontal coordinates and preserve dataset
attributes where possible. Variables without horizontal dimensions are skipped
by default unless the user explicitly requests them through ``variables=``.

History Entries
---------------

Regridding and filtering APIs append a CDO-style ``history`` entry to xarray
outputs. The newest entry is prepended above any existing history:

.. code-block:: text

   Sat May 09 23:30:00 2026: fullpos regrid source_grid=O320 target_grid=O480 method=spectral chunk_size=64 ntrunc=319 variables=t

The history line records:

* source grid
* target grid
* method
* chunk size
* truncation when supplied
* selected variables for Dataset operations

Current Limitations
-------------------

History support is implemented for the xarray regrid and spectral-filter APIs.
Lower-level NumPy functions return arrays and do not attach metadata.

The history entry is intentionally compact. It records what FullPos operation
was requested, not every native internal option used by ECTRANS or FULLPOS.
