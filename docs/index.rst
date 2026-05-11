FullPos Documentation
=====================

FullPos is a Python package that wraps ECMWF/OpenIFS FULLPOS-style
interpolation components for local Python workflows. The package is currently
focused on native ECTRANS/FIAT spectral transforms and selected native FULLPOS
horizontal interpolation kernels.

The public API intentionally avoids silent Python fallbacks for core FULLPOS
features. When a feature is not wired to the native OpenIFS/FULLPOS path yet,
the documentation marks it as not implemented.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   spectral
   filtering
   horizontal
   vertical
   metadata_history
   validation
   implementation_status
   api

Build These Docs
----------------

Install the documentation dependencies and build HTML output from the project
root:

.. code-block:: powershell

   python -m pip install -e .[docs] --no-build-isolation
   python -m sphinx -b html docs docs/_build/html

The generated entry point is ``docs/_build/html/index.html``.
