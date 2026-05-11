Installation
============

Development Environment
-----------------------

Install the package in an environment with a working C/Fortran toolchain and
the required native FIAT/ECTRANS runtime libraries:

.. code-block:: powershell

   python -m pip install -e . --no-build-isolation

For documentation builds:

.. code-block:: powershell

   python -m pip install -e .[docs] --no-build-isolation

Native Dependency Prefix
------------------------

FullPos expects FIAT and ECTRANS runtime libraries to be available. The default
development prefix is:

.. code-block:: text

   extern/fullpos/local

At build time, a different native installation can be selected with the Meson
option ``fullpos_native_prefix``:

.. code-block:: powershell

   python -m pip install -e . --no-build-isolation --no-deps --config-settings=setup-args="-Dfullpos_native_prefix=<path-to-native-prefix>"

At runtime, ``FULLPOS_NATIVE_PREFIX`` can override the configured prefix for
diagnostics and local testing.

Native Libraries
----------------

On Windows, ``.dll`` files are dynamic libraries. The
``_ectrans.cp312-win_amd64.pyd`` file is also a dynamic library, but it is a
Python extension module imported directly by Python.

On Linux, native shared libraries usually use ``.so``. On macOS, they usually
use ``.dylib``.

The current development mode is external dependency mode: FullPos checks that
the runtime libraries are discoverable, but it does not bundle OpenBLAS,
gfortran, FIAT, or ECTRANS runtime libraries into wheels yet.

Runtime Diagnostics
-------------------

Use the doctor command to check whether the native backend can be imported:

.. code-block:: powershell

   python tools/doctor.py

The same information is available from Python:

.. code-block:: python

   from fullpos import backend_info, doctor

   print(backend_info())
   doctor()
