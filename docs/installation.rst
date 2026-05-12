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

Wheel builds copy the FIAT/ECTRANS shared libraries from the selected prefix
into ``fullpos/_native`` by default. This avoids depending on a source-tree
``extern/fullpos/local`` path after installation. Use the Meson option
``-Dfullpos_bundle_native_runtime=false`` for local builds that intentionally
keep all native runtime libraries external.

Native Libraries
----------------

On Windows, ``.dll`` files are dynamic libraries. The
``_ectrans.cp312-win_amd64.pyd`` file is also a dynamic library, but it is a
Python extension module imported directly by Python.

On Linux, native shared libraries usually use ``.so``. On macOS, they usually
use ``.dylib``.

FIAT and ECTRANS are bundled into wheels under ``fullpos/_native``. External
toolchain dependencies such as OpenBLAS and gfortran are handled by the platform
wheel repair tools: ``auditwheel`` on Linux, ``delocate`` on macOS, and
``delvewheel`` on Windows. Editable source installs can still use external
dependency mode through ``FULLPOS_NATIVE_PREFIX`` and the shell loader path.

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
