#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "mpi_serial" for configuration "Release"
set_property(TARGET mpi_serial APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(mpi_serial PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libmpi_serial.a"
  )

list(APPEND _cmake_import_check_targets mpi_serial )
list(APPEND _cmake_import_check_files_for_mpi_serial "${_IMPORT_PREFIX}/lib/libmpi_serial.a" )

# Import target "fiat" for configuration "Release"
set_property(TARGET fiat APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fiat PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libfiat.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libfiat.dll"
  )

list(APPEND _cmake_import_check_targets fiat )
list(APPEND _cmake_import_check_files_for_fiat "${_IMPORT_PREFIX}/lib/libfiat.dll.a" "${_IMPORT_PREFIX}/bin/libfiat.dll" )

# Import target "parkind_sp" for configuration "Release"
set_property(TARGET parkind_sp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(parkind_sp PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libparkind_sp.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libparkind_sp.dll"
  )

list(APPEND _cmake_import_check_targets parkind_sp )
list(APPEND _cmake_import_check_files_for_parkind_sp "${_IMPORT_PREFIX}/lib/libparkind_sp.dll.a" "${_IMPORT_PREFIX}/bin/libparkind_sp.dll" )

# Import target "parkind_dp" for configuration "Release"
set_property(TARGET parkind_dp APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(parkind_dp PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libparkind_dp.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libparkind_dp.dll"
  )

list(APPEND _cmake_import_check_targets parkind_dp )
list(APPEND _cmake_import_check_files_for_parkind_dp "${_IMPORT_PREFIX}/lib/libparkind_dp.dll.a" "${_IMPORT_PREFIX}/bin/libparkind_dp.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
