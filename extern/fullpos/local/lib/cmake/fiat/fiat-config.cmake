# Config file for the fiat package
# Defines the following variables:
#
#  fiat_FEATURES       - list of enabled features
#  fiat_VERSION        - version of the package
#  fiat_GIT_SHA1       - Git revision of the package
#  fiat_GIT_SHA1_SHORT - short Git revision of the package
#


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was project-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

### computed paths
set_and_check(fiat_CMAKE_DIR "${PACKAGE_PREFIX_DIR}/lib/cmake/fiat")
set_and_check(fiat_BASE_DIR "${PACKAGE_PREFIX_DIR}/.")
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(FIAT_CMAKE_DIR ${fiat_CMAKE_DIR})
  set(FIAT_BASE_DIR ${fiat_BASE_DIR})
endif()

### export version info
set(fiat_VERSION           "1.0.0")
set(fiat_GIT_SHA1          "")
set(fiat_GIT_SHA1_SHORT    "")

if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(FIAT_VERSION           "1.0.0" )
  set(FIAT_GIT_SHA1          "" )
  set(FIAT_GIT_SHA1_SHORT    "" )
endif()

### has this configuration been exported from a build tree?
set(fiat_IS_BUILD_DIR_EXPORT OFF)
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(FIAT_IS_BUILD_DIR_EXPORT ${fiat_IS_BUILD_DIR_EXPORT})
endif()

### include the <project>-import.cmake file if there is one
if(EXISTS ${fiat_CMAKE_DIR}/fiat-import.cmake)
  set(fiat_IMPORT_FILE "${fiat_CMAKE_DIR}/fiat-import.cmake")
  include(${fiat_IMPORT_FILE})
endif()

### insert definitions for IMPORTED targets
if(NOT fiat_BINARY_DIR)
  find_file(fiat_TARGETS_FILE
    NAMES fiat-targets.cmake
    HINTS ${fiat_CMAKE_DIR}
    NO_DEFAULT_PATH)
  if(fiat_TARGETS_FILE)
    include(${fiat_TARGETS_FILE})
  endif()
endif()

### include the <project>-post-import.cmake file if there is one
if(EXISTS ${fiat_CMAKE_DIR}/fiat-post-import.cmake)
  set(fiat_POST_IMPORT_FILE "${fiat_CMAKE_DIR}/fiat-post-import.cmake")
  include(${fiat_POST_IMPORT_FILE})
endif()

### handle third-party dependencies
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(FIAT_LIBRARIES         "")
  set(FIAT_TPLS              "" )

  include(${CMAKE_CURRENT_LIST_FILE}.tpls OPTIONAL)
endif()

### publish this file as imported
if( DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT )
  set(fiat_IMPORT_FILE ${CMAKE_CURRENT_LIST_FILE})
  mark_as_advanced(fiat_IMPORT_FILE)
  set(FIAT_IMPORT_FILE ${CMAKE_CURRENT_LIST_FILE})
  mark_as_advanced(FIAT_IMPORT_FILE)
endif()

### export features and check requirements
set(fiat_FEATURES "OMP;WARNINGS")
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(FIAT_FEATURES ${fiat_FEATURES})
endif()
foreach(_f ${fiat_FEATURES})
  set(fiat_${_f}_FOUND 1)
  set(fiat_HAVE_${_f} 1)
  if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
    set(FIAT_HAVE_${_f} 1)
  endif()
endforeach()
check_required_components(fiat)
