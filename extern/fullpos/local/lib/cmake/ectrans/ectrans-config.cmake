# Config file for the ectrans package
# Defines the following variables:
#
#  ectrans_FEATURES       - list of enabled features
#  ectrans_VERSION        - version of the package
#  ectrans_GIT_SHA1       - Git revision of the package
#  ectrans_GIT_SHA1_SHORT - short Git revision of the package
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
set_and_check(ectrans_CMAKE_DIR "${PACKAGE_PREFIX_DIR}/lib/cmake/ectrans")
set_and_check(ectrans_BASE_DIR "${PACKAGE_PREFIX_DIR}/.")
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(ECTRANS_CMAKE_DIR ${ectrans_CMAKE_DIR})
  set(ECTRANS_BASE_DIR ${ectrans_BASE_DIR})
endif()

### export version info
set(ectrans_VERSION           "1.0.0")
set(ectrans_GIT_SHA1          "")
set(ectrans_GIT_SHA1_SHORT    "")

if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(ECTRANS_VERSION           "1.0.0" )
  set(ECTRANS_GIT_SHA1          "" )
  set(ECTRANS_GIT_SHA1_SHORT    "" )
endif()

### has this configuration been exported from a build tree?
set(ectrans_IS_BUILD_DIR_EXPORT OFF)
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(ECTRANS_IS_BUILD_DIR_EXPORT ${ectrans_IS_BUILD_DIR_EXPORT})
endif()

### include the <project>-import.cmake file if there is one
if(EXISTS ${ectrans_CMAKE_DIR}/ectrans-import.cmake)
  set(ectrans_IMPORT_FILE "${ectrans_CMAKE_DIR}/ectrans-import.cmake")
  include(${ectrans_IMPORT_FILE})
endif()

### insert definitions for IMPORTED targets
if(NOT ectrans_BINARY_DIR)
  find_file(ectrans_TARGETS_FILE
    NAMES ectrans-targets.cmake
    HINTS ${ectrans_CMAKE_DIR}
    NO_DEFAULT_PATH)
  if(ectrans_TARGETS_FILE)
    include(${ectrans_TARGETS_FILE})
  endif()
endif()

### include the <project>-post-import.cmake file if there is one
if(EXISTS ${ectrans_CMAKE_DIR}/ectrans-post-import.cmake)
  set(ectrans_POST_IMPORT_FILE "${ectrans_CMAKE_DIR}/ectrans-post-import.cmake")
  include(${ectrans_POST_IMPORT_FILE})
endif()

### handle third-party dependencies
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(ECTRANS_LIBRARIES         "")
  set(ECTRANS_TPLS              "" )

  include(${CMAKE_CURRENT_LIST_FILE}.tpls OPTIONAL)
endif()

### publish this file as imported
if( DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT )
  set(ectrans_IMPORT_FILE ${CMAKE_CURRENT_LIST_FILE})
  mark_as_advanced(ectrans_IMPORT_FILE)
  set(ECTRANS_IMPORT_FILE ${CMAKE_CURRENT_LIST_FILE})
  mark_as_advanced(ECTRANS_IMPORT_FILE)
endif()

### export features and check requirements
set(ectrans_FEATURES "OMP;DOUBLE_PRECISION;TRANSI")
if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
  set(ECTRANS_FEATURES ${ectrans_FEATURES})
endif()
foreach(_f ${ectrans_FEATURES})
  set(ectrans_${_f}_FOUND 1)
  set(ectrans_HAVE_${_f} 1)
  if(DEFINED ECBUILD_2_COMPAT AND ECBUILD_2_COMPAT)
    set(ECTRANS_HAVE_${_f} 1)
  endif()
endforeach()
check_required_components(ectrans)
