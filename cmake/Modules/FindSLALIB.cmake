# - Try to find SLALIB code.
# Variables used by this module:
#  SLALIB_ROOT_DIR     - SLALIB root directory
# Variables defined by this module:
#  SLALIB_FOUND        - system has SLALIB
#  SLALIB_INCLUDE_DIR  - the SLALIB include directory (cached)
#  SLALIB_INCLUDE_DIRS - the SLALIB include directories
#                         (identical to SLALIB_INCLUDE_DIR)
#  SLALIB_LIBRARY      - the SLALIB library (cached)
#  SLALIB_LIBRARIES    - the SLALIB libraries
#                         (identical to SLALIB_LIBRARY)

set(SLALIB_ROOT_DIR $ENV{SLALIB_DIR})

if(NOT DEFINED SLALIB_ROOT_DIR)
	message("-- Warning SLALIB_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED SLALIB_ROOT_DIR)
	message("-- SLALIB_ROOT_DIR = ${SLALIB_ROOT_DIR}")
endif(NOT DEFINED SLALIB_ROOT_DIR)

if(NOT SLALIB_FOUND)

  find_path(SLALIB_INCLUDE_DIR slalib.h
    HINTS ${SLALIB_ROOT_DIR} PATH_SUFFIXES include )
  find_library(SLALIB_LIBRARY sla
    HINTS ${SLALIB_ROOT_DIR} PATH_SUFFIXES lib )
  mark_as_advanced(SLALIB_INCLUDE_DIR SLALIB_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SLALIB DEFAULT_MSG
    SLALIB_LIBRARY SLALIB_INCLUDE_DIR)

  set(SLALIB_INCLUDE_DIRS ${SLALIB_INCLUDE_DIR})
  set(SLALIB_LIBRARIES ${SLALIB_LIBRARY} ${M_LIBRARY})

endif(NOT SLALIB_FOUND)

if (SLALIB_FOUND)
	message("-- Found SLALIB --")
endif (SLALIB_FOUND)


