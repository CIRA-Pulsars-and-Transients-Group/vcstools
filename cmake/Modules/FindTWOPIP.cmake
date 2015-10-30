# - Try to find TWOPIP code.
# Variables used by this module:
#  TWOPIP_ROOT_DIR     - TWOPIP root directory
# Variables defined by this module:
#  TWOPIP_FOUND        - system has TWOPIP
#  TWOPIP_INCLUDE_DIR  - the TWOPIP include directory (cached)
#  TWOPIP_INCLUDE_DIRS - the TWOPIP include directories
#                         (identical to TWOPIP_INCLUDE_DIR)
#  TWOPIP_LIBRARY      - the TWOPIP library (cached)
#  TWOPIP_LIBRARIES    - the TWOPIP libraries
#                         (identical to TWOPIP_LIBRARY)

set(TWOPIP_ROOT_DIR $ENV{TWOPIP_DIR})

if(NOT DEFINED TWOPIP_ROOT_DIR)
	message("-- Warning TWOPIP_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED TWOPIP_ROOT_DIR)
	message("-- TWOPIP_ROOT_DIR = ${TWOPIP_ROOT_DIR}")
endif(NOT DEFINED TWOPIP_ROOT_DIR)

if(NOT TWOPIP_FOUND)

  find_path(TWOPIP_INCLUDE_DIR twopip_header.h
    HINTS ${TWOPIP_ROOT_DIR}/trunk/include PATH_SUFFIXES include include/twopip)
  find_library(TWOPIP_LIBRARY twopip
    HINTS ${TWOPIP_ROOT_DIR}/trunk/lib PATH_SUFFIXES lib)
  find_library(M_LIBRARY m)
  mark_as_advanced(TWOPIP_INCLUDE_DIR TWOPIP_LIBRARY M_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(TWOPIP DEFAULT_MSG
    TWOPIP_LIBRARY M_LIBRARY TWOPIP_INCLUDE_DIR)

  set(TWOPIP_INCLUDE_DIRS ${TWOPIP_INCLUDE_DIR})
  set(TWOPIP_LIBRARIES ${TWOPIP_LIBRARY} ${M_LIBRARY})

endif(NOT TWOPIP_FOUND)

if (TWOPIP_FOUND)
	message("-- Found TWOPIP --")
endif (TWOPIP_FOUND)


