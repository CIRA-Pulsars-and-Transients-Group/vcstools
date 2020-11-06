# - Try to find HYPERBEAM code.
# Variables used by this module:
#  HYPERBEAM_ROOT     - HYPERBEAM root directory
# Variables defined by this module:
#  HYPERBEAM_FOUND        - system has HYPERBEAM
#  HYPERBEAM_INCLUDE_DIR  - the HYPERBEAM include directory (cached)
#  HYPERBEAM_INCLUDE_DIRS - the HYPERBEAM include directories
#                         (identical to HYPERBEAM_INCLUDE_DIR)
#  HYPERBEAM_LIB          - the HYPERBEAM library (cached)
#  HYPERBEAM_LIBRARIES    - the HYPERBEAM libraries
#  HYPERBEAM_HDF5         - the full embedded element pattern HDF5 file

message("Finding HYPERBEAM")

set(HYPERBEAM_ROOT $ENV{HYPERBEAM})

if(NOT DEFINED HYPERBEAM_ROOT)
    message(STATUS "Warning HYPERBEAM_ROOT not set: will try and find it ")
else(NOT DEFINED HYPERBEAM_ROOT)
    message(STATUS "HYPERBEAM_ROOT = ${HYPERBEAM_ROOT}")
endif(NOT DEFINED HYPERBEAM_ROOT)

if(NOT HYPERBEAM_FOUND)

  find_path(HYPERBEAM_INCLUDE_DIR mwa_hyperbeam.h
    HINTS ${HYPERBEAM_ROOT} PATH_SUFFIXES /include /include/mwa_hyperbeam)
  find_library(HYPERBEAM_LIB mwa_hyperbeam
    HINTS ${HYPERBEAM_ROOT} PATH_SUFFIXES lib )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(HYPERBEAM DEFAULT_MSG
    HYPERBEAM_LIB HYPERBEAM_INCLUDE_DIR CFITSIO_LIBRARY M_LIBRARY)

endif(NOT HYPERBEAM_FOUND)

if(NOT HYPERBEAM_HDF5)
    find_file(HYPERBEAM_HDF5 mwa_full_embedded_element_pattern.h5
        HINTS ${HYPERBEAM_HDF5_DIR} ${HYPERBEAM_ROOT} ${HYPERBEAM_INCLUDE_DIR})
endif(NOT HYPERBEAM_HDF5)

if(HYPERBEAM_HDF5)
  message (STATUS "HYPERBEAM_HDF5 found: ${HYPERBEAM_HDF5}")
else(HYPERBEAM_HDF5)
    message (STATUS "Warning HYPERBEAM_HDF5 not set: only analytic beam will be available")
endif(HYPERBEAM_HDF5)