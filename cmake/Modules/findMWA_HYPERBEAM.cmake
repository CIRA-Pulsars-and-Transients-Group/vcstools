# - Try to find MWA_HYPERBEAM code.
# Variables used by this module:
#  MWA_HYPERBEAM_ROOT_DIR     - MWA_HYPERBEAM root directory
# Variables defined by this module:
#  MWA_HYPERBEAM_FOUND        - system has MWA_HYPERBEAM
#  MWA_HYPERBEAM_INCLUDE_DIR  - the MWA_HYPERBEAM include directory (cached)
#  MWA_HYPERBEAM_INCLUDE_DIRS - the MWA_HYPERBEAM include directories
#                         (identical to MWA_HYPERBEAM_INCLUDE_DIR)
#  MWA_HYPERBEAM_LIBRARY      - the MWA_HYPERBEAM library (cached)
#  MWA_HYPERBEAM_LIBRARIES    - the MWA_HYPERBEAM libraries

message("Finding MWA_HYPERBEAM")

set(MWA_HYPERBEAM_ROOT_DIR $ENV{MWA_HYPERBEAM})

if(NOT DEFINED MWA_HYPERBEAM_ROOT_DIR)
    message(STATUS "Warning MWA_HYPERBEAM_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED MWA_HYPERBEAM_ROOT_DIR)
    message(STATUS "MWA_HYPERBEAM_ROOT_DIR = ${MWA_HYPERBEAM_ROOT_DIR}")
endif(NOT DEFINED MWA_HYPERBEAM_ROOT_DIR)

if(NOT MWA_HYPERBEAM_FOUND)

  find_path(MWA_HYPERBEAM_INCLUDE_DIR mwa_hyperbeam.h
    HINTS ${MWA_HYPERBEAM_ROOT_DIR} PATH_SUFFIXES /include /include/mwa_hyperbeam)
  find_library(MWA_HYPERBEAM_LIBRARY mwa_hyperbeam
    HINTS ${MWA_HYPERBEAM_ROOT_DIR} PATH_SUFFIXES lib )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MWA_HYPERBEAM DEFAULT_MSG
    MWA_HYPERBEAM_LIBRARY MWA_HYPERBEAM_INCLUDE_DIR CFITSIO_LIBRARY M_LIBRARY)

  set(MWA_HYPERBEAM_INCLUDE_DIRS ${MWA_HYPERBEAM_INCLUDE_DIR})
  set(MWA_HYPERBEAM_LIBRARIES ${MWA_HYPERBEAM_LIBRARY}) #May rquire hdf5 lib

endif(NOT MWA_HYPERBEAM_FOUND)

if (MWA_HYPERBEAM_FOUND)
    message (STATUS "Found MWA_HYPERBEAM (${MWA_HYPERBEAM_LIBRARIES})")
endif (MWA_HYPERBEAM_FOUND)
