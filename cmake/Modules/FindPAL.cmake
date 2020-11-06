# - Try to find PAL code.
# Variables used by this module:
#  PAL_ROOT_DIR     - PAL root directory
# Variables defined by this module:
#  PAL_FOUND        - system has PAL
#  PAL_INCLUDE_DIR  - the PAL include directory (cached)
#  PAL_INCLUDE_DIRS - the PAL include directories
#                         (identical to PAL_INCLUDE_DIR)
#  PAL_LIBRARY      - the PAL library (cached)
#  PAL_LIBRARIES    - the PAL libraries
#                         (identical to PAL_LIBRARY)

message("Finding PAL")

set(PAL_ROOT_DIR $ENV{PAL_DIR})

if(NOT DEFINED PAL_ROOT_DIR)
	message(STATUS "Warning PAL_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED PAL_ROOT_DIR)
	message(STATUS "PAL_ROOT_DIR = ${PAL_ROOT_DIR}")
endif(NOT DEFINED PAL_ROOT_DIR)

if(NOT PAL_FOUND)

#  find_path(PAL_INCLUDE_DIR pal.h
#      HINTS ${PAL_ROOT_DIR} PATH_SUFFIXES include NO_SYSTEM_ENVIRONMENT_PATH)
#  find_library(PAL_LIBRARY pal
#      HINTS ${PAL_ROOT_DIR} PATH_SUFFIXES lib)
#  mark_as_advanced(PAL_INCLUDE_DIR PAL_LIBRARY NO_SYSTEM_ENVIRONMENT_PATH)

  set(PAL_INCLUDE_DIR ${PAL_ROOT_DIR}/include)
  set(PAL_LIBRARY ${PAL_ROOT_DIR}/lib/libpal.so)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PAL DEFAULT_MSG
    PAL_LIBRARY PAL_INCLUDE_DIR)

  set(PAL_INCLUDE_DIRS ${PAL_INCLUDE_DIR})
  set(PAL_LIBRARIES ${PAL_LIBRARY} ${M_LIBRARY})

endif(NOT PAL_FOUND)

if (PAL_FOUND)
    message(STATUS "Found PAL (${PAL_LIBRARIES})")
endif (PAL_FOUND)


