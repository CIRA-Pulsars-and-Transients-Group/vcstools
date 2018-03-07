# - Try to find XGPU.
# Variables used by this module:
#  XGPU_ROOT_DIR     - XGPU root directory
# Variables defined by this module:
#  XGPU_FOUND        - system has XGPU
#  XGPU_INCLUDE_DIR  - the XGPU include directory (cached)
#  XGPU_INCLUDE_DIRS - the XGPU include directories
#                         (identical to XGPU_INCLUDE_DIR)
#  XGPU_LIBRARY      - the XGPU library (cached)
#  XGPU_LIBRARIES    - the XGPU libraries
#                         (identical to XGPU_LIBRARY)

message("Finding XGPU")

if(NOT XGPU_FOUND)

  find_path(XGPU_INCLUDE_DIR xgpu.h
    HINTS ${XGPU_ROOT_DIR} PATH_SUFFIXES include /usr/local/include include/xgpu)
  find_library(XGPU_LIBRARY xgpu
    HINTS ${XGPU_ROOT_DIR} PATH_SUFFIXES lib/)
  mark_as_advanced(XGPU_INCLUDE_DIR XGPU_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(XGPU DEFAULT_MSG
    XGPU_LIBRARY XGPU_INCLUDE_DIR)

  set(XGPU_INCLUDE_DIRS ${XGPU_INCLUDE_DIR})
  set(XGPU_LIBRARIES ${XGPU_LIBRARY})

endif(NOT XGPU_FOUND)

if (XGPU_FOUND)
    message(STATUS "Found XGPU (${XGPU_LIBRARIES})")
endif (XGPU_FOUND)

