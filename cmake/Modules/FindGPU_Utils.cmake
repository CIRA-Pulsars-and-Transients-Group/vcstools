# - Try to find GPU_UTILS code.
# Variables used by this module:
#  GPU_UTILS_ROOT_DIR     - GPU_UTILS root directory
# Variables defined by this module:
#  GPU_UTILS_FOUND        - system has GPU_UTILS
#  GPU_UTILS_INCLUDE_DIR  - the GPU_UTILS include directory (cached)
#  GPU_UTILS_INCLUDE_DIRS - the GPU_UTILS include directories
#                         (identical to GPU_UTILS_INCLUDE_DIR)
#  GPU_UTILS_LIBRARY      - the GPU_UTILS library (cached)
#  GPU_UTILS_LIBRARIES    - the GPU_UTILS libraries
#                         (identical to GPU_UTILS_LIBRARY)

set(GPU_UTILS_ROOT_DIR $ENV{GPU_UTILS_DIR})

if(NOT DEFINED GPU_UTILS_ROOT_DIR)
	message("-- Warning GPU_UTILS_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED GPU_UTILS_ROOT_DIR)
	message("-- GPU_UTILS_ROOT_DIR = ${GPU_UTILS_ROOT_DIR}")
endif(NOT DEFINED GPU_UTILS_ROOT_DIR)

if(NOT GPU_UTILS_FOUND)

  find_path(GPU_UTILS_INCLUDE_DIR gpu_utils.h
    HINTS ${GPU_UTILS_ROOT_DIR}/include PATH_SUFFIXES include include/gpu_utils)
  find_library(GPU_UTILS_LIBRARY gpu_utils
    HINTS ${GPU_UTILS_ROOT_DIR}/trunk/lib PATH_SUFFIXES lib)
  find_library(M_LIBRARY m)
  mark_as_advanced(GPU_UTILS_INCLUDE_DIR GPU_UTILS_LIBRARY M_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(GPU_UTILS DEFAULT_MSG
    GPU_UTILS_LIBRARY M_LIBRARY GPU_UTILS_INCLUDE_DIR)

  set(GPU_UTILS_INCLUDE_DIRS ${GPU_UTILS_INCLUDE_DIR})
  set(GPU_UTILS_LIBRARIES ${GPU_UTILS_LIBRARY} ${M_LIBRARY})

endif(NOT GPU_UTILS_FOUND)

if (GPU_UTILS_FOUND)
	message("-- Found GPU_UTILS --")
endif (GPU_UTILS_FOUND)

