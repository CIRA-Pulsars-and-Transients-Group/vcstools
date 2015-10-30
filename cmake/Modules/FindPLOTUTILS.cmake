# - Try to find PLOTUTILS code.
# Variables used by this module:
#  PLOTUTILS_ROOT_DIR     - PLOTUTILS root directory
# Variables defined by this module:
#  PLOTUTILS_FOUND        - system has PLOTUTILS
#  PLOTUTILS_INCLUDE_DIR  - the PLOTUTILS include directory (cached)
#  PLOTUTILS_INCLUDE_DIRS - the PLOTUTILS include directories
#                         (identical to PLOTUTILS_INCLUDE_DIR)
#  PLOTUTILS_LIBRARY      - the PLOTUTILS library (cached)
#  PLOTUTILS_LIBRARIES    - the PLOTUTILS libraries
#  PGPLOT_LIBRARY	- the <FORTRAN> PGPLOT library

set(PLOTUTILS_ROOT_DIR $ENV{PLOTUTILS_DIR})

if(NOT DEFINED PLOTUTILS_ROOT_DIR)
	message("-- Warning PLOTUTILS_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED PLOTUTILS_ROOT_DIR)
	message("-- PLOTUTILS_ROOT_DIR = ${PLOTUTILS_ROOT_DIR}")
endif(NOT DEFINED PLOTUTILS_ROOT_DIR)

if(NOT PLOTUTILS_FOUND)

  find_path(PLOTUTILS_INCLUDE_DIR plot.h
    HINTS ${PLOTUTILS_ROOT_DIR} PATH_SUFFIXES include include/pgplot)
  find_library(PLOTUTILS_LIBRARY plot
    HINTS ${PLOTUTILS_ROOT_DIR} PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PLOTUTILS DEFAULT_MSG
    PLOTUTILS_LIBRARY PLOTUTILS_INCLUDE_DIR )

  set(PLOTUTILS_INCLUDE_DIRS ${PLOTUTILS_INCLUDE_DIR})
  set(PLOTUTILS_LIBRARIES ${PLOTUTILS_LIBRARY})

endif(NOT PLOTUTILS_FOUND)
if (PLOTUTILS_FOUND)
	message("-- Found PlotUtils --")
endif(PLOTUTILS_FOUND)
