# - Try to find CPGPLOT code.
# Variables used by this module:
#  CPGPLOT_ROOT_DIR     - CPGPLOT root directory
# Variables defined by this module:
#  CPGPLOT_FOUND        - system has CPGPLOT
#  CPGPLOT_INCLUDE_DIR  - the CPGPLOT include directory (cached)
#  CPGPLOT_INCLUDE_DIRS - the CPGPLOT include directories
#                         (identical to CPGPLOT_INCLUDE_DIR)
#  CPGPLOT_LIBRARY      - the CPGPLOT library (cached)
#  CPGPLOT_LIBRARIES    - the CPGPLOT libraries
#  PGPLOT_LIBRARY	- the <FORTRAN> PGPLOT library

set(CPGPLOT_ROOT_DIR $ENV{PGPLOT_DIR})

if(NOT DEFINED CPGPLOT_ROOT_DIR)
	message("-- Warning CPGPLOT_ROOT_DIR not set: will try and find it ")
else(NOT DEFINED CPGPLOT_ROOT_DIR)
	message("-- CPGPLOT_ROOT_DIR = ${CPGPLOT_ROOT_DIR}")
endif(NOT DEFINED CPGPLOT_ROOT_DIR)

if(NOT CPGPLOT_FOUND)

  find_path(CPGPLOT_INCLUDE_DIR cpgplot.h
    HINTS ${CPGPLOT_ROOT_DIR} PATH_SUFFIXES include include/pgplot)
  find_library(CPGPLOT_LIBRARY cpgplot
    HINTS ${CPGPLOT_ROOT_DIR} PATH_SUFFIXES lib)
  find_library(PGPLOT_LIBRARY pgplot
    HINTS ${CPGPLOT_ROOT_DIR} PATH_SUFFIXES lib)
  find_library(GFORTRAN_LIBRARY gfortran
    HINTS PATH_SUFFIXES gcc/x86_64-redhat-linux/4.4.4/ )
  find_library(X11_LIBRARY X11)
  find_library(M_LIBRARY m)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(CPGPLOT DEFAULT_MSG
    CPGPLOT_LIBRARY M_LIBRARY CPGPLOT_INCLUDE_DIR PGPLOT_LIBRARY GFORTRAN_LIBRARY X11_LIBRARY)

  set(CPGPLOT_INCLUDE_DIRS ${CPGPLOT_INCLUDE_DIR})
  set(CPGPLOT_LIBRARIES ${CPGPLOT_LIBRARY} ${PGPLOT_LIBRARY} ${GFORTRAN_LIBRARY} ${X11_LIBRARY} ${M_LIBRARY})

endif(NOT CPGPLOT_FOUND)

if (CPGPLOT_FOUND)
	message ("--Found cpgplot --")
endif (CPGPLOT_FOUND)

