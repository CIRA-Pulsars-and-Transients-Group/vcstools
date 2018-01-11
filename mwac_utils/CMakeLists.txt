cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project(mwac_utils)

include_directories(/usr/local/include ./ ${FFTW3_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIRS} ${PLOTUTILS_INCLUDE_DIR})

add_library(mwac_utils STATIC mwac_utils.c antenna_mapping.c ascii_header.c uvfits.c vdifio.c )
target_link_libraries(mwac_utils ${FFTW3F_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY})
set_target_properties(mwac_utils PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")

#install(FILES mwac_utils.h antenna_mapping.h uvfits.h vdifio.h ascii_header.h mwa_header.h DESTINATION include)
#install(TARGETS mwac_utils ARCHIVE DESTINATION lib)

