cmake_minimum_required (VERSION 2.6)
project(mwac_utils)

include_directories(/usr/local/include ./ ${FFTW3_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIRS} ${PLOTUTILS_INCLUDE_DIR})

add_library(MWAC_Utils STATIC mwac_utils.c antenna_mapping.c ascii_header.c uvfits.c vdifio.c )
target_link_libraries(MWAC_Utils ${FFTW3F_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY})
set_target_properties(MWAC_Utils PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")

#install(FILES mwac_utils.h antenna_mapping.h uvfits.h vdifio.h ascii_header.h mwa_header.h DESTINATION include)
#install(TARGETS MWAC_Utils ARCHIVE DESTINATION lib)

