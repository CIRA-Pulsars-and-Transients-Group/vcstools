cmake_minimum_required (VERSION 2.6)
project(mwac_utils)

find_package(PLOTUTILS)
find_package(CFITSIO)
find_package(SLALIB)
find_package(FFTW3 COMPONENTS single)

include_directories(/usr/local/include ${FFTW3_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIRS} ${PLOTUTILS_INCLUDE_DIR})

add_library(MWACUtils SHARED mwac_utils.c antenna_mapping.c ascii_header.c uvfits.c vdifio.c )
target_link_libraries(MWACUtils ${FFTW3F_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY})
set_target_properties(MWACUtils PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")

install(FILES mwac_utils.h antenna_mapping.h uvfits.h vdifio.h ascii_header.h mwa_header.h DESTINATION include)
install(TARGETS MWACUtils LIBRARY DESTINATION lib)

