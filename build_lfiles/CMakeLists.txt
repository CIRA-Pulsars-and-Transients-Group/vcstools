
cmake_minimum_required (VERSION 2.6)
project (build_lfiles)
add_executable(build_lfiles build_lfiles.c)
add_executable(uvcompress uvcompress.cpp compress.cpp)
set_target_properties(build_lfiles PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")
set_target_properties(uvcompress PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")
find_package(CFITSIO)
find_package(PLOTUTILS)
find_package(FFTW3 COMPONENTS single)

include_directories(${CFITSIO_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/../mwac_utils ${CMAKE_INSTALL_PREFIX}/include)
target_link_libraries(build_lfiles -L${CMAKE_INSTALL_PREFIX}/lib MWACUtils ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} m)
target_link_libraries(uvcompress ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} m)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/build_lfiles DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/uvcompress DESTINATION bin)
