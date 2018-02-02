cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (tests)

include_directories(${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../make_beam_small)

#file(GLOB sources "*.c")
set(sources "tests.c;${PROJECT_SOURCE_DIR}/../make_beam_small/beam_vdif.c;${PROJECT_SOURCE_DIR}/../make_beam_small/beam_common.c;${PROJECT_SOURCE_DIR}/../make_beam_small/filter.c")

add_executable(tests ${sources})
target_link_libraries(tests ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${FFTW3_OMP_LIBRARY} ${FFTW3_LIBRARY} mwac_utils)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/tests DESTINATION ${PROJECT_SOURCE_DIR})
