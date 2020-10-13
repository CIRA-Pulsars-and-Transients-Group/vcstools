cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (make_beam LANGUAGES C)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    message("No OpenMP")
endif()

set(sources "${PROJECT_SOURCE_DIR}/test_make_beam.c;${PROJECT_SOURCE_DIR}/../make_beam/get_delays_small.c;${PROJECT_SOURCE_DIR}/../make_beam/beam_common.c")
add_executable(test_make_beam ${sources})

set(sources "${PROJECT_SOURCE_DIR}/beam_model_jones_comparison.c;${PROJECT_SOURCE_DIR}/../make_beam/get_delays_small.c;${PROJECT_SOURCE_DIR}/../make_beam/beam_common.c")
add_executable(beam_model_jones_comparison ${sources})

set(sources "${PROJECT_SOURCE_DIR}/sky_model_jones_comparison.c;${PROJECT_SOURCE_DIR}/../make_beam/get_delays_small.c;${PROJECT_SOURCE_DIR}/../make_beam/beam_common.c")
add_executable(sky_model_jones_comparison ${sources})

if(CUDA_FOUND)
    set_target_properties(test_make_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -DHYPERBEAM_HDF5=\\\"${HYPERBEAM_HDF5}\\\"")
    set_target_properties(beam_model_jones_comparison PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -DHYPERBEAM_HDF5=\\\"${HYPERBEAM_HDF5}\\\"")
endif(CUDA_FOUND)

set_target_properties(sky_model_jones_comparison PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHYPERBEAM_HDF5=\\\"${HYPERBEAM_HDF5}\\\"")
include_directories(${FFTW3_INCLUDE_DIR} ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} ${HYPERBEAM_INCLUDE_DIR})

target_link_libraries(test_make_beam -fopenmp ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${HYPERBEAM_LIB})
target_link_libraries(beam_model_jones_comparison -fopenmp ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${HYPERBEAM_LIB})
target_link_libraries(sky_model_jones_comparison -fopenmp ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${HYPERBEAM_LIB})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/test_make_beam DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/beam_model_jones_comparison DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/sky_model_jones_comparison DESTINATION bin)

