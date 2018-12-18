cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (make_beam LANGUAGES C CUDA)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    message("No OpenMP")
endif()

if(CUDA_FOUND)

    #file(GLOB sources "*.c" "*.cu")
    set(sources "${PROJECT_SOURCE_DIR}/beam_common.c;${PROJECT_SOURCE_DIR}/beam_psrfits.c;${PROJECT_SOURCE_DIR}/beam_vdif.c;${PROJECT_SOURCE_DIR}/filter.c;${PROJECT_SOURCE_DIR}/form_beam.cu;${PROJECT_SOURCE_DIR}/get_delays_small.c;${PROJECT_SOURCE_DIR}/ipfb.cu;make_beam.c;${PROJECT_SOURCE_DIR}/ascii_header.c;${PROJECT_SOURCE_DIR}/vdifio.c")
    add_executable(make_beam ${sources})
    set_target_properties(make_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -lineinfo")
    include_directories(${FFTW3_INCLUDE_DIR} ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
    target_link_libraries(make_beam -fopenmp ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY})

else ()

    file(GLOB sources "*.c")
    add_executable(make_beam ${sources})
    include_directories(${FFTW3_INCLUDE_DIR} ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR})
    target_link_libraries(make_beam ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${FFTW3_OMP_LIBRARY} ${FFTW3_LIBRARY})

endif(CUDA_FOUND)    

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam DESTINATION bin)

