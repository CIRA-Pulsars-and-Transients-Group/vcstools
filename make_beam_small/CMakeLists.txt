cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (make_beam_small LANGUAGES C CUDA)

if(CUDA_FOUND)

    include_directories(${FFTW3_INCLUDE_DIR} ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS})
    #file(GLOB sources "*.c" "*.cu")
    set(sources "${PROJECT_SOURCE_DIR}/beam_common.c;${PROJECT_SOURCE_DIR}/beam_psrfits.c;${PROJECT_SOURCE_DIR}/beam_vdif.c;${PROJECT_SOURCE_DIR}/filter.c;${PROJECT_SOURCE_DIR}/form_beam.cu;${PROJECT_SOURCE_DIR}/get_delays_small.c;${PROJECT_SOURCE_DIR}/ipfb.cu;make_beam_small.c;${PROJECT_SOURCE_DIR}/ascii_header.c;${PROJECT_SOURCE_DIR}/vdifio.c")
    add_executable(make_beam_small ${sources})
    target_link_libraries(make_beam_small ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${FFTW3_OMP_LIBRARY} ${FFTW3_LIBRARY})
    set_target_properties(make_beam_small PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA")

else ()

    include_directories(${FFTW3_INCLUDE_DIR} ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR})
    file(GLOB sources "*.c")
    add_executable(make_beam_small ${sources})
    target_link_libraries(make_beam_small ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} ${FFTW3_OMP_LIBRARY} ${FFTW3_LIBRARY})

endif(CUDA_FOUND)    

#if(CUDA_FOUND)
#    if (CRAY_CUDATOOLKIT_DIR)
#        set_target_properties(make_beam_small PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA ${CRAY_CUDA_TOOLKIT_INCLUDE_OPTS}")
#    else ()
#        set_target_properties(make_beam_small PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -I/usr/local/cuda/samples/common/inc/")
#    endif ()    
#endif(CUDA_FOUND)    

#if (CUDA_FOUND)
#    include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CUDA_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX}/include ${SLALIB_INCLUDE_DIR})
#    
#else ()
#    include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CMAKE_INSTALL_PREFIX}/include ${SLALIB_INCLUDE_DIR})
#endif (CUDA_FOUND)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam_small DESTINATION bin)

