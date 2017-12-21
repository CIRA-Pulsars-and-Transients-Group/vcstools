cmake_minimum_required (VERSION 2.6)
project (make_beam_small)

include_directories(${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${SLALIB_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/../mwac_utils)

file(GLOB sources "*.c")
add_executable(make_beam_small ${sources})
target_link_libraries(make_beam_small ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${SLALIB_LIBRARY} ${M_LIBRARY} mwac_utils)



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

#install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam_small DESTINATION bin)

