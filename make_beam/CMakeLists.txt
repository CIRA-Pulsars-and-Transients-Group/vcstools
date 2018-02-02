
cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (make_beam)
add_executable(make_beam make_beam)

SET(BEAM_CXX_FLAGS "" )
SET(BEAM_C_FLAGS "" )
if(MPI_COMPILE_FLAGS) 
  set_target_properties(make_beam PROPERTIES 
    COMPILE_FLAGS "{CMAKE_C_FLAGS} ${BEAM_C_FLAGS} ${MPI_COMPILE_FLAGS}") 
endif() 

if(CUDA_FOUND)
    if (CRAY_CUDATOOLKIT_DIR)
        set_target_properties(make_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA ${CRAY_CUDA_TOOLKIT_INCLUDE_OPTS}")
    else ()
        set_target_properties(make_beam PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -I/usr/local/cuda/samples/common/inc/")
    endif ()    
endif(CUDA_FOUND)    

if(MPI_LINK_FLAGS) 
  set_target_properties(make_beam PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}") 
endif() 

if (CUDA_FOUND)
    include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CUDA_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX}/include)
    
else ()
    include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CMAKE_INSTALL_PREFIX}/include)
endif (CUDA_FOUND)

if (MPI_C_LIBRARIES)
    if (CUDA_FOUND)
        target_link_libraries(make_beam mwac_utils MWAC_GPU_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_C_LIBRARIES} m)
    else ()
        target_link_libraries(make_beam mwac_utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_C_LIBRARIES} m)
    endif (CUDA_FOUND)
else()
    if (CUDA_FOUND)
        target_link_libraries(make_beam mwac_utils MWAC_GPU_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_LIBRARIES} m)
    else ()
        target_link_libraries(make_beam mwac_utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_LIBRARIES} m)
    endif (CUDA_FOUND)
endif()

if (CUDA_FOUND)
    #    set_target_properties(make_beam PROPERTIES LINK_FLAGS "-Wl,-rpath,/Developer/NVIDIA/CUDA-7.5/lib") 
    target_link_libraries(make_beam ${CUDA_CUFFT_LIBRARIES} m)
    target_link_libraries(make_beam ${CUDA_LIBRARIES} m)
endif ()

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (MPI_EXTRA_LIBRARY)
    target_link_libraries(make_beam ${MPI_EXTRA_LIBRARY} m)
endif()

add_definitions(${MPI_C_COMPILE_FLAGS})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam DESTINATION bin)

