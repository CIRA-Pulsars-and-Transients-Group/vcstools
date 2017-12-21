
cmake_minimum_required (VERSION 2.6)
project (make_beam_small)
add_executable(make_beam_small make_beam_small get_delays_small beam_common beam_psrfits beam_vdif)
#target_link_libraries(get_delays_small -L${CMAKE_INSTALL_PREFIX}/lib MWAC_Utils)
#target_link_libraries(get_delays_small ${PSRFITS_UTILS_LIBRARY} )
#target_link_libraries(get_delays_small ${SLALIB_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} -lm)

find_package(CPGPLOT)
find_package(SLALIB)
find_package(PSRFITS_UTILS)
find_package(CFITSIO)
find_package(MPI REQUIRED)
find_package(FFTW3 COMPONENTS single)
find_package(CUDA)
#find_package(OpenMP)


SET(BEAM_CXX_FLAGS "" )
SET(BEAM_C_FLAGS "" )
if(MPI_COMPILE_FLAGS) 
  set_target_properties(make_beam_small PROPERTIES 
    COMPILE_FLAGS "{CMAKE_C_FLAGS} ${BEAM_C_FLAGS} ${MPI_COMPILE_FLAGS}") 
endif() 

if(CUDA_FOUND)
    if (CRAY_CUDATOOLKIT_DIR)
        set_target_properties(make_beam_small PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA ${CRAY_CUDA_TOOLKIT_INCLUDE_OPTS}")
    else ()
        set_target_properties(make_beam_small PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA -I/usr/local/cuda/samples/common/inc/")
    endif ()    
endif(CUDA_FOUND)    

if(MPI_LINK_FLAGS) 
  set_target_properties(make_beam_small PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}") 
endif() 

if (CUDA_FOUND)
    include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CUDA_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX}/include ${SLALIB_INCLUDE_DIR})
    
else ()
    include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CMAKE_INSTALL_PREFIX}/include ${SLALIB_INCLUDE_DIR})
endif (CUDA_FOUND)

if (MPI_C_LIBRARIES)
    if (CUDA_FOUND)
        target_link_libraries(make_beam_small MWAC_Utils MWAC_GPU_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_C_LIBRARIES} m)
    else ()
        target_link_libraries(make_beam_small MWAC_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_C_LIBRARIES} m)
    endif (CUDA_FOUND)
else()
    if (CUDA_FOUND)
        target_link_libraries(make_beam_small MWAC_Utils MWAC_GPU_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_LIBRARIES} m)
    else ()
        target_link_libraries(make_beam_small MWAC_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_LIBRARIES} m)
    endif (CUDA_FOUND)
endif()

if (CUDA_FOUND)
    #    set_target_properties(make_beam_small PROPERTIES LINK_FLAGS "-Wl,-rpath,/Developer/NVIDIA/CUDA-7.5/lib") 
    target_link_libraries(make_beam_small ${CUDA_CUFFT_LIBRARIES} m)
    target_link_libraries(make_beam_small ${CUDA_LIBRARIES} m)
endif ()

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

message ("MPI LIBS ${MPI_LIBRARIES}")
message ("MPI COMPILER ${MPI_COMPILER}")
message ("CUDA ${CUDA_CUFFT_LIBRARIES}")

if (MPI_EXTRA_LIBRARY)
    target_link_libraries(make_beam_small ${MPI_EXTRA_LIBRARY} m)
endif()

add_definitions(${MPI_C_COMPILE_FLAGS})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam_small DESTINATION bin)

