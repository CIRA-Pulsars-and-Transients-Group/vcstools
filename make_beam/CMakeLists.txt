
cmake_minimum_required (VERSION 2.6)
project (make_beam)
add_executable(make_beam make_beam)

find_package(CPGPLOT)
find_package(PSRFITS_UTILS)
find_package(CFITSIO)
find_package(MPI REQUIRED)
find_package(FFTW3 COMPONENTS single)
find_package(CUDA)

#set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_C_FLAGS "-O3 ")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_RELEASE "-O3 ")
set(CMAKE_CXX_FLAGS "-O3 ")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 ")


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

include_directories({MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/../mwac_utils ${PROJECT_SOURCE_DIR}/../gpu_utils ${PSRFITS_UTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${MPI_INCLUDE_PATH} ${CUDA_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX}/include)
if (MPI_C_LIBRARIES)
    target_link_libraries(make_beam MWAC_Utils MWAC_GPU_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_C_LIBRARIES})
else()
    target_link_libraries(make_beam MWAC_Utils MWAC_GPU_Utils ${PSRFITS_UTILS_LIBRARY} ${CFITSIO_LIBRARY} ${FFTW3F_LIBRARY} ${MPI_LIBRARIES})
endif()

if (CUDA_FOUND)
    #    set_target_properties(make_beam PROPERTIES LINK_FLAGS "-Wl,-rpath,/Developer/NVIDIA/CUDA-7.5/lib") 
    target_link_libraries(make_beam ${CUDA_CUFFT_LIBRARIES})
    target_link_libraries(make_beam ${CUDA_LIBRARIES})
endif ()


message ("MPI LIBS ${MPI_LIBRARIES}")
message ("MPI COMPILER ${MPI_COMPILER}")
message ("CUDA ${CUDA_CUFFT_LIBRARIES}")

if (MPI_EXTRA_LIBRARY)
    target_link_libraries(make_beam ${MPI_EXTRA_LIBRARY})
endif()

add_definitions(${MPI_C_COMPILE_FLAGS})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_beam DESTINATION bin)

