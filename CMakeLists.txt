cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project(vcstools VERSION 1.4.1)

# Set up version number
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VERSION_GIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(VCSTOOLS_VERSION "${PROJECT_VERSION}_${VERSION_GIT}")

add_definitions( -DVERSION_BEAMFORMER="${VCSTOOLS_VERSION}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message("Finding CUDA")
find_package(CUDA)
if (CUDA_FOUND)
    if(DEFINED ENV{PAWSEY_CLUSTER})
        if("$ENV{PAWSEY_CLUSTER}" STREQUAL "athena")
            set(CUDA_ARCH "sm_60")
        else
            set(CUDA_ARCH "sm_35")
        endif()
    endif()
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=${CUDA_ARCH}")
endif()

message("Finding OpenMP")
include(FindOpenMP)

find_package(SLALIB)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
find_package(FFTW3 COMPONENTS single openmp)
find_package(XGPU)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g -fPIC -Wall -Wextra ${OpenMP_C_FLAGS}")

# Add the subdirectories for the three target projects
if (XGPU_FOUND)
    add_subdirectory(offline_correlator)
endif()
add_subdirectory(make_beam)
add_subdirectory(make_psrfits)
add_subdirectory(galaxy-scripts)

INCLUDE(InstallRequiredSystemLibraries)
