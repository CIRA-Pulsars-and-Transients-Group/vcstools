cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project(vcstools)

# Set up version number
execute_process(
  COMMAND bash -c "git describe --tags --long | sed 's/-/./' | sed 's/-/_/'"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE VCSTOOLS_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

#set(VCSTOOLS_VERSION "${PROJECT_VERSION}_${VERSION_GIT}")
message($VCSTOOLS_VERSION)

add_definitions( -DVERSION_BEAMFORMER="${VCSTOOLS_VERSION}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message("Finding CUDA")
find_package(CUDA)

message("Finding OpenMP")
include(FindOpenMP)

find_package(SLALIB)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
find_package(FFTW3 COMPONENTS single openmp)
find_package(XGPU)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -O2 -fPIC -Wall -Wextra ${OpenMP_C_FLAGS}")

# Add the subdirectories for the three target projects
if (XGPU_FOUND)
    add_subdirectory(offline_correlator)
endif()
add_subdirectory(make_beam)
add_subdirectory(make_psrfits)
add_subdirectory(read_mwac)

INCLUDE(InstallRequiredSystemLibraries)
