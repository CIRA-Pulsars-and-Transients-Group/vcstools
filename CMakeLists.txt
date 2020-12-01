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


if (${CMAKE_BEAMFORMER})
  message("Finding CUDA")
  find_package(CUDA)
endif()

message("Finding OpenMP")
include(FindOpenMP)

find_package(PAL)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
find_package(FFTW3 COMPONENTS single openmp)
find_package(XGPU)
find_package(HYPERBEAM)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fPIC -Wall -Wextra ${OpenMP_C_FLAGS}")

EXECUTE_PROCESS( COMMAND gcc -dumpversion GCC_VERSION )
if(GCC_VERSION EQUAL 7.30)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
endif()

# Add the subdirectories for the three target projects
if (XGPU_FOUND)
    add_subdirectory(offline_correlator)
endif()

if (${CMAKE_BEAMFORMER})
  message("Compiling make_beam")
  add_subdirectory(make_beam)
else()
  message("Not compiling make_beam")
endif()

if (${CTESTS})
  message("Compiling C tests")
  add_subdirectory(tests)
else()
    message("Not compiling C tests")
endif()

add_subdirectory(make_psrfits)
add_subdirectory(read_mwac)

INCLUDE(InstallRequiredSystemLibraries)
