cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project(MWA_Beamformer VERSION 1.3.2)

add_definitions( -DVERSION_BEAMFORMER="${PROJECT_VERSION}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

message("Finding CUDA")
find_package(CUDA)
if (CUDA_FOUND)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_35")
endif()

message("Finding OpenMP")
include(FindOpenMP)

find_package(SLALIB)
find_package(CFITSIO)
find_package(PSRFITS_UTILS)
#find_package(FFTW3 COMPONENTS single openmp)
find_package(XGPU)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g -fPIC -Wall -Wextra ${OpenMP_C_FLAGS}")

# Add the subdirectories for the three target projects
if (XGPU_FOUND)
    add_subdirectory(offline_correlator)
endif()
add_subdirectory(make_beam_small)
add_subdirectory(make_psrfits)
add_subdirectory(galaxy-scripts)

## The following are deprecated:
#add_subdirectory(mwac_utils)
#add_subdirectory(read_mwac)
#add_subdirectory(read_pfb)
#add_subdirectory(get_delays)
#add_subdirectory(gpu_utils)    
#add_subdirectory(make_beam)

INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MWA Correlator Utilities")
SET(CPACK_PACKAGE_VENDOR "Stephen M. Ord")

SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
SET(CPACK_PACKAGE_EXECUTABLES "read_mwac" "make_beam" "read_pfb" "get_delays" "make_beam_small")
SET(CPACK_GENERATOR "RPM")
SET(CPACK_SET_DESTDIR ON)

INCLUDE(CPack)

