cmake_minimum_required (VERSION 2.6)

cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (MWA_Beamformer VERSION 1.2.0)

add_definitions( -DVERSION_BEAMFORMER="${PROJECT_VERSION}" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(SLALIB)
find_package(PSRFITS_UTILS)
find_package(CFITSIO)
find_package(FFTW3 COMPONENTS single)
#find_package(Postgres)
#find_package(CUDA)

message("Finding OpenMP")
include(FindOpenMP)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -g -fPIC -Wall -Wextra ${OpenMP_C_FLAGS}")

add_subdirectory(mwac_utils)
add_subdirectory(read_mwac)
add_subdirectory(read_pfb)
add_subdirectory(make_psrfits)
add_subdirectory(get_delays)
#add_subdirectory(gpu_utils)    
add_subdirectory(make_beam)
add_subdirectory(make_beam_small)

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

