
cmake_minimum_required (VERSION 2.6)
project (MWA_Utilities)

set(CMAKE_C_FLAGS "-O3 -g -fPIC")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -fPIC")

set (MWA_Utilities_MAJOR_VERSION 1)
set (MWA_Utilities_MINOR_VERSION 1)
set (MWA_Utilities_PATCH_VERSION 0)
set (MWA_Utilities_VERSION
    ${MWA_Utilities_MAJOR_VERSION}.${MWA_Utilities_MINOR_VERSION}.${MWA_Utilities_PATCH_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(USE_PROJECT_CMAKE_MODULE_PATH "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}")
set(USE_PROJECT_INSTALL_PREFIX "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")

if (APPLE)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BEAM_C_FLAGS} ${MPI_C_COMPILE_FLAGS} ")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BEAM_CXX_FLAGS} ")
else ()
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BEAM_C_FLAGS} ${MPI_C_COMPILE_FLAGS} -fopenmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BEAM_CXX_FLAGS} -fopenmp")
endif ()

find_package(CFITSIO)
find_package(PLOTUTILS)

if (CFITSIO_FOUND AND PLOTUTILS_FOUND)    
    add_subdirectory(mwac_utils)
endif()

find_package(PSRFITS_UTILS)

find_package(SLALIB)
find_package(Postgres)
find_package(MPI)

if (CFITSIO_FOUND)
    add_subdirectory(read_mwac)
    add_subdirectory(read_pfb)
    if (PSRFITS_UTILS_FOUND)  
        add_subdirectory(make_psrfits)
    endif (PSRFITS_UTILS_FOUND)
endif (CFITSIO_FOUND)

if (SLALIB_FOUND)
    if (PSRFITS_UTILS_FOUND) 
        add_subdirectory(get_delays)
    endif (PSRFITS_UTILS_FOUND)
endif (SLALIB_FOUND)


find_package(CUDA)
if (CUDA_FOUND)
    add_subdirectory(gpu_utils)    
endif (CUDA_FOUND)

if (PSRFITS_UTILS_FOUND)
    add_subdirectory(make_beam)
    add_subdirectory(make_beam_small)
endif (PSRFITS_UTILS_FOUND)


INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MWA Correlator Utilities")
SET(CPACK_PACKAGE_VENDOR "Stephen M. Ord")

SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake
${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
SET(CPACK_PACKAGE_EXECUTABLES "read_mwac" "make_beam" "read_pfb" "get_delays" "make_beam_small")
SET(CPACK_GENERATOR "RPM")
SET(CPACK_SET_DESTDIR ON)

INCLUDE(CPack)

