cmake_minimum_required(VERSION 2.6)
project(read_pfb)
find_package(FFTW3 COMPONENTS single)
add_executable(read_pfb read_pfb)
set_target_properties(read_pfb PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")

SET(PFB_CXX_FLAGS " " )
SET(PFB_C_FLAGS " " )
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PFB_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PFB_CXX_FLAGS}")

include_directories(${CMAKE_INSTALL_PREFIX}/include ${PROJECT_SOURCE_DIR}/mwac_utils)
target_link_libraries(read_pfb -L${CMAKE_INSTALL_PREFIX}/lib MWAC_Utils ${FFTW3F_LIBRARY})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/read_pfb DESTINATION bin)

