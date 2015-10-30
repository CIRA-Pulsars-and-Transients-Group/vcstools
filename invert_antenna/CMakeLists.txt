
cmake_minimum_required (VERSION 2.6)
project (invert_antenna)
add_executable(invert_antenna invert_antenna)
find_package(FFTW3 COMPONENTS single)
include_directories(${CMAKE_INSTALL_PREFIX}/include 
 ${FFTW3_INCLUDE_DIR})
target_link_libraries(invert_antenna  
${FFTW3F_LIBRARY} )

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/invert_antenna DESTINATION bin)

