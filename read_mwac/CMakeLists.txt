cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (read_mwac)
add_executable(read_mwac read_mwac.c)
#set_target_properties(read_mwac PROPERTIES LINK_FLAGS "-Wl,-no_compact_unwind")
include_directories(${CFITSIO_INCLUDE_DIR} ${CMAKE_INSTALL_PREFIX}/include/ )
target_link_libraries(read_mwac -L${CMAKE_INSTALL_PREFIX}/lib MWAC_Utils ${CFITSIO_LIBRARY}  m)

#install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/read_mwac DESTINATION bin)
