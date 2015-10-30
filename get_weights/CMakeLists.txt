
cmake_minimum_required (VERSION 2.6)
project (get_weights)
add_executable(get_weights get_weights)
include_directories(${CMAKE_INSTALL_PREFIX}/include)
target_link_libraries(get_weights -lm)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/get_weights DESTINATION bin)

