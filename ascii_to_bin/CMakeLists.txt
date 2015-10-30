
cmake_minimum_required (VERSION 2.6)
project (ascii_to_bin)
add_executable(ascii_to_bin ascii_to_bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ascii_to_bin DESTINATION bin)

