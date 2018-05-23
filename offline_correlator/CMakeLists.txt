cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (offline_correlator LANGUAGES C CUDA)

file(GLOB sources "*.cu" "*.c")
add_executable(offline_correlator ${sources})
set_target_properties(offline_correlator PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA")
include_directories(${CFITSIO_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} ${XGPU_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/../make_beam /usr/include)
target_link_libraries(offline_correlator ${CFITSIO_LIBRARY} ${M_LIBRARY} ${XGPU_LIBRARY})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/offline_correlator DESTINATION bin)

