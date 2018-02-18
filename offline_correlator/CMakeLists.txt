cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (offline_correlator LANGUAGES C CUDA)

include_directories(${CFITSIO_INCLUDE_DIR} ${CUDA_INCLUDE_DIRS} ${XGPU_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/../make_beam_small /usr/include)
file(GLOB sources "*.cc" "*.cu" "*.c")
add_executable(offline_correlator ${sources})
target_link_libraries(offline_correlator ${CFITSIO_LIBRARY} ${M_LIBRARY} ${XGPU_LIBRARY})
set_target_properties(offline_correlator PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DHAVE_CUDA")

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/offline_correlator DESTINATION bin)

