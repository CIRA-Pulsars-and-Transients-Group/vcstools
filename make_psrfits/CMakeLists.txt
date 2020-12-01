cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)

project (make_psrfits)

add_executable(make_psrfits make_psrfits)
add_executable(splice_psrfits splice_psrfits)

include_directories(${CMAKE_INSTALL_PREFIX}/include ${PAL_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${PSRFITS_UTILS_INCLUDE_DIR})

target_link_libraries(make_psrfits ${PSRFITS_UTILS_LIBRARY} ${PAL_LIBRARY} ${CFITSIO_LIBRARY} ${M_LIBRARY})
target_link_libraries(splice_psrfits ${PSRFITS_UTILS_LIBRARY} ${PAL_LIBRARY} ${CFITSIO_LIBRARY} ${M_LIBRARY})

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/make_psrfits DESTINATION bin)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/splice_psrfits DESTINATION bin)

