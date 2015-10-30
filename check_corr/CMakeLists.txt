
add_executable(check_corr check_corr)
find_package(CPGPLOT)
find_package(FFTW3 COMPONENTS single double)
find_package(PLOTUTILS)
find_package(CFITSIO)
include_directories(${CMAKE_INSTALL_PREFIX}/include 
${TWOPIP_INCLUDE_DIR} ${CPGPLOT_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} 
${PLOTUTILS_INCLUDE_DIR} ${CFITSIO_INCLUDE_DIR} ${TWOPIP_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/mwac_utils)
target_link_libraries(check_corr MWACUtils ${TWOPIP_LIBRARIES}
${CPGPLOT_LIBRARIES} ${FFTW3F_LIBRARY} ${PLOTUTILS_LIBRARY} ${CFITSIO_LIBRARY} -L/opt/X11/lib -lX11 -lXext -lXt -lXaw -lpng)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/check_corr DESTINATION bin)

