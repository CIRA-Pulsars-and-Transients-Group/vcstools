cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)


project(utils)

set(UTILITY_TARGETS "${PROJECT_SOURCE_DIR}/zapchan.py;${PROJECT_SOURCE_DIR}/calc_ephem.py;${PROJECT_SOURCE_DIR}/check_disk_usage.sh;${PROJECT_SOURCE_DIR}/check_quota.sh;${PROJECT_SOURCE_DIR}/mdir.py;${PROJECT_SOURCE_DIR}/mwa_metadb_utils.py;${PROJECT_SOURCE_DIR}/job_submit.py;${PROJECT_SOURCE_DIR}/plotFlatTileBeam.py;${PROJECT_SOURCE_DIR}/plotPolarTileBeam.py;${PROJECT_SOURCE_DIR}/plotTiedArrayBeam.py;${PROJECT_SOURCE_DIR}/get_MWAobsInfo.py;${PROJECT_SOURCE_DIR}/aocal.py")

file(WRITE ${CMAKE_INSTALL_PREFIX}/bin/version.py "__version__ = \"${VCSTOOLS_VERSION}\"")

install(PROGRAMS ${UTILITY_TARGETS} DESTINATION bin)

