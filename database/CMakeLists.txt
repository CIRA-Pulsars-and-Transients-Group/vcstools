cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)


project(database)

set(DATABASE_TARGETS "${PROJECT_SOURCE_DIR}/submit_to_database.py;${PROJECT_SOURCE_DIR}/database_vcs.py;${PROJECT_SOURCE_DIR}/known_repeating_FRBs.csv")

file(WRITE ${CMAKE_INSTALL_PREFIX}/bin/version.py "__version__ = \"${VCSTOOLS_VERSION}\"")

install(PROGRAMS ${DATABASE_TARGETS} DESTINATION bin)

