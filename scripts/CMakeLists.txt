cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)


project(scripts)

set(SCRIPT_TARGETS "${PROJECT_SOURCE_DIR}/bf_adjust_flags.py;${PROJECT_SOURCE_DIR}/checks.py;${PROJECT_SOURCE_DIR}/calibrate_vcs.py;${PROJECT_SOURCE_DIR}/create_psrfits.sh;${PROJECT_SOURCE_DIR}/file_maxmin.py;${PROJECT_SOURCE_DIR}/find_pulsar_in_obs.py;${PROJECT_SOURCE_DIR}/fits_db.py;${PROJECT_SOURCE_DIR}/fix_offset.py;${PROJECT_SOURCE_DIR}/get_meta.py;${PROJECT_SOURCE_DIR}/get_voltage_obs.py;${PROJECT_SOURCE_DIR}/hdr.fits;${PROJECT_SOURCE_DIR}/obs_query.py;${PROJECT_SOURCE_DIR}/process_vcs.py;${PROJECT_SOURCE_DIR}/recombine.py;${PROJECT_SOURCE_DIR}/rename_corr_output.py;${PROJECT_SOURCE_DIR}/reorder_chans.py;${PROJECT_SOURCE_DIR}/rts2ao.py;${PROJECT_SOURCE_DIR}/run_rts.sh;${PROJECT_SOURCE_DIR}/untar.sh;${PROJECT_SOURCE_DIR}/cleanup.py;${PROJECT_SOURCE_DIR}/create_ics_psrfits.py;${PROJECT_SOURCE_DIR}/splice.sh")

file(WRITE ${CMAKE_INSTALL_PREFIX}/bin/version.py "__version__ = \"${VCSTOOLS_VERSION}\"")

install(PROGRAMS ${SCRIPT_TARGETS} DESTINATION bin)

