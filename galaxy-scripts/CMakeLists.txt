cmake_minimum_required (VERSION 2.6)
cmake_policy(SET CMP0012 NEW)
cmake_policy(SET CMP0048 NEW)


project(galaxy-scripts)

set(s "${PROJECT_SOURCE_DIR}/scripts")
set(d "${PROJECT_SOURCE_DIR}/database")
set(u "${PROJECT_SOURCE_DIR}/utils")

set(SCRIPT_TARGETS "${s}/bf_adjust_flags.py;${s}/checks.py;${s}/calibrate_vcs.py;${s}/create_psrfits.sh;${s}/file_maxmin.py;${s}/find_pulsar_in_obs.py;${s}/fits_db.py;${s}/fix_offset.py;${s}/get_meta.py;${s}/get_voltage_obs.py;${s}/hdr.fits;${s}/obs_query.py;${s}/process_vcs.py;${s}/recombine.py;${s}/rename_corr_output.py;${s}/reorder_chans.py;${s}/rts2ao.py;${s}/run_rts.sh;${s}/untar.sh")

set(DATABASE_TARGETS "${d}/submit_to_database.py;${d}/database_vcs.py")

set(UTILITY_TARGETS "${u}/zapchan.py;${u}/calc_ephem.py;${u}/check_disk_usage.sh;${u}/check_quota.sh;${u}/plot_tilebeams.py;${u}/mdir.py;${u}/mwa_metadb_utils.py")


install(PROGRAMS ${SCRIPT_TARGETS} DESTINATION bin)
install(PROGRAMS ${DATABASE_TARGETS} DESTINATION bin)
install(PROGRAMS ${UTILITY_TARGETS} DESTINATION bin)

