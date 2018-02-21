INSTALL_DIR = /group/mwaops/PULSAR/bin/

SCRIPT_TARGETS = bf_adjust_flags.py \
          build_mwacutils.sh \
		  calibrate_vcs.py \
          checks.py \
          create_psrfits.sh \
          fetch_mwacutils.sh \
          file_maxmin.py \
          find_pulsar_in_obs.py \
          fits_db.py \
          fix_offset.py \
          get_meta.py \
          get_voltage_obs.py \
          hdr.fits \
          import_mwacutils.sh \
          obs_query.py \
          plot_BPcal_128T.py \
          prepare.py \
          process_all.py \
          process_vcs.py \
          recombine.py \
          rename_corr_output.py \
          reorder_chans.py \
          rts2ao.py \
          run_rts.sh \
          test_mwacutils.sh \
          untar.sh \
          vcs_obs_w_files.txt

DATABASE_TARGETS = submit_to_database.py\
         database_vcs.py

UTILITY_TARGETS = zapchan.py\
	calc_ephem.py \
	check_disk_usage.sh \
	check_quota.sh \
	plot_tilebeams.py \
	mdir.py \
	mwa_metadb_utils.py

TARGETS = $(addprefix scripts/, $(SCRIPT_TARGETS)) $(addprefix database/, $(DATABASE_TARGETS)) $(addprefix utils/, $(UTILITY_TARGETS))

install: $(TARGETS)
	cp $^ $(INSTALL_DIR)

clean: $(addprefix $(INSTALL_DIR), $(SCRIPT_TARGETS)) $(addprefix $(INSTALL_DIR), $(DATABASE_TARGETS)) $(addprefix $(INSTALL_DIR), $(UTILITY_TARGETS))
	$(RM) $^
