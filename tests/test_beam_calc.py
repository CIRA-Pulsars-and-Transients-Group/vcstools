#! /usr/bin/env python
"""
Tests the radiometer_equation.py script
"""
import os
from numpy.testing import assert_almost_equal
from vcstools.metadb_utils import get_common_obs_metadata
import psrqpy

from vcstools import data_load
from vcstools.beam_calc import source_beam_coverage, source_beam_coverage_and_times

import logging
logger = logging.getLogger(__name__)


def test_source_beam_coverage_and_times():
    """
    Tests the source_beam_coverage_and_times function
    """
    print("source_beam_coverage_and_times")
    obsid = 1223042480
    pulsar = "J2234+2114"
    test_cases = []
    #                  files_beg, files_end,
    test_cases.append((1223042487, 1223043086,
    #        expected: dect_beg, dect_end, dect_beg_norm, dect_end_norm, files_beg_norm, files_end_norm, obs_beg, obs_end, obs_dur
                       1223042487.0, 1223042844.4768934, 0.0, 0.07274662057244538, 0.0, 0.5957948223749796, 1223042487, 1223047400, 4914))

    for files_beg, files_end, exp_dect_beg, exp_dect_end, exp_dect_beg_norm, exp_dect_end_norm,\
        exp_files_beg_norm, exp_files_end_norm, exp_obs_beg, exp_obs_end, exp_obs_dur in test_cases:
        dect_beg, dect_end, dect_beg_norm, dect_end_norm, files_beg_norm, files_end_norm,\
        obs_beg, obs_end, obs_dur = source_beam_coverage_and_times(obsid, pulsar,
                                   p_ra="22:34:56.64", p_dec="+21:14:18.80",
                                   files_beg=files_beg, files_end=files_end)
        dect_dur = dect_end - dect_beg + 1
        assert_almost_equal(dect_beg, exp_dect_beg, decimal=4)
        assert_almost_equal(dect_end, exp_dect_end, decimal=4)
        assert_almost_equal(dect_beg_norm,  exp_dect_beg_norm,  decimal=4)
        assert_almost_equal(dect_end_norm,  exp_dect_end_norm,  decimal=4)
        assert_almost_equal(files_beg_norm, exp_files_beg_norm, decimal=4)
        assert_almost_equal(files_end_norm, exp_files_end_norm, decimal=4)
        assert_almost_equal(obs_end, exp_obs_end, decimal=4)
        assert_almost_equal(obs_beg, exp_obs_beg, decimal=4)
        assert_almost_equal(obs_dur, exp_obs_dur, decimal=4)