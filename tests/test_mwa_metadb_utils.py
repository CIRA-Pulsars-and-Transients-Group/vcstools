#! /usr/bin/env python3
"""
Tests the mwa_metadb_utils.py script
"""
from mwa_metadb_utils import mwa_alt_az_za, getmeta, get_common_obs_metadata, get_obs_array_phase 
from numpy.testing import assert_approx_equal, assert_almost_equal

def test_get_obs_array_phase():
    """Test FWHM calculation"""
    for obsid, expect_ans in [(1117101752, 'P1'),
                              (1252177744, 'P2C'),
                              (1247832024, 'P2E'),
                              (1188439520, 'OTH')]:
        ans = get_obs_array_phase(obsid)
        if ans != expect_ans:
            raise AssertionError()


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
