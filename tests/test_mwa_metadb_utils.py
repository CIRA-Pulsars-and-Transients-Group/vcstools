#! /usr/bin/env python3
"""
Tests the mwa_metadb_utils.py script
"""
from vcstools.metadb_utils import mwa_alt_az_za, getmeta, get_obs_array_phase
from numpy.testing import assert_almost_equal


def test_mwa_alt_az_za():
    """Test the mwa_alt_az_za function"""
    # obsid, alt, az, za
    tests = [(1117101752, -71.10724927808731, 145.74310748819693, 161.1072492780873),
             (1252177744, -14.709910184536241, 264.22976419794514, 104.70991018453624),
             (1247832024, 68.90133642304133, 161.50105995238945, 21.09866357695867),
             (1188439520, 60.78396503767497, 161.03537536398974, 29.216034962325033)]

    for obsid, exp_alt, exp_az, exp_za in tests:
        alt, az, za = mwa_alt_az_za(obsid)
        assert_almost_equal(alt, exp_alt, decimal=3)
        assert_almost_equal(az,  exp_az,  decimal=3)
        assert_almost_equal(za,  exp_za,  decimal=3)

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
