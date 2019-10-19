#! /usr/bin/env python3
"""
Tests the find_pulsar_in_obs.py script
"""
#hack to force the iers download from a url that works
from astropy.utils import iers
iers.IERS_A_URL = 'https://astroconda.org/aux/astropy_mirror/iers_a_1/finals2000A.all'
from astroplan import download_IERS_A
download_IERS_A()


import process_vcs

def test_gps_to_utc():
    """test the gps_to_utc function"""
    tests = [(1111111111, '2015-03-23T01:58:15'), (1222222222, '2018-09-29T02:10:04')]
    for obs, exp_utc in tests:
        utc = process_vcs.gps_to_utc(obs)
        if utc != exp_utc:
            raise AssertionError


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
