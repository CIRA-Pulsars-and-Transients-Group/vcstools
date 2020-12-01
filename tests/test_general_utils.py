#! /usr/bin/env python3
"""Tests the find_pulsar_in_obs.py script."""
from vcstools.general_utils import gps_to_utc, is_number

def test_is_number():
    """test the is_number function."""
    tests = [(1111111111, True), ('abc', False)]
    for input_test, expect in tests:
        ans = is_number(input_test)
        if ans != expect:
            raise AssertionError

def test_gps_to_utc():
    """test the gps_to_utc function."""
    tests = [(1111111111, '2015-03-23T01:58:15'), (1222222222, '2018-09-29T02:10:04')]
    for obs, exp_utc in tests:
        utc = gps_to_utc(obs)
        if utc != exp_utc:
            raise AssertionError


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
