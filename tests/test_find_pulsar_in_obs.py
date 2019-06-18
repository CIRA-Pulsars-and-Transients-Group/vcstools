#! /usr/bin/env python3
"""
Tests the find_pulsar_in_obs.py script
"""
import find_pulsar_in_obs as fpio
from numpy.testing import assert_approx_equal, assert_almost_equal

def test_calcFWHM():
    """Test FWHM calculation"""
    for freq, expect_ans in [(150., 0.599584),
                             (180., 0.499654)]:
        ans = fpio.calcFWHM(freq)
        assert_almost_equal(ans, expect_ans, decimal=6)


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
