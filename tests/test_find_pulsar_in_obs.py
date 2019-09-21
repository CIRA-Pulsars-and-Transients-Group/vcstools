#! /usr/bin/env python3
"""
Tests the find_pulsar_in_obs.py script
"""
import find_pulsar_in_obs as fpio
from numpy.testing import assert_approx_equal, assert_almost_equal

def test_get_psrcat_ra_dec():
    """Test the psrqpy query and the max DM of get_psrcat_ra_dec"""
    ans = fpio.get_psrcat_ra_dec(pulsar_list=['J0002+6216','J0006+1834','J1056-6258'], max_dm=300)
    # Removes J0002+6216 because it has no DM and J1056-6258 becausethe DM is over 300
    if ans != [['J0006+1834', '00:06:04.8', '+18:34:59']]:
        raise AssertionError()

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
