#! /usr/bin/env python3
"""Tests the mwa_metadb_utils.py script."""
from vcstools.metadb_utils import mwa_alt_az_za, get_obs_array_phase, get_common_obs_metadata,\
                                  get_files, calc_ta_fwhm, get_channels, obs_max_min
from numpy.testing import assert_almost_equal

# Can't test singles_source_search or find_obsids_meta_pages because it's time consuming and always updating

def test_get_obs_array_phase():
    """Test FWHM calculation."""
    for obsid, expect_ans in [(1117101752, 'P1'),
                              (1252177744, 'P2C'),
                              (1247832024, 'P2E'),
                              (1188439520, 'OTH')]:
        ans = get_obs_array_phase(obsid)
        if ans != expect_ans:
            raise AssertionError()


def test_mwa_alt_az_za():
    """Test the mwa_alt_az_za function."""
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

def test_get_common_obs_metadata():
    """Test the get_common_obs_metadata function."""
    expect = [1117101752, 141.662872328864, 10.6540169083522, 4800,
              [[21, 19, 17, 15, 16, 14, 12, 10, 11, 9, 7, 5, 6, 4, 2, 0],
               [21, 19, 17, 15, 16, 14, 12, 10, 11, 9, 7, 5, 6, 4, 2, 0]],
              184.96, [133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
                       145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156]]
    ans = get_common_obs_metadata(1117101752)
    if ans != expect:
        raise AssertionError

def test_get_files():
    """Test the get_files function."""
    if get_files(1266823286) != ['1266823286_metafits_ppds.fits']:
        raise AssertionError

def test_calc_ta_fwhm():
    """Test the calc_ta_fwhm function."""
    tests = [(1133775752, 4.208367616629838e-08),
             (1164110416, 4.098704979920858e-08),
             (1194350120, 3.9949300287565104e-08)]

    for obsid, expect_fwhm in tests:
        ans = calc_ta_fwhm(obsid)
        assert_almost_equal(ans, expect_fwhm, decimal=10)

def test_get_channels():
    """Test the get_channels function."""
    if get_channels(1133775752) != [133, 134, 135, 136, 137, 138,
                                    139, 140, 141, 142, 143, 144,
                                    145, 146, 147, 148, 149, 150,
                                    151, 152, 153, 154, 155, 156]:
        raise AssertionError

def test_obs_max_min():
    """Test the obs_max_min function."""
    if obs_max_min(1133775752) != (1133775759, 1133780672):
        raise AssertionError

if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
