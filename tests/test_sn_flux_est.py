#! /usr/bin/env python3
"""
Tests the sn_flux_est.py script
"""
from numpy.testing import assert_almost_equal
import mwa_metadb_utils
import sn_flux_est as snfe

def test_pulsar_beam_coverage():
    """
    Tests the pulsar_beam_coverage funtion
    """
    obsid = 1223042480
    pulsar = "J2234+2114"
    test_cases = []
    test_cases.append((None, None, 0.0, 0.072731816627943369))
    test_cases.append((1223042487, 1223042587, 0.0, 1.0))
    test_cases.append((-1e10, -1e4, 0.0, 0.072731816627943369))
    test_cases.append((0, 0, -1, -1))
    test_cases.append((1e20, 1e30, -1, -1))
       
    for beg, end, exp_enter, exp_exit in test_cases:
        out_enter, out_exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg, end)
        assert_almost_equal(exp_enter, out_enter, decimal=6)
        assert_almost_equal(exp_exit, out_exit, decimal=6)
 
def test_find_times():
    """
    Tests the find_times function
    """
    obsid = 1223042480
    pulsar = "J2234+2114"
    test_cases = []
    test_cases.append((None, None, 1223042487.0, 1223042844.3314152, 4913))
    test_cases.append((1223042487, 1223042587, 1223042487, 1223042587, 100.0))
    test_cases.append((-1e10, -1e4, 1223042487, 1223047400, 357.33141509308575))
    test_cases.append((0, 0, 1223042487, 1223047400, 0.0))
    test_cases.append((1e20, 1e30, 1e+20, 1e+30, 0.0))

    for beg, end, exp_beg, exp_end, exp_int in test_cases:
        out_beg, out_end, out_int = snfe.find_times(obsid, pulsar, beg=beg, end=end)
        assert_almost_equal(exp_beg, out_beg, decimal=6)
        assert_almost_equal(exp_end, out_end, decimal=6)
        assert_almost_equal(exp_int, out_int, decimal=6)

        
def test_find_t_sys_gain():
    """
    Tests the find_t_sys_gain function
    """
    obsid = 1223042480
    pulsar = "J2234+2114" 
    md = mwa_metadb_utils.get_common_obs_metadata(obsid) #this is to speed up the tests

    test_cases = []
    test_cases.append((None, None, None, None,\
                325.17339020027805, 6.5034678040055613, 0.060789018928146123, 0.047390137981582696))
    test_cases.append((1223042887, None, None, md,\
                 325.17339020027805, 6.5034678040055613, 0.049205077008620583, 0.038512861880721581))
    test_cases.append((None, "22:00:00", "20:00:00", md,\
                325.17339020027805, 6.5034678040055613, 0.014811201472754599, 0.011868437146105981))

    for beg, p_ra, p_dec, metadata, exp_t_sys, exp_t_sys_err, exp_gain, exp_gain_err in test_cases:
        t_sys, t_sys_err, gain, gain_err = snfe.find_t_sys_gain(\
                pulsar, obsid, beg=beg, p_ra=p_ra, p_dec=p_dec, obs_metadata=metadata)
        assert_almost_equal(exp_t_sys, t_sys, decimal=6)
        assert_almost_equal(exp_t_sys_err, t_sys_err, decimal=6)
        assert_almost_equal(exp_gain, gain, decimal=6)
        assert_almost_equal(exp_gain_err, gain_err, decimal=6)


def test_find_pulsar_w50():
    """
    Tests the find_pulsar_w50 function
    """
    test_cases=[]
    test_cases.append(("J2241-5236", 7.0000000000000007e-05, 3.5000000000000004e-06))
    test_cases.append(("J0206-4028", 0.0049000000000000007, 0.00024500000000000005)) 
    test_cases.append(("J2222-0137", 0.00056999999999999998, 5.0000000000000004e-06))

    for psr, exp_w50, exp_w50_err in test_cases:
        w50, w50_err = snfe.find_pulsar_w50(psr)
        assert_almost_equal(exp_w50, w50, decimal=6)
        assert_almost_equal(exp_w50_err, w50_err, decimal=6)


def test_est_pulsar_flux():
    """
    Tests the est_pulsar_flux function
    """
    
    test_cases = []
    test_cases.append(("J2241-5236", 1224252736, None, 0.19431618769079156, 0.02246695016871639))
    test_cases.append(("J2241-5236", None, 1.5e8, 0.20695827573194084, 0.026340991400397213))
    test_cases.append(("J2330-2005", 1226062160, None, 0.061795706713084243, 0.035885388925934944))
    test_cases.append(("J0152-1637", None, 1.5e8, 0.062465964922226662, 0.009151549839825322))
    
    for psr, obsid, f_mean, exp_flux, exp_flux_err in test_cases:
        flux, flux_err = snfe.est_pulsar_flux(psr, obsid=obsid, f_mean=f_mean)
        assert_almost_equal(exp_flux, flux, decimal=6)
        assert_almost_equal(exp_flux_err, flux_err, decimal=6)

def test_est_pulsar_sn():
    """
    I have tested this function various pulsar/observation combinations and compared to the estimated S/N from bestrof files. Here are the results:

    1225462936, -b 1225462943 -e 1225463543 psr J0152-1637: The function estimated a S/N of 99.39 +/- 23.55. I compared this to what was estimated from the profile analyse_pulse_prof() and attained a S/N of 75.26 +/- 1.70
    
    1226062160, -b 1226062167 -e 1226062767 psr J2330-2005: Function estimated a S/N of 62.21 +/- 38.63. Compared this to the profile estimation and attained S/N of 51.51 +/- 1.16

    1225713560, -b 1223042487 -e 1223043087 psr J2241-5236: Function estimated a S/N of . Compared this to the profile estimation and attained S/N of 17.36 +/- 2.41
    
    """
    
    test_cases=[]
    #Has 3 fluxes on database
    test_cases.append(("J2241-5236", 1225713560, None, None, None, None,\
                    232.29490006709051, 87.173008519118738))
    #Has 7 fluxes on database
    test_cases.append(("J0152-1637", 1225462936, 1225462943, 1225463543, None, None,\
                    99.3856180501139, 23.555069097506475))
    #Has 7 fluxes on database
    test_cases.append(("J2145-0750", 1221832280, 1221832287, 1221832887, None, None,\
                    120.19178597146058, 43.267175447498218))
    test_cases.append(("J2145-0750", 1221832280, 0, 1, None, None,\
                    0.0, 0.0))
    #Has 6 fluxes on database
    test_cases.append(("J2330-2005", 1226062160, None, None, "23:00:00", "-20:00:00",\
                    25.403853956665806, 16.326028900697739))
    #Has no fluxes on database
    test_cases.append(("J2234-0944", 1223042480, 1223042487, 1223043087, None, None,\
                    ))
    for psr, obsid, beg, end, ra, dec, exp_sn, exp_sn_err in test_cases:
        sn, sn_err = snfe.est_pulsar_sn(psr, obsid, beg=beg, end=end, p_ra=ra, p_dec=dec)
        assert_almost_equal(exp_sn, sn, decimal=6)
        assert_almost_equal(exp_sn_err, sn_err, decimal=6)

if __name__ == "__main__":
    """    
    Tests the relevant functions in sn_flux_est.py 
    Uses psrcat version 1.59. Values may change for different versions
    """

    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
