#! /usr/bin/env python3
"""
Tests the sn_flux_est.py script
"""
from numpy.testing import assert_almost_equal
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
    test_cases.append((0, 0, 1223042487.0, 1.0))
    test_cases.append((1e20, 1e30, 0.0, -1.0000000000877696e-10))
       
    for beg, end, exp_enter, exp_exit in test_cases:
        out_enter, out_exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg, end)
        assert_almost_equal(exp_enter, out_enter, decimal=6)
        assert_almost_equal(exp_exit, out_exit, decimal=6)
 
def test_find_times()
    """
    Tests the find_times function
    """
    obsid = 1223042480
    pulsar = "J2234+2114"
    test_cases = []
    test_cases.append((None, None, 1223042487.0, 1223042844.3314152, 4913))
    test_cases.append((1223042487, 1223042587, 1223042487, 1223042587, 100.0))
    test_cases.append((-1e10, -1e4, , , ))
    test_cases.append((0, 0, , , ))
    test_cases.append((1e20, 1e30, , , ))

    for beg, end, exp_beg, exp_end, exp_int in test_cases:
        out_beg, out_end, out_int = snfe.find_times(obisd, pulsar, beg=beg, end=end)
        assert_almost_equal(exp_beg, out_beg, decimal=6)
        assert_almost_equal(exp_end, out_end, decimal=6)
        assert_almost_equal(exp_int, out_int, decimal=6)

        
def test_find_t_sys_gain():
    pass

if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
