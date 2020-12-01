#! /usr/bin/env python3
"""
Tests the find_pulsar_in_obs.py script
"""
import find_pulsar_in_obs as fpio
from numpy.testing import assert_approx_equal, assert_almost_equal

def test_get_psrcat_ra_dec():
    """Test the psrqpy query and the max DM of get_psrcat_ra_dec"""
    ans = fpio.get_psrcat_ra_dec(pulsar_list=['J0007+7303','J0006+1834','J1056-6258'], max_dm=300)
    # Removes J0007+7303 because it has no DM and J1056-6258 becausethe DM is over 300
    if ans != [['J0006+1834', '00:06:04.8', '+18:34:59']]:
        raise AssertionError()

def test_get_source_alog():
    """Test get_source_alog"""
    
    #source_type, pulsar_list, include_dm, answer
    tests = [['Pulsar','J2313+4253'   , False, [['J2313+4253', '23:13:08.6209', '+42:53:13.043']]],
             ['Pulsar','J2313+4253'   , True,  [['J2313+4253', '23:13:08.6209', '+42:53:13.043', 17.27693]]],
             ['FRB'   ,'FRB171019'    , False, [['FRB171019', '22:17.5', '-08:40']]],
             ['FRB'   ,'FRB171019'    , True,  [['FRB171019', '22:17.5', '-08:40', '460.8']]],
             ['rFRB'  ,'FRB171019'    , False, [['FRB171019', '22:17:30', '-08:40']]],
             ['rFRB'  ,'FRB171019'    , True,  [['FRB171019', '22:17:30', '-08:40', '460.8']]],
             ['RRATs' ,'J1913+1330'   , False, [['J1913+1330', '19:13:17', '13:30:32.8']]],
             ['RRATs' ,'J1913+1330'   , True,  [['J1913+1330', '19:13:17', '13:30:32.8', '175.64']]]]
             # Removing because can't test them on travis without making the candidates public
             #['Fermi' ,'J2219.7-6837' , False, [['J2219.7-6837', '22:19:47.256', '-68:37:02.28']]],
             #['Fermi' ,'J2219.7-6837' , True,  [['J2219.7-6837', '22:19:47.256', '-68:37:02.28', 2.43]]]]

    for test in tests:
        stype, name, dm, expected_ans = test
        ans = fpio.grab_source_alog(source_type=stype, pulsar_list=[name], include_dm=dm)
        if ans != expected_ans:
            raise AssertionError()




if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
