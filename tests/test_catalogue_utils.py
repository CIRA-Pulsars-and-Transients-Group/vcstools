#! /usr/bin/env python3
"""
Tests the find_pulsar_in_obs.py script
"""
import pytest
from vcstools.catalogue_utils import get_psrcat_ra_dec, grab_source_alog, get_rFRB_info

def test_get_psrcat_ra_dec():
    """Test the psrqpy query and the max DM of get_psrcat_ra_dec."""
    ans = get_psrcat_ra_dec(pulsar_list=['J0007+7303','J0006+1834','J1056-6258'], max_dm=300)
    # Removes J0007+7303 because it has no DM and J1056-6258 becausethe DM is over 300
    if ans != [['J0006+1834', '0:06:04.8', '18:34:59']]:
        raise AssertionError()

def test_get_source_alog():
    """Test get_source_alog."""
    #source_type, pulsar_list, include_dm, answer
    tests = [
             ['Pulsar','J2313+4253'   , False, [['J2313+4253', '23:13:08.6209', '42:53:13.043']]],
             ['Pulsar','J2313+4253'   , True,  [['J2313+4253', '23:13:08.6209', '42:53:13.043', 17.27693]]],
             # Below are removed until we fix issue #214
             #['FRB'   ,'FRB_20171019A', False, [['FRB_20171019A', '22:17:30.000', '-08:40:00.00']]],
             #['FRB'   ,'FRB_20171019A', True,  [['FRB_20171019A', '22:17:30.000', '-08:40:00.00', '460.8']]],
             #['FRB'   ,'FRB_20210630A', False, [['FRB_20210630A', '17:23:07.409', '+07:51:41.85']]],
             #['FRB'   ,'FRB_20210630A', True, [['FRB_20210630A', '17:23:07.409', '+07:51:41.85', '943.7']]],
             #['rFRB'  ,'FRB171019'    , False, [['FRB171019', '22:17:30', '-08:40']]],
             #['rFRB'  ,'FRB171019'    , True,  [['FRB171019', '22:17:30', '-08:40', '460.8']]],
             ['RRATs' ,'J1913+1330'   , False, [['J1913+1330', '19:13:17', '13:30:32.8']]],
             ['RRATs' ,'J1913+1330'   , True,  [['J1913+1330', '19:13:17', '13:30:32.8', '175.64']]],
             # Removing because can't test them on travis without making the candidates public
             #['Fermi' ,'J2219.7-6837' , False, [['J2219.7-6837', '22:19:47.256', '-68:37:02.28']]],
             #['Fermi' ,'J2219.7-6837' , True,  [['J2219.7-6837', '22:19:47.256', '-68:37:02.28', 2.43]]]
            ]

    for test in tests:
        stype, name, dm, expected_ans = test
        ans = grab_source_alog(source_type=stype, pulsar_list=[name], include_dm=dm)
        if ans != expected_ans:
            raise AssertionError()

@pytest.mark.skip(reason="need to reevaluate access to FRB data via TNS")
def test_get_rFRB_info():
    """Test get_rFRB_info."""
    expected_all = [['FRB171019', '22:17:30', '-08:40', '460.8', '1.1\n'],
                    ['FRB121102', '05:32:09', '33:05:13', '557', '2\n'],
                    ['FRB180814.J1554+74', '15:54', '+74:01', '238.32', '0.01\n'],
                    ['FRB180916.J0158+65', '01:58', '65:44', '349.2', '0.4\n'],
                    ['FRB181030.J1054+73', '10:54', '73:44', '103.5', '0.7\n'],
                    ['FRB181128.J0456+63', '04:56', '63:23', '450.2', '0.3\n'],
                    ['FRB181119.J12+65', '12:42', '65:08', '364.2', '1\n'],
                    ['FRB190116.J1249+27', '12:49', '27:09', '444', '0.6\n'],
                    ['FRB181017.J1705+68', '17:05', '68:17', '1281.9', '0.4\n'],
                    ['FRB190209.J0937+77', '09:37', '77:40', '424.6', '0.6\n'],
                    ['FRB190222.J2052+69', '20:52', '69:50', '460.6', '0.1\n']]
    ans_all = get_rFRB_info()
    if ans_all != expected_all:
        raise AssertionError()

    expected_one = [['FRB171019', '22:17:30', '-08:40', '460.8', '1.1\n']]
    ans_one = get_rFRB_info(name='FRB171019')
    if ans_one != expected_one:
        raise AssertionError()


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
