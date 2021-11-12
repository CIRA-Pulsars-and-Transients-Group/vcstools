#!/usr/bin/env python3
import logging
import argparse

from vcstools.metadb_utils import get_obs_array_phase, obs_max_min, calc_ta_fwhm,\
                                  get_best_cal_obs, get_common_obs_metadata, files_available, getmeta
from vcstools.beam_calc import field_of_view

logger = logging.getLogger(__name__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Returns information on an input Obs ID""")
    parser.add_argument("obsid", type=int, help="Input Observation ID")
    parser.add_argument("-c", "--cal_best", action="store_true", help="If this option is used it will list "
                        "calibration observations within 2 days that have the same observing channels and "
                        "list them from closest in time to furthest.")
    args = parser.parse_args()

    beam_common_data, data_dict = get_common_obs_metadata(args.obsid, return_all=True)
    channels = data_dict["rfstreams"]["0"]["frequencies"]
    centre_freq = ( min(channels) + max(channels) ) / 2. * 1.28
    array_phase = get_obs_array_phase(args.obsid)
    fov = field_of_view(args.obsid, common_metadata=beam_common_data)

    print("-------------------------    Obs Info    -------------------------")
    print("Obs Name:           {}".format(data_dict["obsname"]))
    print("Creator:            {}".format(data_dict["rfstreams"]["0"]["creator"]))
    print("Array phase:        {}".format(array_phase))
    print("RA Pointing (deg):  {:6.2f}".format(data_dict["metadata"]["ra_pointing"]))
    print("DEC Pointing (deg): {:6.2f}".format(data_dict["metadata"]["dec_pointing"]))
    print("Centrefreq (MHz):   {}".format(centre_freq))
    print("Channels:           {}".format(data_dict["rfstreams"]["0"]["frequencies"]))
    if array_phase != 'OTH':
        print("~FWHM (arcminute):  {:4.2f}".format(calc_ta_fwhm(centre_freq,
                                                    array_phase=array_phase)*60.))
    print("FoV (square deg):   {:5.1f}".format(fov))


    # Perform file metadata calls
    files_meta_data = getmeta(servicetype='metadata', service='data_files', params={'obs_id':str(args.obsid)})
    start, stop = obs_max_min(args.obsid, files_meta_data=files_meta_data)
    # Check available files
    available_files, all_files = files_available(args.obsid, files_meta_data=files_meta_data)
    # Split into raw and combined files to give a clearer idea of files available
    available_comb = []
    available_raw  = []
    available_ics  = []
    all_comb = []
    all_raw  = []
    all_ics  = []
    for file_name in all_files:
        if 'combined' in file_name:
            all_comb.append(file_name)
            if file_name in available_files:
                available_comb.append(file_name)
        elif 'ics' in file_name:
            all_ics.append(file_name)
            if file_name in available_files:
                available_ics.append(file_name)
        else:
            all_raw.append(file_name)
            if file_name in available_files:
                available_raw.append(file_name)

    if start is None or stop is None:
        print("No VCS files found")
    else:
        print("Start time:         {}".format(start))
        print("Stop time:          {}".format(stop))
        print("Duration (s):       {}".format(stop-start))
    if len(all_raw) != 0:
        print("Raw  files avail:   {:5.1f}%  ({}/{})".format(len(available_raw)/len(all_raw)*100,
                                                                len(available_raw), len(all_raw)))
    if len(all_ics) != 0:
        print("ICS  files avail:   {:5.1f}%  ({}/{})".format(len(available_ics)/len(all_ics)*100,
                                                                len(available_ics), len(all_ics)))
    if len(all_ics) != 0:
        print("Comb files avail:   {:5.1f}%  ({}/{})".format(len(available_comb)/len(all_comb)*100,
                                                                len(available_comb), len(all_comb)))

    if args.cal_best:
        all_cals = get_best_cal_obs(args.obsid)
        print()
        print("{:14}|{:8}|{}".format("Calibration ID", "Hrs away", "Cal Target"))
        for cal in all_cals:
            print("{:14}|{:8.2f}|{}".format(cal[0], cal[1]/60., cal[2]))