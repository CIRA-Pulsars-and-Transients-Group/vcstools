#! /usr/bin/env python3
"""
Author: Nicholas Swainston
Creation Date: /03/2016

Some of the orignal code was created by Bradley Meyers

This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.

Two most comon uses for this code is to search a single observation ID for pulsars like so:
> find_pulsar_in_obs.py -o 1099414416
which will output a text file called 1099414416_analytic_beam.txt
or to search all observation IDs for a single pulsar like so:
> find_pulsar_in_obs.py -p J0534+2200 --obs_for_source
which will output a text file called J0534+2200_analytic_beam.txt
"""

__author__ = 'Nicholas Swainston'
__date__ = '2016-03-21'

import os
import sys
import argparse
import numpy as np
import csv

from vcstools.pointing_utils import sex2deg, format_ra_dec
from vcstools.radiometer_equation import est_pulsar_sn, est_pulsar_flux, multi_psr_snfe
from vcstools.metadb_utils import get_obs_array_phase, singles_source_search,\
                                  find_obsids_meta_pages, obs_max_min
from vcstools.catalogue_utils import grab_source_alog
from vcstools.beam_calc import find_sources_in_obs
from vcstools.general_utils import setup_logger


import logging
logger = logging.getLogger(__name__)

class NoSourcesError(Exception):

    """Raise when no sources are found for any reason"""
    pass


def yes_no(answer):
    yes = set(['Y','yes','y', 'ye', ''])
    no = set(['N','no','n'])

    while True:
        choice = input(answer).lower()
        if choice in yes:
           return True
        elif choice in no:
           return False
        else:
           logger.warning("Please respond with 'yes' or 'no'\n")


def cal_on_database_check(obsid):
    from mwa_pulsar_client import client
    web_address = 'https://mwa-pawsey-volt01.pawsey.org.au'
    auth = ('mwapulsar','veovys9OUTY=')
    detection_list = client.detection_list(web_address, auth)

    cal_used = False
    cal_avail = False
    for d in detection_list:
        if int(d[u'observationid']) == int(obsid):
            if d[u'calibrator'] is not None:
                cal_used = True
                #TODO add a check if there is a cal file option

    #No cal
    check_result = 'N'
    if cal_avail:
        #Cal available
        check_result = 'A'
    elif cal_used:
        #Cal used
        check_result = 'U'

    return check_result


def write_output_source_files(output_data,
                              beam='analytic', min_power=0.3, cal_check=False,
                              SN_est=False, flux_est=False, plot_est=False,
                              min_time=0):
    """
    Writes an ouput file using the output of find_sources_in_obs when obs_for_source is true.
    """
    for source in output_data:
        out_name = "{0}_{1}_beam.txt".format(source, beam)
        with open(out_name,"w") as output_file:
            output_file.write('#All of the observation IDs that the {0} beam model '
                              'calculated a power of {1} or greater for the source: '
                              '{2}\n'.format(beam, min_power, source))
            output_file.write('#Column headers:\n')
            output_file.write('#Obs ID: Observation ID\n')
            output_file.write('#Dur:    The duration of the observation in seconds\n')
            output_file.write('#Enter:  The fraction of the observation when '
                                        'the source entered the beam\n')
            output_file.write('#Exit:   The fraction of the observation when '
                                        'the source exits the beam\n')
            output_file.write('#Power:  The maximum zenith normalised power of the source.\n')
            output_file.write("#OAP:    The observation's array phase where P1 is the "
                                        "phase 1 array, P2C is the phase compact array "
                                        "and P2E is the phase 2 extended array.\n")
            output_file.write("#Freq:   The centre frequency of the observation in MHz\n")
            output_file.write("#Band:   Bandwidth of the observation in MHz. If it is greater "
                                        "than 30.72 than it is a picket fence observation\n")

            if SN_est:
                output_file.write("#S/N Est: An estimate of the expected signal to noise using ANTF flux desnities\n")
                output_file.write("#S/N Err: The uncertainty of S/N Est\n")
            if flux_est:
                output_file.write("#Flux Est: An estimate of the expected flux density (mJy) using ANTF flux desnities\n")
                output_file.write("#Flux Err: The uncertainty of Flux Est (mJy)\n")
            if cal_check:
                output_file.write('#Cal ID: Observation ID of an available '+\
                                            'calibration solution\n')
            output_file.write('#Obs ID   |Dur |Enter|Exit |Power| OAP | Freq | Band ')
            if SN_est:
                output_file.write("|S/N Est|S/N Err")
            if flux_est:
                output_file.write("|Flux Est|Flux Err")

            if cal_check:
                output_file.write("|Cal ID\n")
            else:
                output_file.write('\n')
            for data in output_data[source]:
                obsid, duration, enter, leave, max_power, freq, band = data
                if duration > min_time:
                    if SN_est:
                        beg, end = obs_max_min(obsid)
                    oap = get_obs_array_phase(obsid)
                    output_file.write('{} {:4d} {:1.3f} {:1.3f} {:1.3f}  {:.3}  {:6.2f} {:6.2f}'.\
                            format(obsid, duration, enter, leave, max_power, oap, freq, band))
                    if SN_est:
                        pulsar_sn, pulsar_sn_err, _, _ = est_pulsar_sn(source, obsid, beg, end, plot_flux=plot_est)
                        if pulsar_sn is None:
                            output_file.write('   None    None')
                        else:
                            output_file.write(' {:9.2f} {:9.2f}'.format(pulsar_sn, pulsar_sn_err))
                    if flux_est:
                        pulsar_flux, pulsar_flux_err = est_pulsar_flux(source, obsid, plot_flux=plot_est)
                        if pulsar_flux is None:
                            output_file.write('   None    None')
                        else:
                            output_file.write(' {:8.2f} {:8.2f}'.format(pulsar_flux*1000, pulsar_flux_err*1000))

                    if cal_check:
                        #checks the MWA Pulsar Database to see if the obsid has been
                        #used or has been calibrated
                        logger.info("Checking the MWA Pulsar Databse for the obsid: {0}".format(obsid))
                        cal_check_result = cal_on_database_check(obsid)
                        output_file.write("   {0}\n".format(cal_check_result))
                    else:
                        output_file.write("\n")
    return


def write_output_obs_files(output_data, obsid_meta,
                           beam='analytic', min_power=0.3,
                           cal_check=False,
                           SN_est=False, flux_est=False, plot_est=False,
                           min_time=0):
    """
    Writes an ouput file using the output of find_sources_in_obs when obs_for_source is false.
    """

    for on, obsid in enumerate(output_data):
        if SN_est or flux_est:
            beg, end = obs_max_min(obsid)
            psr_list = [el[0] for el in output_data[obsid]]
            sn_dict = multi_psr_snfe(psr_list, obsid, beg, end,
                                     min_z_power=min_power, plot_flux=plot_est)

        oap = get_obs_array_phase(obsid)
        out_name = "{0}_{1}_beam.txt".format(obsid, beam)
        with open(out_name,"w") as output_file:
            output_file.write('#All of the sources that the {0} beam model calculated a power'
                              'of {1} or greater for observation ID: {2}\n'.format(beam,
                              min_power, obsid))
            output_file.write('#Observation data :RA(deg): {0} DEC(deg): {1} Duration(s): '
                              '{2} Array Phase: {3}\n'.format(obsid_meta[on][1],
                                  obsid_meta[on][2], obsid_meta[on][3], oap))
            if cal_check:
                #checks the MWA Pulsar Database to see if the obsid has been
                #used or has been calibrated
                logger.info("Checking the MWA Pulsar Databse for the obsid: {0}".format(obsid))
                cal_check_result = cal_on_database_check(obsid)
                output_file.write("#Calibrator Availability: {0}\n".format(cal_check_result))
            output_file.write('#Column headers:\n')
            output_file.write('#Source: Pulsar Jname\n')
            output_file.write('#Enter:  The fraction of the observation when '+\
                                        'the source entered the beam\n')
            output_file.write('#Exit:   The fraction of the observation when '+\
                                        'the source exits the beam\n')
            output_file.write('#Power:  The maximum zenith normalised power of the source.\n')
            if SN_est:
                output_file.write("#S/N Est: An estimate of the expected signal to noise using ANTF flux desnities\n")
                output_file.write("#S/N Err: The uncertainty of S/N Est\n")
            if flux_est:
                output_file.write("#Flux Est: An estimate of the expected flux density (mJy) using ANTF flux desnities\n")
                output_file.write("#Flux Err: The uncertainty of Flux Est (mJy)\n")

            output_file.write('#Source    |Enter|Exit |Power')
            if SN_est:
                output_file.write('| S/N Est | S/N Err')
            if flux_est:
                output_file.write("|Flux Est|Flux Err")
            output_file.write('\n')

            for data in output_data[obsid]:
                pulsar, enter_beam, exit_beam, max_power = data
                if (exit_beam - enter_beam) * obsid_meta[on][3] > min_time:
                    output_file.write('{:11} {:1.3f} {:1.3f} {:1.3f} '.format(pulsar,
                                    enter_beam, exit_beam, max_power))
                    if SN_est:
                        pulsar_sn, pulsar_sn_err, _, _ = sn_dict[pulsar]
                        if pulsar_sn is None:
                            output_file.write('   None    None')
                        else:
                            output_file.write('{:9.2f} {:9.2f}'.format(pulsar_sn, pulsar_sn_err))
                    if flux_est:
                        _, _, pulsar_flux, pulsar_flux_err = sn_dict[pulsar]
                        if pulsar_flux is None:
                            output_file.write('    None     None')
                        else:
                            output_file.write('{:8.2f} {:8.2f}'.format(pulsar_flux*1000, pulsar_flux_err*1000))
                    output_file.write('\n')

    return


if __name__ == "__main__":
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    beam_models = ['analytic', 'advanced', 'full_EE', 'hyperbeam']
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
    This code is used to list the sources within the beam of observations IDs or using --obs_for_source list all the observations for each source. The sources can be input serval ways: using a list of pulsar names (--pulsar), using a complete catalogue file of pulsars (--dl_PSRCAT) or RRATs (--RRAT and --dl_RRAT), using a compatable catalogue (--in_cat with the help of --names and --coordstype) or using a RA and DEC coordinate (--coords). The observation IDs can be input (--obsid) or gathered from a directory (--FITS_dir). The default is to search all observation IDs from http://mwa-metadata01.pawsey.org.au/metadata/ that have voltages and list every known pulsar from PSRCAT in each observation ID.
    """)
    parser.add_argument('--obs_for_source', action='store_true',
                        help='Instead of listing all the sources in each observation it will list all of the observations for each source. For increased efficiency it will only search OBSIDs within the primary beam.')
    parser.add_argument('-b', '--beam', type=str, default = 'analytic',
                        help='Decides the beam approximation that will be used. Options: "analytic" the analytic beam model (2012 model, fast and reasonably accurate), "advanced" the advanced beam model (2014 model, fast and slighty more accurate) or "full_EE" the full EE model (2016 model, slow but accurate). " Default: "analytic"')
    parser.add_argument('-m', '--min_power', type=float, default=0.3,
                        help='The minimum fraction of the zenith normalised power that a source needs to have to be recorded. Default 0.3')
    parser.add_argument('--sn_est', action='store_true',
                        help='Make a expected signal to noise calculation using the flux densities from the ANTF pulsar catalogue and include them in the output file. Default: False.')
    parser.add_argument('--flux_est', action='store_true',
                        help='Make a expected flux density calculation using the flux densities from the ANTF pulsar catalogue and include them in the output file. Default: False.')
    parser.add_argument('--plot_est', action='store_true',
                        help='If used, will output flux estimation plots while sn_est arg is true. Default: False.')
    parser.add_argument('--output', type=str, default = './',
                        help='Chooses a file for all the text files to be output to. The default is your current directory')
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                                    choices=loglevels.keys(), default="INFO")
    parser.add_argument("-V", "--version", action="store_true", help="Print version and quit")

    #source options
    sourargs = parser.add_argument_group('Source options', 'The different options to control which sources are used. Default is all known pulsars.')
    sourargs.add_argument('-p', '--pulsar', type=str, nargs='*', default = None,
                          help='Searches for all known pulsars. This is the default. To search for individual pulsars list their Jnames in the format " -p J0534+2200 J0630-2834"')
    sourargs.add_argument('--max_dm', type=float, default = 250.,
                          help='The maximum DM for pulsars. All pulsars with DMs higher than the maximum will not be included in output files. Default=250.0')
    sourargs.add_argument('--source_type', type=str, default='Pulsar',
                          help="An astronomical source type from ['Pulsar', 'FRB', 'rFRB', 'GC', 'RRATs', Fermi] to search for all sources in their respective web catalogue.")
    sourargs.add_argument('--in_cat', type=str,
                          help='Location of source catalogue, must be a csv where each line is in the format "source_name, hh:mm:ss.ss, +dd:mm:ss.ss".')
    sourargs.add_argument('-c', '--coords', type=str, nargs='*',
                          help='String containing the source\'s coordinates to be searched for in the format "RA_DEC" "RA_DEC". Must be enterered as either: "hh:mm:ss.ss_+dd:mm:ss.ss" or "deg_-deg". Please only use one format.')
    #finish above later and make it more robust to incclude input as sex or deg and perhaps other coordinte systmes

    #observation options
    obargs = parser.add_argument_group('Observation ID options', 'The different options to control which observation IDs are used. Default is all observation IDs with voltages.')
    obargs.add_argument('--FITS_dir', type=str,
                        help='Instead of searching all OBS IDs, only searchs for the obsids in the given directory. Does not check if the .fits files are within the directory. Default = /group/mwavcs/vcs')
    obargs.add_argument('-o','--obsid', type=int, nargs='*',
                        help='Input several OBS IDs in the format " -o 1099414416 1095506112". If this option is not input all OBS IDs that have voltages will be used')
    parser.add_argument('--min_time', type=float, default=0,
                        help='The minimum observation duration to include in output files. Default 0')
    parser.add_argument('--freq_chan', type=int,
                        help='Only use observations that include this frequency channel.')
    parser.add_argument('--contig', action='store_true',
                        help='Only use observations that have contiguous frequency channels.')
    obargs.add_argument('--all_volt', action='store_true',
                        help='Includes observation IDs even if there are no raw voltages in the archive. Some incoherent observation ID files may be archived even though there are raw voltage files. The default is to only include files with raw voltage files.')
    obargs.add_argument('--cal_check', action='store_true',
                        help='Check the MWA Pulsar Database to check if the obsid has every succesfully detected a pulsar and if it has a calibration solution.')
    args=parser.parse_args()



    if args.version:
        try:
            import version
            print(version.__version__)
            sys.exit(0)
        except ImportError as ie:
            print("Couldn't import version.py - have you installed vcstools?")
            print("ImportError: {0}".format(ie))
            sys.exit(0)

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    #Parse options
    if args.in_cat and args.coords:
        logger.error("Can't use --in_cat and --coords. Please input your cooridantes "
                     "using one method. Exiting.")
        quit()
    if args.obsid and args.FITS_dir:
        logger.error("Can't use --obsid and --FITS_dir at the same time. Exiting.")
        quit()
    if args.beam not in beam_models:
        logger.error("Unknown beam model. Please use one of {0}. Exiting.".format(beam_models))
        quit()

    logger.info("Gathering sources")
    degrees_check = False
    if args.in_cat:
        names_ra_dec = []
        with open(args.in_cat,"r") as input_catalogue:
            reader = csv.reader(input_catalogue)
            for row in reader:
                names_ra_dec.append(row)
    elif args.coords:
        names_ra_dec = []
        for c in args.coords:
            names_ra_dec.append([c, c.split('_')[0], c.split('_')[1]])
            if ":" not in c:
                degrees_check = True
    else:
        names_ra_dec = grab_source_alog(source_type=args.source_type,
                                        pulsar_list=args.pulsar,
                                        max_dm=args.max_dm)
        if len(names_ra_dec) == 0:
            raise NoSourcesError(f"""No sources found in catalogue with:
                                    source type:    {args.source_type}
                                    pulsars:        {args.pulsar}
                                    max dm:         {args.max_dm}""")

    #format ra and dec
    if not degrees_check:
        names_ra_dec = format_ra_dec(names_ra_dec, ra_col=1, dec_col=2)
    names_ra_dec = np.array(names_ra_dec)

    #Check if the user wants to use --obs for source
    if (len(names_ra_dec) ==1) and (not args.obs_for_source):
        args.obs_for_source = yes_no('You are only searching for one pulsar so it is '+\
                                     'recommened that you use --obs_for_source. Would '+\
                                     'you like to use --obs_for_source. (Y/n)')
        if args.obs_for_source:
            logger.info("Using option --obs_for_source")
        else:
            logger.info("Not using option --obs_for_source")

    #get obs IDs
    logger.info("Gathering observation IDs")
    # A dictionary of constraints used to search for suitable observations as explained here:
    # https://wiki.mwatelescope.org/display/MP/Web+Services#WebServices-Findobservations
    params = {'mode':'VOLTAGE_START'}
    if args.freq_chan:
        # Update params to only use obs that contain this frequency channel
        params.update({'anychan':args.freq_chan})
        logger.debug("params: {}".format(params))
    if args.contig:
        # Update params to only use obs that have contiguous channels
        params.update({'contigfreq':1})
        logger.debug("params: {}".format(params))

    if args.obsid:
        obsid_list = args.obsid
    elif args.FITS_dir:
        obsid_list = os.walk(args.FITS_dir).next()[1]
    elif len(names_ra_dec) == 1:
        #if there is a single pulsar simply use nearby obsids
        if degrees_check:
            ob_ra = names_ra_dec[0][1]
            ob_dec = names_ra_dec[0][2]
        else:
            ob_ra, ob_dec = sex2deg(names_ra_dec[0][1], names_ra_dec[0][2])
        obsid_list = singles_source_search(ob_ra, params=params)
    else:
        #use all obsids
        obsid_list = find_obsids_meta_pages(params)


    if args.beam == 'full_EE':
        dt = 300
    else:
        dt = 100

    logger.debug("names_ra_dec:{}".format(names_ra_dec))
    logger.debug("obsid:{}".format(obsid_list))
    logger.info("Getting observation metadata and calculating the tile beam")
    output_data, obsid_meta = find_sources_in_obs(obsid_list, names_ra_dec,
                                obs_for_source=args.obs_for_source, dt_input=dt,
                                beam=args.beam, min_z_power=args.min_power,
                                cal_check=args.cal_check, all_volt=args.all_volt,
                                degrees_check=degrees_check)

    logger.info("Writing data to files")
    if args.obs_for_source:
        write_output_source_files(output_data,
                                  beam=args.beam, min_power=args.min_power,
                                  cal_check=args.cal_check,
                                  SN_est=args.sn_est, flux_est=args.flux_est, plot_est=args.plot_est,
                                  min_time=args.min_time)
    else:
        write_output_obs_files(output_data, obsid_meta,
                               beam=args.beam, min_power=args.min_power,
                               cal_check=args.cal_check,
                               SN_est=args.sn_est, flux_est=args.flux_est, plot_est=args.plot_est,
                               min_time=args.min_time)
