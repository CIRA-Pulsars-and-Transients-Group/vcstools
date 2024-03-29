#!/usr/bin/env python3
import logging
import argparse
import sys

from vcstools.radiometer_equation import est_pulsar_sn
from vcstools.pulsar_spectra import ATNF_spectral_data_plot
from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)


def snfe_main(kwargs):
    if kwargs["mode"] == "SNFE":
        raj = None
        decj = None
        if kwargs["pointing"]:
            raj = kwargs["pointing"].split("_")[0]
            decj = kwargs["pointing"].split("_")[1]
        SN, SN_e, flux, flux_e = est_pulsar_sn(kwargs["pulsar"], kwargs["obsid"],
                                plot_flux=kwargs["plot_est"], p_ra=raj, p_dec=decj)
        logger.info(f"{kwargs['pulsar']} Flux: {flux:.2f} +/- {flux_e:.2f} Jy")
        logger.info(f"{kwargs['pulsar']} S/N: {SN:.2f} +/- {SN_e:.2f}")
    elif kwargs["mode"] == "ATNF":
        ATNF_spectral_data_plot(kwargs["pulsar"])
    else:
        logger.error("Valid mode not selected. Please refer to documentation for options")
        sys.exit(1)


if __name__ == "__main__":

    loglevels = dict(
        DEBUG=logging.DEBUG,
        INFO=logging.INFO,
        WARNING=logging.WARNING,
        ERROR=logging.ERROR
    )

    parser = argparse.ArgumentParser(description="""A utility file for estimating the S/N of a pulsar in an obsid.
    Example command: sn_flux_est.py -o 1221399680 -p J2129-5721 --pointing 21:29:00_-57:21:00 --plot_est""")

    parser.add_argument("-o", "--obsid", type=int, required=True, help="The Observation ID (e.g. 1221399680)")
    parser.add_argument("-p", "--pulsar", type=str, required=True, help="The pulsar's name (e.g. J2241-5236).\
                        Takes a single argument for flux estimation or multiple for ATNF plotting")
    parser.add_argument("-b", "--beg", type=int, help="The beginning time of observation.\
                        If None, will use beginning given by a metadata call. Default: None")
    parser.add_argument("-e", "--end", type=int, help="The end time of observation.\
                        If None, will use the end given by a metadata call. Default: None")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbostity level. Default: INFO")
    parser.add_argument("--pointing", type=str, default=None, help="The pointing of the target in the format '12:34:56_98:76:54'.\
                        If None, will obtain from a call to ATNF. Default: None")
    parser.add_argument("--min_z_power", type=float, default=0.3, help="The minimum zenith normalised power used to determine if the pulsar\
                        is in the beam or not")
    parser.add_argument("--plot_est", action="store_true", help="Use this tag to create a plot of flux estimation.")
    parser.add_argument("--mode", type=str, default="SNFE",
                        help="MODES: 'SNFE' = Estimate S/N and flux for a single pulsar in an obsid.\
                             'ATNF' = Plot the spectral energy distribution for any number of pulsars using data from ATNF.")
    args = parser.parse_args()

    # set up the logger for stand-alone execution
    logger = setup_logger(logger, log_level=loglevels[args.loglvl])

    kwargs=vars(args)
    snfe_main(kwargs)
