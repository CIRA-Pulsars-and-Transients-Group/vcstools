#!/usr/bin/env python3
"""
Script to work out flagging of channels in PRESTO and PSRCHIVE formats.

Written by: Bradley Meyers
Originally created: 6 Nov 2016
Updated: 10 Nov 2017 (BWM)
"""

import sys
import argparse
import logging
import numpy as np
from itertools import groupby, count

from vcstools.general_utils import setup_logger

logger = logging.getLogger(__name__)


class ZapchanError(Exception):
    pass


def write_paz_file(ctf, fname):
    # small function to handle the file-writing side of things
    # checks if we can actually write, otherwise raises an error
    try:
        with open(fname, "w") as fileobj:
            for chan in ctf:
                if chan is not None:
                    fileobj.write("{0}\n".format(chan))
    except IOError:
        errmsg = "Couldn't open file to write to: {0}".format(fname)
        logger.erro(errmsg)
        raise ZapchanError(errmsg)


def zap(r, p, chans, zapedges=False, nzap=20, zapmid=False, zaps=None):
    """
    Main function to take in the user's requested masking information and create a
    usable string and/or file containing the channels that need to be flagged.
    """

    if chans[1] is not None:
        coarse = chans[0]
        fine = chans[1]
        total = coarse * fine
    else:
        fine = 128
        total = chans[0]
        coarse = int(total / fine)

    if zaps is not None:
        zaps = np.array(zaps, dtype=int)

    logger.debug("number of coarse channels: {0}".format(coarse))
    logger.debug("number of fine channels: {0}".format(fine))
    logger.debug("total number of channels: {0}".format(total))
    logger.debug("user-specified channels to zap: {0}".format(zaps))
    logger.debug("number of edge channels to zap on each side: {0}".format(nzap))
    logger.debug("zapping middle channel? {0}".format(zapmid))

    # Let's create a list of the channels, grouped into adjacent subbands, that need to be flagged
    # NOTE: This is very similar to what we do when figuring out which channels to automatically flag when producing
    # the RTS configuration files
    allchans = np.arange(total)  # list all fine channels
    splitchans = np.array_split(allchans, coarse)  # split it into 24 coarse channel chunks

    # now figure out which channels we want to flag and discard the rest
    chans_to_flag = []
    nzapped = 0
    for coarse_channel in splitchans:
        center_chan = coarse_channel[fine // 2]
        left_edge_chans = coarse_channel[:nzap]
        right_edge_chans = coarse_channel[-nzap:]
        if zapmid and zapedges:
            # zap middle and edge channels
            nzapped = 2 * nzap + 1  # actual number of channels zapped
            d = [left_edge_chans, center_chan, right_edge_chans]
            chans_to_flag.append(d)
        elif zapedges and (zapmid is False):
            # just zap edge channels
            nzapped = 2 * nzap
            d = [left_edge_chans, np.array([None]), right_edge_chans]
            chans_to_flag.append(d)
        elif zapmid and (zapedges is False):
            # just zap center channels
            nzapped = 1
            d = [np.array([None]), center_chan, np.array([None])]
            chans_to_flag.append(d)

        # the only other option is you asked for no edge of center channels to flagged,
        # so chans_to_flag should be empty

    # For the paz file: need to just flatten everything out into one list
    cchan_flags = []
    if chans_to_flag:
        for c in range(coarse):
            cchan_flags.append(np.concatenate([x.flatten() for x in chans_to_flag[c]]))

    # collapse it all into 1D
    p_ctf = np.array(cchan_flags).flatten()
    # now check to see if we need to include any of the user-define channels to remove
    nuniq = 0
    if zaps is not None:
        for z in zaps:
            if z in p_ctf:
                continue
            else:
                nuniq += 1
                p_ctf = np.append(p_ctf, z)

    # regardless of if we added something, remove the None values and convert to a sorted list
    # NOTE: filter(None, p_ctf) doesn't work here because it also removes 0's
    p_ctf = sorted([x for x in p_ctf if x is not None])

    # For the rfifind string: trickier formatting, so we need to account for the user-defined channel zaps first
    # We should now start from the p_ctf list, because that has ALL required channels, and we can re-split based
    # on whether the channels are consecutive or not
    groups = [list(g) for k, g in groupby(p_ctf, lambda i, j=count(): next(j)-i)]
    r_ctf = ""
    for g in groups:
        if len(g) > 1:
            r_ctf += "{0}:{1},".format(str(min(g)), str(max(g)))
        else:
            r_ctf += "{0},".format(str(g[0]))

    # and remove the last character as it'll be an extraneous comma
    r_ctf = r_ctf[:-1]

    logger.debug("total number of edge and/or center channels zapped per coarse channel: {0}".format(nzapped))
    logger.debug("number of user-defined channels zapped not within edges: {0}".format(nuniq))
    logger.debug("total number of zapped channels: {0}".format(len(p_ctf)))
    logger.debug("fraction of band masked: {0:.1f}%".format(100 * len(p_ctf) / float(total)))

    # We now have a list of channels that need to be flagged in both formats
    if p and not r:
        logger.debug("p and not r")
        write_paz_file(p_ctf, p)
    elif r and not p:
        logger.debug("r and not p")
        print(r_ctf)
    elif r and p:
        logger.debug("r and p")
        write_paz_file(p_ctf, p)
        print(r_ctf)
    else:
        logger.debug("not r and not p")
        # didn't specify r or p, so print nothing
        print("")


if __name__ == '__main__':
    # SETUP OPTIONS #

    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)

    parser = argparse.ArgumentParser(description="""Output the correctly formated MWA band edge channels and
                                    user-defined channels to remove for aliasing/RFI excision.
                                    If no options are given, returns an empty string.""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # there should be two options for the channels, either:
    # give me the total number of fine channels (and I'll assume there's 128 fine channels per coarse channel)
    # OR
    # give me the number of coarse channels and the number of fine channels and I'll use that information to then
    # calculate the intervals appropriate.
    parser.add_argument("-N", "--nchans", type=int, metavar='Nchan',
                        help="TOTAL number of frequency channels in data set. Assumes 128 fine channels  per "
                             "coarse channel.", default=None)

    parser.add_argument("-C", "--n_coarse_fine", type=int, nargs=2, metavar="Nchan",
                        help="Number of coarse channels followed by the number of fine channels per coarse channel. "
                             "This options overrides --nchans.", default=None)

    # zap the edge fine channels?
    parser.add_argument("-Z", "--zapedges", action='store_true',
                        help="Zap the fine channel edges of the coarse channels.")

    # if so, how many?
    parser.add_argument("-n", "--nzap", type=int, metavar="Nedge",
                        help="Number of edge channels to remove from each side of a coarse channel. "
                             "If given, but --zapedges is not then this argument is ignored.", default=20)

    # zap the middle channels?
    parser.add_argument("-m", "--middle", action='store_true',
                        help="Flag the center fine channels for each coarse channel?")

    # user-defined channels to zap (will be prioritised over edge channel zapping)
    parser.add_argument("-z", type=str, metavar="chan", nargs="+", help="Individual channels to zap")

    # what version of output do you want? PRESTO and/or PSRCHIVE
    parser.add_argument("-r", "--rfifind", action='store_true',
                        help="Output to screen the edge channels in a format readable by the PRESTO rfifind routine",
                        default=False)

    parser.add_argument("-p", "--paz", metavar="filename", type=str,
                        help="Output channels for PSRCHIVE paz routine into the given file, using paz's '-k filename' "
                             "option")

    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: WARNING",
                        choices=loglevels.keys(), default="WARNING")

    parser.add_argument("-V", "--version", action='store_true', help="Print version and quit")
    args = parser.parse_args()

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

    # check to make sure that the channels variable is formatted correctly,
    # no matter the input (should be a tuple/iterable)
    if args.nchans and not args.n_coarse_fine:
        channels = (args.nchans, None)
        logger.debug("assuming 128 fine channels per coarse channel")
    else:
        logger.debug("assuming values from -C argument and ignoring -N value")
        channels = args.n_coarse_fine

    # TODO: should look at moving some of these exit failures to within the function and raise a ZapchanError
    # check to make sure we can actually figure out how many channels we have and
    # how they are split up
    if channels is None:
        logger.warning("You didn't provide any way to compute the number of channels. "
                           "Please use either the -N or -C options")
        sys.exit(0)

    # check to make sure that if nzap is given but user has not elected
    # to zap edges, we zap 0 edge channels
    if args.nzap and not args.zapedges:
        args.nzap = 0

    if args.zapedges is False and args.middle is False and args.z is None:
        logger.warning("No flagging requested (i.e. none of -Z, -z or -m were given). Aborting.")
        sys.exit(0)

    # DO THE THING! #
    zap(args.rfifind, args.paz, channels, args.zapedges, args.nzap, args.middle, args.z)
