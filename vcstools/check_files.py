import subprocess
import os
import numpy as np
import traceback
import logging
from vcstools.metadb_utils import getmeta, get_files
from vcstools.config import load_config_file

logger = logging.getLogger(__name__)

def check_download(obsID, directory=None, startsec=None, n_secs=None, data_type='raw'):
    '''
    Checks that the number of files in directory (default is /astro/mwavcs/vcs/[obsID]/raw/) is the same
    as that found on the archive and also checks that all files have the same size (253440000 for raw, 7864340480 for recombined tarballs by default).
    '''
    comp_config = load_config_file()
    if not data_type in ['raw', 'tar_ics', 'ics']:
        logger.error("Wrong data type given to download check.")
        return True
    if not directory:
        directory = os.path.join(comp_config['base_data_dir'], str(obsID), "raw") if data_type == 'raw' else os.path.join(comp_config['base_data_dir'], str(obsID), "combined")
    base = "\n Checking file size and number of files for obsID {0} in {1} for ".format(obsID, directory)
    n_secs = n_secs if n_secs else 1
    logger.info(base + "gps times {0} to {1}".format(startsec, startsec+n_secs-1) if startsec else base + "the whole time range.")

    # put files in
    try:
        files, suffix, required_size = get_files_and_sizes(obsID, data_type, mintime=startsec, maxtime=startsec + n_secs)
    except:
        return True

    if not startsec:
        n_files_expected = len(files)
        command = "ls -l %s/*%s | ((tee /dev/fd/5 | wc -l >/dev/fd/4) 5>&1 | " %(directory, suffix) + \
            "awk '($5!=%s){print \"file \" $9 \" has size \" $5 \" (expected %s)\"}' >> %s/%s_all.txt) 4>&1;" %(required_size, required_size,directory, obsID) + \
            "cat %s/%s_all.txt; rm -rf %s/%s_all.txt" %(directory, obsID, directory, obsID)
        output = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).stdout
    else:
        n_files_expected = 0
        #remove stray metafits from list that causes int errors
        files = [ x for x in files if "metafits" not in x ]
        times = [int(time[11:21]) for time in files]
        for sec in range(startsec,startsec+n_secs):
            n_files_expected += times.count(sec)
        output = subprocess.Popen(["count=0;for sec in `seq -w %s %s `;do let count=${count}+`ls -l %s/*${sec}*%s | " %(startsec, startsec+n_secs-1, directory, suffix) + \
                                       "((tee /dev/fd/5 | wc -l >/dev/fd/4) 5>&1 | awk '($5!=%s) " %(required_size) + \
                                       "{print \"file \" $9 \" has size \" $5 \" (expected %s)\"}' >> %s/errors_%s.txt) 4>&1`;done;" %(required_size,directory,startsec) +\
                                       "echo ${count}; cat %s/errors_%s.txt;rm -rf %s/errors_%s.txt" %(directory,startsec,directory,startsec)],
                                  stdout=subprocess.PIPE, shell=True).stdout
    output = output.readlines()
    files_in_dir = int(output[0].strip())

    error = False

    # in case we're checking for downloaded tarballs also need to check ics-files.
    if data_type == 'tar_ics':
        logger.info("Now checking ICS files")
        error, n_ics = check_recombine_ics(directory=directory,
                                           startsec=startsec,
                                           n_secs=n_secs,#n_files_expected,
                                           obsID=obsID)
        n_files_expected *= 2
        files_in_dir += n_ics


    if not files_in_dir == n_files_expected:
        logger.error("We have {0} files but expected {1}".format(files_in_dir, n_files_expected))
        error = True
    for line in output[1:]:
        if b'file' in line:
            logger.error(line)
            error = True
    if not error:
        logger.info("We have all {0} {1} files as expected.".format(files_in_dir, data_type))
    return error

def check_recombine(obsID, directory=None, required_size=327680000, \
                        required_size_ics=30720000, startsec=None, n_secs=None):
    '''
    Checks that the number of files in directory (/astro/mwavcs/vcs/[obsID]/combined/) is ....
    as that found on the archive and also checks that all files have the same size (327680000 by default).
    '''
    comp_config = load_config_file()
    if not directory:
        directory = os.path.join(comp_config['base_data_dir'], str(obsID), "combined")
    base = "\n Checking file size and number of files for obsID {0} in {1} for ".format(obsID, directory)
    n_secs = n_secs if n_secs else 1
    logger.info(base + "gps times {0} to {1}".format(startsec, startsec+n_secs-1) if startsec else base + "the whole time range.")
    required_size = required_size
    # we need to get the number of unique seconds from the file names
    files = np.array(get_files(obsID))
    mask = np.array(['.dat' in file for file in files])
    if not startsec:
        times = [time[11:21] for time in files[mask]]
        n_secs = len(set(times))
        command = "ls -l %s/*ch*.dat | ((tee /dev/fd/5 | wc -l >/dev/fd/4) 5>&1 | " %(directory) + \
            "awk '($5!=%s){print $9}' | tee >> %s/%s_all.txt | xargs rm -rf) 4>&1;" %(required_size, directory, obsID) + \
            "cat %s/%s_all.txt; rm -rf %s/%s_all.txt" %(directory, obsID, directory, obsID)
        output = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).stdout
    else:
        output = subprocess.Popen(["count=0;for sec in `seq -w %s %s `;do let count=${count}+`ls -l %s/*${sec}*ch*.dat | " %(startsec, startsec+n_secs-1, directory) + \
                                       "((tee /dev/fd/5 | wc -l >/dev/fd/4) 5>&1 | awk '($5!=%s) " %(required_size) + \
                                       "{print $9}' | tee >> %s/errors_%s.txt | xargs rm -rf) 4>&1`;done;" %(directory,startsec) +\
                                       "echo ${count}; cat %s/errors_%s.txt;rm -rf %s/errors_%s.txt" %(directory,startsec,directory,startsec)],
                                  stdout=subprocess.PIPE, shell=True).stdout

    output = output.readlines()
    files_in_dir = int(output[0].strip())

    expected_files = n_secs * 25
    error = False
    error, n_ics = check_recombine_ics(directory=directory, \
                                           startsec=startsec, n_secs=n_secs, required_size=required_size_ics)
    files_in_dir += n_ics
    if not files_in_dir == expected_files:
        logger.error("We have {0} files but expected {1}".format(files_in_dir, expected_files))
        error = True
    for line in output[1:]:
        if b'dat' in line:
            logger.warning("Deleted {0} due to wrong size.".format(line.strip()))
            error = True
    if not error:
        logger.info("We have all {0} files as expected.".format(files_in_dir))
    return error

def check_recombine_ics(directory=None, startsec=None, n_secs=None, required_size=None, obsID=None):
    if not required_size:
        try:
            _, _, required_size = get_files_and_sizes(obsID, 'ics', mintime=startsec, maxtime=startsec + n_secs)
        except:
            traceback.print_exc()
            return True, 0

    if not startsec:
        #output = subprocess.Popen(["ls -ltr %s/*ics.dat | awk '($5!=%s){print \"file \" $9 \" has size \" $5 \" (expected %s)\"}'" %(directory, required_size, required_size)],
        #                          stdout=subprocess.PIPE, shell=True).communicate()[0]
        command = "ls -l %s/*ics.dat | ((tee /dev/fd/5 | wc -l >/dev/fd/4) 5>&1 | " %(directory) + \
            "awk '($5!=%s){print $9}' | tee >> %s/ics_all.txt | xargs rm -rf) 4>&1;" %(required_size, directory) + \
            "cat %s/ics_all.txt; rm -rf %s/ics_all.txt" %(directory, directory)
        output = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True).stdout
    else:
        n_secs = n_secs if n_secs else 1
        output = subprocess.Popen(["count=0;for sec in `seq -w %s %s `;do let count=${count}+`ls -l %s/*${sec}*ics.dat | " %(startsec, startsec+n_secs-1, directory) + \
                                       "((tee /dev/fd/5 | wc -l >/dev/fd/4) 5>&1 | awk '($5!=%s) " %(required_size) + \
                                       "{print $9}' | tee >> %s/errors_%s.txt | xargs rm -rf) 4>&1`;done;" %(directory,startsec) +\
                                       "echo ${count}; cat %s/errors_%s.txt;rm -rf %s/errors_%s.txt" %(directory,startsec,directory,startsec)],
                                  stdout=subprocess.PIPE, shell=True).stdout

    output = output.readlines()
    files_in_dir = int(output[0].strip())
    error = False
    if not files_in_dir == n_secs:
        logger.error("We have {0} ics-files but expected {1}".format(files_in_dir, n_secs))
        error = True
    for line in output[1:]:
        if b'dat' in line:
            error = True
            line = line.strip().decode()
            logger.error("Deleted {0} due to wrong size.".format(line))
            dat_files = line.replace('_ics.dat','*.dat')
            rm_cmd = "rm -rf {0}".format(dat_files)
            logger.warning("Also running {0} to make sure ics files are rebuilt.".format(rm_cmd))
            subprocess.Popen(rm_cmd, stdout=subprocess.PIPE, shell=True)
    if error == False:
        logger.info("We have all {0} ICS files as expected.".format(files_in_dir))
    return error, files_in_dir


def get_files_and_sizes(obsID, mode, mintime=0, maxtime=2000000000):
    """
    Get files and sizes from the MWA metadata server and check that they're all the same size

    Parameters:
    -----------
    obsID: int
        The MWA observation ID
    mode: str
        The typ of file from 'raw', 'tar_ics' and 'ics'
    mintime: int
        The minimum GPS time of observations to check (inclusive, >=)  Default: 0
    maxtime: int
        The maximum GPS time of observations to check (exculsive, <)  Default: 2000000000

    Returns:
    --------
    files_masked, suffix, sizes[0]: list
        files_masked: list of the files with the input mode/suffix
        suffix:       '.dat', '.tar' or '_ics.dat' depnding on the input mode
        sizes[0]:     size of files in bytes
    """
    if mode == 'raw':
        suffix = '.dat'
    elif mode == 'tar_ics':
        suffix = '.tar'
    elif mode == 'ics':
        suffix = '_ics.dat'
    else:
        logger.error("Wrong mode supplied. Options are raw, tar_ics, and ics")
        return
    logger.info("Retrieving file info from MWA database for all {0} files...".format(suffix))
    files_meta = getmeta(service='data_files', params={'obs_id':obsID, 'nocache':1, 'mintime':mintime, 'maxtime':maxtime})
    # 'nocache' is used above so we get don't use the cached metadata as that could
    # be out of data so we force it to get up to date values
    files = np.array(list(files_meta.keys()))
    files_masked = []
    sizes = []
    for f in files:
        if suffix in f:
            sizes.append(files_meta[f]['size'])
            files_masked.append(f)
    logger.info("...Done. Expect all on database to be {0} bytes in size...".format(sizes[0]))

    size_check = True
    for s in sizes:
        if not s == sizes[0]:
            size_check = False
    if size_check:
        logger.info("...yep they are. Now checking on disk.")
        return files_masked, suffix, sizes[0]
    else:
        logger.error("Not all files have the same size. Check your data!")
        logger.error("{0}".format(np.vstack((files,sizes)).T))
        return