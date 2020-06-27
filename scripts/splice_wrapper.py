#! /usr/bin/env python3

import subprocess
import os
import argparse
import time
import mwa_metadb_utils as meta
import glob
from shutil import copy
from itertools import groupby, count

parser = argparse.ArgumentParser(description="""
Wraps the splice_psrfits.sh script to automate it. Should be run from the foulder containing the files.
""")
parser.add_argument('-o','--observation',type=str,help='The observation ID to be used.')
parser.add_argument('-i','--incoh',action='store_true',help='Use this check if there are and incoh files from the beamformer.')
parser.add_argument('-d','--delete',action='store_true',help='This will cause the script to remove the unspliced files if splice_psrfits.sh succeeds (error code of 0).')
parser.add_argument('-w','--work_dir',type=str,help='Working directory of the vcs files.', default="./")
parser.add_argument('-c', '--channels', type=int, nargs=24, help='A list of the observations channel IDs for example "-c 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132". If this option is not used a metadata call will be used to find the channel IDs.')
args=parser.parse_args()


obsid = args.observation

# Check if already spliced
if glob.glob('{0}/{1}*fits'.format(args.work_dir, args.observation)) and \
   not glob.glob('{0}/*_{1}*fits'.format(args.work_dir, args.observation)):
    print('All files are already spliced so exiting')
    exit()

# Get frequency channels
if args.channels:
    channels = args.channels
else:
    print("Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obsid))
    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
print("Chan order: {}".format(channels))


# Move into working dir
old_dir = os.getcwd()
if args.work_dir:
    os.chdir(args.work_dir)


# Getting number of files list
n_fits_list = []
fits_len = []
for ch in channels:
    if args.incoh:
        fits_glob = glob.glob('{}/*{}_incoh_ch{:03d}_*fits'.format(args.work_dir, obsid, ch))
    else:
        fits_glob = glob.glob('{}/*{}*_ch{:03d}_*fits'.format(args.work_dir, obsid, ch))
    n_fits = []
    for file_name in fits_glob:
        n_fits.append(int(file_name[-9:-5]))
    n_fits.sort()
    n_fits_list.append(n_fits)
    fits_len.append(len(n_fits))
print('Total fits number: {}'.format(max(fits_len)))

# Grab pointing from unspliced files
if args.incoh:
    pointing = 'incoh'
else:
    pointing = "{0}_{1}".format(fits_glob[0].split("_")[2], fits_glob[0].split("_")[3])

# Check all files are there
max_fits_order = n_fits_list[fits_len.index(max(fits_len))]
missing_files = []
for ch in channels:
    for file_num in max_fits_order:
        if args.incoh:
            if not glob.glob('{}/*{}_incoh_ch{:03d}_{:04d}.fits'.format(args.work_dir, obsid, ch, file_num)):
                missing_files.append('{}/*{}_incoh_ch{:03d}_{:04d}.fits'.format(args.work_dir, obsid, ch, file_num))
        else:
            if not glob.glob('{}/*{}*_ch{:03d}_{:04d}.fits'.format(args.work_dir, obsid, ch, file_num)):
                missing_files.append('{}/*{}*_ch{:03d}_{:04d}.fits'.format(args.work_dir, obsid, ch, file_num))
if len(missing_files) != 0:
    print("Missing the following files:")
    for mf in missing_files:
        print(mf)
    print("Exiting")
    sys.exit(1)

print('Fits number order: {}'.format(max_fits_order))

# Put channels into consecutive lists to properly splice picket fence data
con_channels = []
for _, g in groupby(channels, key=lambda n, c=count(): n-next(c)):
    con_channels.append(list(g))

for n in max_fits_order:
    for channels in con_channels:
        # List unspliced files
        unspliced_files = []
        for ch in channels:
            if args.incoh:
                unspliced_files.append(glob.glob('{}/*{}_incoh_ch{:03d}_{:04d}.fits'.format(args.work_dir,
                                                    obsid, ch, n))[0])
            else:
                unspliced_files.append(glob.glob('{}/*{}*_ch{:03d}_{:04d}.fits'.format(args.work_dir,
                                                    obsid, ch, n))[0])

        # Create splice command and submit
        submit_line = 'splice_psrfits '
        for us_file in unspliced_files:
            submit_line += '{} '.format(us_file)
        submit_line += 'temp_{}'.format(n)

        print(submit_line)
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        out_lines = submit_cmd.stdout
        for l in out_lines:
            print(l.decode()[:-1])
        time.sleep(1)


        print('Finished {}_{}_ch{:03d}-{:03d}_{:04d}.fits'.format(obsid, pointing, min(channels), max(channels), n))
        os.rename('temp_{}_0001.fits'.format(n),
                  '{}_{}_ch{:03d}-{:03d}_{:04d}.fits'.format(obsid, pointing, min(channels), max(channels), n))

        #wait to get error code
        (output, err) = submit_cmd.communicate()
        p_status = submit_cmd.wait()

        print("exit code: " + str(submit_cmd.returncode))
        if args.delete and int(submit_cmd.returncode) == 0:
            for us_file in unspliced_files:
                print("Deleting: " + str(us_file))
                if os.path.isfile(us_file):
                    os.remove(us_file)
