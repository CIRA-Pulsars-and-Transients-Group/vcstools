#!/usr/bin/python
import sys
def options (options):

    print "\noptions:\n"
    print "--mode {0}".format(options.mode)
    print "-B [1/0]\t Submit download jobs to the copyq - at the moment this mode will only download and will perform <NO> subsequent processing [%d] \n" % (opts['batch_download'])
    print "-b:\t GPS/UNIX time of the beginning [%d]]\n" % (opts['begin'])
    print "-c:\t Coarse channel count (how many to process) [%d]\n" % (opts['ncoarse_chan'])
    print "-d:\t Number of parallel downloads to envoke if using '-g' [%d]\n" % (opts['parallel_dl'])
    print "-e:\t GPS/UNIX time of the end [%d]\n" % (opts['end'])
 #   print "-g:\t Get the data? (True/False) add this to get fresh data from the archive [%s]\n" % (opts['get_data'])
    print "-i:\t Increment in seconds (how much we process at once) [%d]\n" % (opts['inc'])
    print "-j:\t [corrdir] Use Jones matrices from the RTS [%s,%s]\n" % (opts['useJones'],opts['corrdir'])
    print "-m:\t Beam forming mode (0 == NO BEAMFORMING 1==PSRFITS, 2==VDIF) [%d]\n" % (opts['mode'])
    print "-n:\t Number of fine channels per coarse channel [%d]\n" % (opts['nchan'])
    print "-o:\t obsid [%s]\n" % opts['obsid']
    print "-p:\t beam pointing [%s]\n" % opts['pointing']
    print "-s:\t single step (only process one increment and this is it (-1 == do them all) [%d]\n" % opts['single_step']
#    print "-r:\t [corrdir] Run the offline correlator - this will submit a job to process the .dat files into visibility sets into the specified directory. These are needed if you want an RTS calibration solution [%s]\n" % opts['corrdir']
    print "-G:\t Submit the beamformer/correlator job [Do it = %s]\n" % opts['Go']
#   print "-R:\t New VCS mode - requires the recombine operation [runRECOMBINE = %s]\n" % opts['runRECOMBINE']
    print "-w:\t Working root directory [%s]\n" % opts['root']
#    print "-z:\t Add to switch off PFB formation/testing [runPFB = %s]\n" % opts['runPFB']


def vcs_download():
    print "Downloading files from archive"


def vcs_recombine():
    print "Running recombine on files"


def vcs_correlate():
    print "Correlating files"


def make_pfb_files():
    print "Creating PFB files"


def coherent_beam():
    print "Forming coherent beam"



if __name__ == '__main__':

    modes=['download','recombine','correlate','make_pfb','beamform']
    jobs_per_node = 8
    chan_list_full=["ch01","ch02","ch03","ch04","ch05","ch06","ch07","ch08","ch09","ch10","ch11","ch12","ch13","ch14","ch15","ch16","ch17","ch18","ch19","ch20","ch21","ch22","ch23","ch24"]
    chan_list = []


    from optparse import OptionParser, OptionGroup

    parser=OptionParser(description="Use this damn thing already!")

    group_download = OptionGroup(parser, 'Download Options')
    group_download.add_option("-B", "--copyq", action="store_true", default=False, help="Submit download jobs to the copyq [default=%default]")
    group_download.add_option("--format", type="choice", choices=['11','12'], default='11', help="Voltage data type (Raw = 11, Recombined Raw = 12) [default=%default]")
    group_download.add_option("-d", "--parallel_dl", type="int", default=3, help="Number of parallel downloads to envoke [default=%default]")

    group_recombine = OptionGroup(parser, 'Recombine Options')

    group_correlate = OptionGroup(parser, 'Correlator Options')
    group_correlate.add_option("--ft_res", metavar="FREQ RES,TIME RES", type="int", nargs=2, default=(40,1), help="Frequency (kHz) and Time (s) resolution to run correlator at. [default=%default]")

    group_pfb = OptionGroup(parser, 'PFB Creation Options')

    group_beamform = OptionGroup(parser, 'Beamforming Options')
    group_beamform.add_option("-p", "--pointing", nargs=2, help="R.A. and Dec. of pointing")
    group_beamform.add_option("--bf_mode", type="choice", choices=['0','1','2'], help="Beam forming mode (0 == NO BEAMFORMING 1==PSRFITS, 2==VDIF)")
    group_beamform.add_option("-j", "--useJones", action="store_true", default=False, help="Use Jones matrices from the RTS [default=%default]")

    parser.add_option("-m", "--mode", type="choice", choices=['download','recombine','correlate','make_pfb','beamform'], help="Mode you want to run. {0}".format(modes))
    parser.add_option("-o", "--obs", metavar="OBS ID", type="int", help="Observation ID you want to process [no default]")
    parser.add_option("-b", "--begin", type="int", help="First GPS time to process [no default]")
    parser.add_option("-e", "--end", type="int", help="Last GPS time to process [no default]")
    parser.add_option("-i", "--increment", type="int", default=200, help="Increment in seconds (how much we process at once) [default=%default]")
    parser.add_option("-s", action="store_true", default=False, help="Single step (only process one increment and this is it (False == do them all) [default=%default]")
    parser.add_option("-w", "--work_dir", metavar="DIR", default="/scratch/mwaops/vcs/", help="Base directory you want to run from. This will create a folder for the Obs. ID if it doesn't exist [default=%default]")
    parser.add_option("-c", "--ncoarse_chan", type="int", default=24, help="Coarse channel count (how many to process) [default=%default]")
    parser.add_option("-n", "--nfine_chan", type="int", default=128, help="Number of fine channels per coarse channel [default=%default]")
    parser.add_option("-G", "--Go", action="store_true", default=False, help="Include this option to run script [default=%default]")
    parser.add_option_group(group_download)
#    parser.add_option_group(group_recombine)
    parser.add_option_group(group_correlate)
#   parser.add_option_group(group_pfb)
    parser.add_option_group(group_beamform)

    (opts, args) = parser.parse_args()

    if not opts.mode:
      print "Mode required {0}. Please specify with -m or --mode.".format(modes)
      quit()

    if not opts.obs:
        print "Observation ID required, please put in with -o or --obs"
        quit()

 #   options(opts)

    if opts.mode == 'download':
        print opts.mode
        vcs_download()
    elif opts.mode == 'recombine':
        print opts.mode
        vcs_recombine()
    elif opts.mode == 'correlate':
        print opts.mode
        vcs_correlate()
    elif opts.mode == 'make_pfb':
        print opts.mode
        make_pfb_files()
    elif opts.mode == 'beamformer':
        print opts.mode
        coherent_beam()
    else:
        print "Somehow your non-standard mode snuck through. Try again with one of {0}".format(modes)
        quit()


