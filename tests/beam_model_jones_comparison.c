/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

// TODO: Remove superfluous #includes
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include "star/pal.h"
#include "star/palmac.h"
#include "../make_beam/ascii_header.h"
#include "../make_beam/mwa_header.h"
#include <glob.h>
#include <fcntl.h>
#include <assert.h>
#include "../make_beam/beam_common.h"
#include "../make_beam/beam_psrfits.h"
#include "../make_beam/beam_vdif.h"
#include "../make_beam/make_beam.h"
#include "../make_beam/vdifio.h"
#include "../make_beam/filter.h"
#include "psrfits.h"
#include "../make_beam/mycomplex.h"
#include "../make_beam/form_beam.h"
#include <omp.h>


int main(int argc, char **argv)
{
    #ifndef HAVE_CUDA
    // Initialise FFTW with OpenMP
    fftw_init_threads();
    fftw_plan_with_nthreads( omp_get_max_threads() );
    #endif

    // A place to hold the beamformer settings
    struct make_beam_opts opts;

    /* Set default beamformer settings */

    // Variables for required options
    opts.obsid       = NULL; // The observation ID
    opts.begin       = 0;    // GPS time -- when to start beamforming
    opts.end         = 0;    // GPS time -- when to stop beamforming
    opts.time_utc    = NULL; // utc time string "yyyy-mm-ddThh:mm:ss"
    opts.pointings   = NULL; // list of pointings "dd:mm:ss_hh:mm:ss,dd:mm:ss_hh:mm:ss"
    opts.datadir     = NULL; // The path to where the recombined data live
    opts.metafits    = NULL; // filename of the metafits file
    opts.rec_channel = NULL; // 0 - 255 receiver 1.28MHz channel
    opts.frequency   = 0;    // = rec_channel expressed in Hz
#ifdef HYPERBEAM_HDF5
    opts.beam_model  = BEAM_FEE2016;
#else
    opts.beam_model  = BEAM_ANALYTIC;
#endif

    // Variables for MWA/VCS configuration
    opts.nstation      = 128;    // The number of antennas
    opts.nchan         = 128;    // The number of fine channels (per coarse channel)
    opts.chan_width    = 10000;  // The bandwidth of an individual fine chanel (Hz)
    opts.sample_rate   = 10000;  // The VCS sample rate (Hz)
    opts.custom_flags  = NULL;   // Use custom list for flagging antennas

    // Output options
    opts.out_incoh     = 0;  // Default = PSRFITS (incoherent) output turned OFF
    opts.out_coh       = 0;  // Default = PSRFITS (coherent)   output turned OFF
    opts.out_vdif      = 0;  // Default = VDIF                 output turned OFF
    opts.out_bf        = 1;  // Default = beamform all (non-flagged) antennas
    opts.out_ant       = 0;  // The antenna number (0-127) to write out if out_bf = 0
    opts.synth_filter  = NULL;
    opts.out_summed    = 0;  // Default = output only Stokes I output turned OFF
    opts.max_sec_per_file = 200; // Number of seconds per fits files

    // Variables for calibration settings
    opts.cal.filename          = NULL;
    opts.cal.bandpass_filename = NULL;
    opts.cal.chan_width        = 40000;
    opts.cal.nchan             = 0;
    opts.cal.cal_type          = NO_CALIBRATION;
    opts.cal.offr_chan_num     = 0;

    // GPU options
    opts.gpu_mem               = -1.0;

    // Parse command line arguments
    make_beam_parse_cmdline( argc, argv, &opts );

    // Create "shorthand" variables for options that are used frequently
    int nstation             = opts.nstation;
    int nchan                = opts.nchan;
    const int npol           = 2;   // (X,Y)

    //float vgain = 1.0; // This is re-calculated every second for the VDIF output


    // Calculate the number of files
    int nfiles = opts.end - opts.begin + 1;
    if (nfiles <= 0) {
        fprintf(stderr, "Cannot beamform on %d files (between %lu and %lu)\n", nfiles, opts.begin, opts.end);
        exit(EXIT_FAILURE);
    }

    // Parse input pointings
    int max_npointing = 120; // Could be more
    char RAs[max_npointing][64];
    char DECs[max_npointing][64];
    int npointing = sscanf( opts.pointings,
            "%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,],%[^_]_%[^,]," ,
                            RAs[0],  DECs[0],  RAs[1],  DECs[1],  RAs[2],  DECs[2],
                            RAs[3],  DECs[3],  RAs[4],  DECs[4],  RAs[5],  DECs[5],
                            RAs[6],  DECs[6],  RAs[7],  DECs[7],  RAs[8],  DECs[8],
                            RAs[9],  DECs[9],  RAs[10], DECs[10], RAs[11], DECs[11],
                            RAs[12], DECs[12], RAs[13], DECs[13], RAs[14], DECs[14],
                            RAs[15], DECs[15], RAs[16], DECs[16], RAs[17], DECs[17],
                            RAs[18], DECs[18], RAs[19], DECs[19], RAs[20], DECs[20],
                            RAs[21], DECs[21], RAs[22], DECs[22], RAs[23], DECs[23],
                            RAs[24], DECs[24], RAs[25], DECs[25], RAs[26], DECs[26],
                            RAs[27], DECs[27], RAs[28], DECs[28], RAs[29], DECs[29],
                            RAs[30], DECs[30], RAs[31], DECs[31], RAs[32], DECs[32],
                            RAs[33], DECs[33], RAs[34], DECs[34], RAs[35], DECs[35],
                            RAs[36], DECs[36], RAs[37], DECs[37], RAs[38], DECs[38],
                            RAs[39], DECs[39], RAs[40], DECs[40], RAs[41], DECs[41],
                            RAs[42], DECs[42], RAs[43], DECs[43], RAs[44], DECs[44],
                            RAs[45], DECs[45], RAs[46], DECs[46], RAs[47], DECs[47],
                            RAs[48], DECs[48], RAs[49], DECs[49], RAs[50], DECs[50],
                            RAs[51], DECs[51], RAs[52], DECs[52], RAs[53], DECs[53],
                            RAs[54], DECs[54], RAs[55], DECs[55], RAs[56], DECs[56],
                            RAs[57], DECs[57], RAs[58], DECs[58], RAs[59], DECs[59],
                            RAs[60], DECs[60], RAs[61], DECs[61], RAs[62], DECs[62],
                            RAs[63], DECs[63], RAs[64], DECs[64], RAs[65], DECs[65],
                            RAs[66], DECs[66], RAs[67], DECs[67], RAs[68], DECs[68],
                            RAs[69], DECs[69], RAs[70], DECs[70], RAs[71], DECs[71],
                            RAs[72], DECs[72], RAs[73], DECs[73], RAs[74], DECs[74],
                            RAs[75], DECs[75], RAs[76], DECs[76], RAs[77], DECs[77],
                            RAs[78], DECs[78], RAs[79], DECs[79], RAs[80], DECs[80],
                            RAs[81], DECs[81], RAs[82], DECs[82], RAs[83], DECs[83],
                            RAs[84], DECs[84], RAs[85], DECs[85], RAs[86], DECs[86],
                            RAs[87], DECs[87], RAs[88], DECs[88], RAs[89], DECs[89],
                            RAs[90], DECs[90], RAs[91], DECs[91], RAs[92], DECs[92],
                            RAs[93], DECs[93], RAs[94], DECs[94], RAs[95], DECs[95],
                            RAs[96], DECs[96], RAs[97], DECs[97], RAs[98], DECs[98],
                            RAs[99], DECs[99], RAs[100], DECs[100], RAs[101], DECs[101],
                            RAs[102], DECs[102], RAs[103], DECs[103], RAs[104], DECs[104],
                            RAs[105], DECs[105], RAs[106], DECs[106], RAs[107], DECs[107],
                            RAs[108], DECs[108], RAs[109], DECs[109], RAs[110], DECs[110],
                            RAs[111], DECs[111], RAs[112], DECs[112], RAs[113], DECs[113],
                            RAs[114], DECs[114], RAs[115], DECs[115], RAs[116], DECs[116],
                            RAs[117], DECs[117], RAs[118], DECs[118], RAs[119], DECs[119] );

    if (npointing%2 == 1)
    {
        fprintf(stderr, "Number of RAs do not equal the number of Decs given. Exiting\n");
        fprintf(stderr, "npointings : %d\n", npointing);
        fprintf(stderr, "RAs[0] : %s\n", RAs[0]);
        fprintf(stderr, "DECs[0] : %s\n", DECs[0]);
        exit(0);
    }
    else
        npointing /= 2; // converting from number of RAs and DECs to number of pointings

    char pointing_array[npointing][2][64];
    int p;
    for ( p = 0; p < npointing; p++)
    {
       strcpy( pointing_array[p][0], RAs[p] );
       strcpy( pointing_array[p][1], DECs[p] );
    }

    // Allocate memory
    ComplexDouble  ****complex_weights_array = create_complex_weights( npointing, nstation, nchan, npol );
    ComplexDouble  ****invJi                 = create_invJi( nstation, nchan, npol );

    // Read in info from metafits file
    struct metafits_info mi;
    get_metafits_info( opts.metafits, &mi, opts.chan_width );

    // If using bandpass calibration solutions, calculate number of expected bandpass channels
    if (opts.cal.cal_type == RTS_BANDPASS)
        opts.cal.nchan = (nchan * opts.chan_width) / opts.cal.chan_width;

    // If a custom flag file has been provided, use that instead of the metafits flags
    int i;
    if (opts.custom_flags != NULL)
    {
        // Reset the weights to 1
        for (i = 0; i < nstation*npol; i++)
            mi.weights_array[i] = 1.0;

        // Open custom flag file for reading
        FILE *flagfile = fopen( opts.custom_flags, "r" );
        if (flagfile == NULL)
        {
            fprintf( stderr, "error: couldn't open flag file \"%s\" for "
                             "reading\n", opts.custom_flags );
            exit(EXIT_FAILURE);
        }

        // Read in flags from file
        int nitems;
        int flag, ant;
        while (!feof(flagfile))
        {
            // Read in next item
            nitems = fscanf( flagfile, "%d", &ant );
            if (nitems != 1 && !feof(flagfile))
            {
                fprintf( stderr, "error: couldn't parse flag file \"%s\"\n",
                        opts.custom_flags );
                exit(EXIT_FAILURE);
            }

            // Flag both polarisations of the antenna in question
            flag = ant*2;
            mi.weights_array[flag]   = 0.0;
            mi.weights_array[flag+1] = 0.0;
        }

        // Close file
        fclose( flagfile );
    }

    // Issue warnings if any antennas are being used which are flagged in the metafits file
    for (i = 0; i < nstation*npol; i++)
    {
        if (mi.weights_array[i] != 0.0 &&
            mi.flag_array[i]    != 0.0)
        {
            fprintf( stderr, "warning: antenna %3d, pol %d is included even "
                             "though it is flagged in the metafits file\n",
                             i / npol,
                             i % npol );
        }
    }

    int j, file_no, ant, ch, xpol, ypol;
    FEEBeam *beam = NULL;
    char comp_fname[1024];
    for (j=0; j<2; j++)
    {
        int beam_model;
        if (j==0)
        {
            beam_model=BEAM_FEE2016;
            sprintf(comp_fname, "beam_model_comparison_FEE2016.dat");
        }
        else
        {
            beam_model=BEAM_ANALYTIC;
            sprintf(comp_fname, "beam_model_comparison_ANALYTIC.dat");
        }

        // Load the FEE2016 beam, if requested
        if (opts.beam_model == BEAM_FEE2016) {
            new_fee_beam( HYPERBEAM_HDF5, &beam );
        }
        FILE *f = fopen(comp_fname, "w");
        for (file_no = 0; file_no < nfiles; file_no++)
        {
            // Read in data from next file
            // Get the next second's worth of phases / jones matrices, if needed
            get_delays(
                    pointing_array,         // an array of pointings [pointing][ra/dec][characters]
                    npointing,              // number of pointings
                    opts.frequency,         // middle of the first frequency channel in Hz
                    &opts.cal,              // struct holding info about calibration
                    opts.sample_rate,       // = 10000 samples per sec
                    beam_model,             // beam model type
                    beam,                   // Hyperbeam struct
                    opts.time_utc,          // utc time string
                    (double)file_no,        // seconds offset from time_utc at which to calculate delays
                    NULL,                   // Don't update delay_vals
                    &mi,                    // Struct containing info from metafits file
                    complex_weights_array,  // complex weights array (answer will be output here)
                    invJi );                // invJi array           (answer will be output here)



            for (ant = 0; ant < nstation; ant++)
            {
                fprintf(f, "# File number %d Antenna %d\n", file_no, ant);
                for (ch = 0; ch < nchan; ch++)
                {
                    for (xpol = 0; xpol < npol; xpol++)
                    {
                        for (ypol = 0; ypol < npol; ypol++)
                        {
                            fprintf(f, "%lf %lf ", CReald(invJi[ant][ch][xpol][ypol]), CImagd(invJi[ant][ch][xpol][ypol]));
                        }
                    }
                    fprintf(f, "\n");
                }
                fprintf(f, "\n\n");
            }

        }
        fclose(f);
    }
    // Free up memory
    destroy_complex_weights( complex_weights_array, npointing, nstation, nchan );
    destroy_invJi( invJi, nstation, nchan, npol );

    destroy_metafits_info( &mi );

    free( opts.obsid        );
    free( opts.time_utc     );
    free( opts.pointings    );
    free( opts.datadir      );
    free( opts.metafits     );
    free( opts.rec_channel  );
    free( opts.cal.filename );
    free( opts.custom_flags );
    free( opts.synth_filter );

    // Clean up Hyperbeam
    if (opts.beam_model == BEAM_FEE2016) {
        free_fee_beam( beam );
    }

#ifndef HAVE_CUDA
    // Clean up FFTW OpenMP
    fftw_cleanup_threads();
#endif

    return 0;
}


void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: make_beam [OPTIONS]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-o, --obsid=GPSTIME       ");
    fprintf(stderr, "Observation ID (GPS seconds).\n");
    fprintf(stderr, "\t-b, --begin=GPSTIME       ");
    fprintf(stderr, "Begin time of observation, in GPS seconds\n");
    fprintf(stderr, "\t-e, --end=GPSTIME         ");
    fprintf(stderr, "End time of observation, in GPS seconds\n");
    fprintf(stderr, "\t-z, --utc-time=UTCTIME    ");
    fprintf(stderr, "The UTC time that corresponds to the GPS time given by the -b\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "option. UTCTIME must have the format: yyyy-mm-ddThh:mm:ss\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-P, --pointings=hh:mm:ss.s_dd:mm:ss.s,hh:mm:ss.s_dd:mm:ss.s...\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "Right ascension and declinations of multiple pointings\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-m, --metafits-file=FILE  ");
    fprintf(stderr, "FILE is the metafits file pertaining to the OBSID given by the\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr,  "-o option\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-f, --coarse-chan=N       ");
    fprintf(stderr, "Absolute coarse channel number (0-255)\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "MWA/VCS CONFIGURATION OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-a, --antennas=N          ");
    fprintf(stderr, "The number of antennas in the array. For MWA Phase 2, N=128.     ");
    fprintf(stderr, "[default: 128]\n");
    fprintf(stderr, "\t-n, --num-fine-chans=N    ");
    fprintf(stderr, "The number of fine channels per coarse channel.                  ");
    fprintf(stderr, "[default: 128]\n");
    fprintf(stderr, "\t-w, --fine-chan-width=N   ");
    fprintf(stderr, "The bandwidth of an individual fine channel (Hz).                ");
    fprintf(stderr, "[default: 10000]\n");
    fprintf(stderr, "\t-r, --sample-rate=N       ");
    fprintf(stderr, "The VCS sample rate, in Hz. (The sample rate given in the meta-  ");
    fprintf(stderr, "[default: 10000]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "metafits file matches the correlator settings at the time of\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the observation, which is not necessarily the same as that of\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the VCS. Hence the necessity of this option.)\n");
    fprintf(stderr, "\t-F, --custom-flags=file   ");
    fprintf(stderr, "Flag the antennas listed in file instead of those flagged in the ");
    fprintf(stderr, "[default: none]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "metafits file given by the -m option.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "CALIBRATION OPTIONS (RTS)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-J, --dijones-file=PATH   ");
    fprintf(stderr, "The direction-independent Jones matrix file that is output from\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "the RTS. Using this option instructs the beamformer to use the\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "RTS-generated calibration solution. Either -J or -O must be\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "supplied. If both are supplied the one that comes last will\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "override the former.\n");
    fprintf(stderr, "\t-B, --bandpass-file=PATH  ");
    fprintf(stderr, "The bandpass file that is output from the RTS. If this option\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "is given, the RTS calibration solution will be applied to each\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "fine channel. If -J is supplied but -B is not, then the coarse\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "channel solution will be applied to ALL fine channels\n");
    fprintf(stderr, "\t-W, --rts-chan-width      ");
    fprintf(stderr, "RTS calibration channel bandwidth (Hz)                           ");
    fprintf(stderr, "[default: 40000]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "CALIBRATION OPTIONS (OFFRINGA)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-O, --offringa-file=PATH  ");
    fprintf(stderr, "The calibration solution file that is output from the tools\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "made by Andre Offringa. Using this option instructs the beam-\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "former to use the Offringa-style calibration solution. Either\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "-J or -O must be supplied. If both are supplied the one that\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "comes last will override the former.\n");
    fprintf(stderr, "\t-C, --offringa-chan=N     ");
    fprintf(stderr, "The zero-offset position of the coarse channel solution in the   ");
    fprintf(stderr, "[default: 0]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "calibration file given by the -O option.\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OTHER OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-h, --help                ");
    fprintf(stderr, "Print this help and exit\n");
    fprintf(stderr, "\t-g, --gpu-mem=N     ");
    fprintf(stderr, "The maximum amount of GPU memory you want make_beam to use in GB ");
    fprintf(stderr, "[default: -1]\n");
    fprintf(stderr, "\t-V, --version             ");
    fprintf(stderr, "Print version number and exit\n");
    fprintf(stderr, "\n");
}



void make_beam_parse_cmdline(
        int argc, char **argv, struct make_beam_opts *opts )
{
    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"obsid",           required_argument, 0, 'o'},
                {"begin",           required_argument, 0, 'b'},
                {"end",             required_argument, 0, 'e'},
                {"utc-time",        required_argument, 0, 'z'},
                {"pointings",       required_argument, 0, 'P'},
                {"metafits-file",   required_argument, 0, 'm'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"antennas",        required_argument, 0, 'a'},
                {"num-fine-chans",  required_argument, 0, 'n'},
                {"fine-chan-width", required_argument, 0, 'w'},
                {"sample-rate",     required_argument, 0, 'r'},
                {"custom-flags",    required_argument, 0, 'F'},
                {"dijones-file",    required_argument, 0, 'J'},
                {"bandpass-file",   required_argument, 0, 'B'},
                {"rts-chan-width",  required_argument, 0, 'W'},
                {"offringa-file",   required_argument, 0, 'O'},
                {"offringa-chan",   required_argument, 0, 'C'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "a:A:b:B:C:d:e:f:F:g:hHiJ:m:n:o:O:pP:r:sS:t:vVw:W:z:",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'a':
                    opts->nstation = atoi(optarg);
                    break;
                case 'A':
                    opts->out_bf = 0; // Turn off normal beamforming
                    opts->out_ant = atoi(optarg); // 0-127
                    break;
                case 'b':
                    opts->begin = atol(optarg);
                    break;
                case 'B':
                    opts->cal.bandpass_filename = strdup(optarg);
                    opts->cal.cal_type = RTS_BANDPASS;
                    break;
                case 'C':
                    opts->cal.offr_chan_num = atoi(optarg);
                    break;
                case 'd':
                    opts->datadir = strdup(optarg);
                    break;
                case 'e':
                    opts->end = atol(optarg);
                    break;
                case 'f':
                    opts->rec_channel = strdup(optarg);
                    // The base frequency of the coarse channel in Hz
                    opts->frequency = atoi(optarg) * 1.28e6 - 640e3;
                    break;
                case 'F':
                    opts->custom_flags = strdup(optarg);
                    break;
                case 'g':
                    opts->gpu_mem = atof(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'i':
                    opts->out_incoh = 1;
                    break;
                case 'J':
                    opts->cal.filename = strdup(optarg);
                    if (opts->cal.cal_type != RTS_BANDPASS)
                        opts->cal.cal_type = RTS;
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'n':
                    opts->nchan = atoi(optarg);
                    break;
                case 'o':
                    opts->obsid = strdup(optarg);
                    break;
                case 'O':
                    opts->cal.filename = strdup(optarg);
                    opts->cal.cal_type = OFFRINGA;
                    break;
                case 'p':
                    opts->out_coh = 1;
                    break;
                case 'P':
                    opts->pointings = strdup(optarg);
                    break;
                case 'r':
                    opts->sample_rate = atoi(optarg);
                    break;
                case 'S':
                    opts->synth_filter = strdup(optarg);
                    break;
                case 's':
                    opts->out_summed = 1;
                    break;
                case 't':
                    opts->max_sec_per_file = atoi(optarg);
                    break;
                case 'v':
                    opts->out_vdif = 1;
                    break;
                case 'V':
                    fprintf( stderr, "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(0);
                    break;
                case 'w':
                    opts->chan_width = atoi(optarg);
                    break;
                case 'W':
                    opts->cal.chan_width = atoi(optarg);
                    break;
                case 'z':
                    opts->time_utc = strdup(optarg);
                    break;
                default:
                    fprintf(stderr, "error: make_beam_parse_cmdline: "
                                    "unrecognised option '%s'\n", optarg);
                    usage();
                    exit(EXIT_FAILURE);
            }
        }
    }
    else {
        usage();
        exit(EXIT_FAILURE);
    }


    // Check that all the required options were supplied
    assert( opts->obsid        != NULL );
    assert( opts->begin        != 0    );
    assert( opts->end          != 0    );
    assert( opts->time_utc     != NULL );
    assert( opts->pointings    != NULL );
    assert( opts->datadir      != NULL );
    assert( opts->metafits     != NULL );
    assert( opts->rec_channel  != NULL );
    assert( opts->cal.cal_type != NO_CALIBRATION );

    // If neither -i, -p, nor -v were chosen, set -p by default
    if ( !opts->out_incoh && !opts->out_coh && !opts->out_vdif )
    {
        opts->out_coh = 1;
    }
}


ComplexDouble ****create_complex_weights( int npointing, int nstation, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, ant, ch; // Loop variables
    ComplexDouble ****array;

    array = (ComplexDouble ****)malloc( npointing * sizeof(ComplexDouble ***) );

    for (p = 0; p < npointing; p++)
    {
        array[p] = (ComplexDouble ***)malloc( nstation * sizeof(ComplexDouble **) );

        for (ant = 0; ant < nstation; ant++)
        {
            array[p][ant] = (ComplexDouble **)malloc( nchan * sizeof(ComplexDouble *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][ant][ch] = (ComplexDouble *)malloc( npol * sizeof(ComplexDouble) );
        }
    }
    return array;
}


void destroy_complex_weights( ComplexDouble ****array, int npointing, int nstation, int nchan )
{
    int p, ant, ch;
    for (p = 0; p < npointing; p++)
    {
        for (ant = 0; ant < nstation; ant++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][ant][ch] );

            free( array[p][ant] );
        }
        free( array[p] );
    }
    free( array );
}

ComplexDouble ****create_invJi( int nstation, int nchan, int npol )
// Allocate memory for (inverse) Jones matrices
{
    int ant, pol, ch; // Loop variables
    ComplexDouble ****invJi;
    invJi = (ComplexDouble ****)malloc( nstation * sizeof(ComplexDouble ***) );

    for (ant = 0; ant < nstation; ant++)
    {
        invJi[ant] =(ComplexDouble ***)malloc( nchan * sizeof(ComplexDouble **) );

        for (ch = 0; ch < nchan; ch++)
        {
            invJi[ant][ch] = (ComplexDouble **)malloc( npol * sizeof(ComplexDouble *) );

            for (pol = 0; pol < npol; pol++)
                invJi[ant][ch][pol] = (ComplexDouble *)malloc( npol * sizeof(ComplexDouble) );
        }
    }
    return invJi;
}


void destroy_invJi( ComplexDouble ****array, int nstation, int nchan, int npol )
{
    int ant, ch, pol;
    for (ant = 0; ant < nstation; ant++)
    {
        for (ch = 0; ch < nchan; ch++)
        {
            for (pol = 0; pol < npol; pol++)
                free( array[ant][ch][pol] );

            free( array[ant][ch] );
        }
        free( array[ant] );
    }
    free( array );
}


ComplexDouble ****create_detected_beam( int npointing, int nsamples, int nchan, int npol )
// Allocate memory for complex weights matrices
{
    int p, s, ch; // Loop variables
    ComplexDouble ****array;

    array = (ComplexDouble ****)malloc( npointing * sizeof(ComplexDouble ***) );
    for (p = 0; p < npointing; p++)
    {
        array[p] = (ComplexDouble ***)malloc( nsamples * sizeof(ComplexDouble **) );

        for (s = 0; s < nsamples; s++)
        {
            array[p][s] = (ComplexDouble **)malloc( nchan * sizeof(ComplexDouble *) );

            for (ch = 0; ch < nchan; ch++)
                array[p][s][ch] = (ComplexDouble *)malloc( npol * sizeof(ComplexDouble) );
        }
    }
    return array;
}

void destroy_detected_beam( ComplexDouble ****array, int npointing, int nsamples, int nchan )
{
    int p, s, ch;
    for (p = 0; p < npointing; p++)
    {
        for (s = 0; s < nsamples; s++)
        {
            for (ch = 0; ch < nchan; ch++)
                free( array[p][s][ch] );

            free( array[p][s] );
        }

        free( array[p] );
    }

    free( array );
}
