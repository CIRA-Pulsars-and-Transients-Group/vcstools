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
#include "slalib.h"
#include "slamac.h"
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
    // A place to hold the beamformer settings
    struct make_beam_opts opts;

    /* Set default beamformer settings */

    // Variables for required options
    opts.pointings   = NULL; // File containing az_za pointings
    opts.metafits    = NULL; // filename of the metafits file
    opts.frequency   = 0;    // = rec_channel expressed in Hz

    // Variables for MWA/VCS configuration
    opts.nstation      = 128;    // The number of antennas
    opts.nchan         = 128;    // The number of fine channels (per coarse channel)
    opts.chan_width    = 10000;  // The bandwidth of an individual fine chanel (Hz)
    opts.sample_rate   = 10000;  // The VCS sample rate (Hz)

    // Parse command line arguments
    make_beam_parse_cmdline( argc, argv, &opts );

    // Create "shorthand" variables for options that are used frequently
    int nstation             = opts.nstation;
    int nchan                = opts.nchan;

    // YET TO DO: Prepare pointings file for reading


    // Read in info from metafits file
    struct metafits_info mi;
    get_metafits_info( opts.metafits, &mi, opts.chan_width );

    int j, file_no, ant, ch, xpol, ypol;
    int beam_model;
    FEEBeam *beam = new_fee_beam( HYPERBEAM_HDF5 );
    char comp_fname[1024];
    FILE *fF = fopen("sky_model_comparison_FEE2016.dat", "w");
    FILE *fA = fopen("sky_model_comparison_ANALYTIC.dat", "w");
    if (fF == NULL || fA == NULL)
    {
        fprintf( stderr, "Unable to open file(s) for writing\n" );
        exit(EXIT_FAILURE);
    }

        // Load the FEE2016 beam, if requested
        if (opts.beam_model == BEAM_FEE2016) {
        }
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
    fprintf(stderr, "usage: sky_model_jones_comparison [OPTIONS]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-p, --pointing_file=FILE  ");
    fprintf(stderr, "Two-column ascii file containing az and za (in deg)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-m, --metafits-file=FILE  ");
    fprintf(stderr, "FILE is the metafits file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-f, --frequency=FREQ      ");
    fprintf(stderr, "FREQ is the frequency in Hz");
    fprintf(stderr, "\n");

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

    fprintf(stderr, "\n");
    fprintf(stderr, "OTHER OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-h, --help                ");
    fprintf(stderr, "Print this help and exit\n");
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
                {"antennas",        required_argument, 0, 'a'},
                {"frequency",       required_argument, 0, 'f'},
                {"help",            required_argument, 0, 'h'},
                {"metafits-file",   required_argument, 0, 'm'},
                {"num-fine-chans",  required_argument, 0, 'n'},
                {"pointings_file",  required_argument, 0, 'p'},
                {"version",         required_argument, 0, 'V'},
                {"fine-chan-width", required_argument, 0, 'w'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "a:f:hm:n:p:r:Vw:",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'a':
                    opts->nstation = atoi(optarg);
                    break;
                case 'f':
                    opts->frequency = atoi(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'n':
                    opts->nchan = atoi(optarg);
                    break;
                case 'p':
                    opts->pointings = strdup(optarg);
                    break;
                case 'S':
                    opts->synth_filter = strdup(optarg);
                    break;
                case 'V':
                    fprintf( stderr, "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(0);
                    break;
                case 'w':
                    opts->chan_width = atoi(optarg);
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
    assert( opts->pointings    != NULL );
    assert( opts->metafits     != NULL );
    assert( opts->frequency    != 0 );

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
