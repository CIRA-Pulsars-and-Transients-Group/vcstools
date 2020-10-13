/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include "../make_beam/beam_common.h"
#include "mwa_hyperbeam.h"
#include "complex.h"

#define MWA_LAT -26.703319             /* Array latitude. degrees North */
#define DD2R      1.74532925199433e-02 /* = PI/180 */
#define DPIBY2    1.570796326794897    /* = PI/2   */

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

    // Read in info from metafits file
    struct metafits_info mi;
    get_metafits_info( opts.metafits, &mi, opts.chan_width );

    // Always use the first tile/pol's set of delays and assume no dead
    // dipoles
    unsigned int *delays = (unsigned int *)mi->delays[0];
    double *amps = { 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0 };

    // Load the FEE2016 beam model
    FEEBeam *beam = new_fee_beam( HYPERBEAM_HDF5 );

    // Open files for reading and writing
    FILE *fP = fopen( opts.pointings, "r" );
    FILE *fF = fopen( "sky_model_comparison_FEE2016.dat", "w" );
    FILE *fA = fopen( "sky_model_comparison_ANALYTIC.dat", "w" );

    if (fF == NULL || fA == NULL)
    {
        fprintf( stderr, "Unable to open file(s) for writing\n" );
        exit(EXIT_FAILURE);
    }

    complex double  JA[4]; // Beam model Jones -- ANALYTIC
    double         *JF;    // Beam model Jones -- FEE2016

    double az, za; // (az)imuth & (z)enith (a)ngle in radians

    int i; // Generic loop counter

    while (fscanf( fP, "%lf %lf", &az, &za ) != EOF)
    {
        // Calculate jones matrix for the ANALYTIC beam
        calcEjones_analytic( JA,                      // pointer to 4-element (2x2) voltage gain Jones matrix
                opts.frequency,                       // observing freq of fine channel (Hz)
                (MWA_LAT*DD2R),                       // observing latitude (radians)
                mi.tile_pointing_az*DD2R,             // azimuth & zenith angle of tile pointing
                (DPIBY2-(mi.tile_pointing_el*DD2R)),
                az, za );                             // azimuth & zenith angle of pencil beam

        // Calculate jones matrix for the ANALYTIC beam
        JF = calc_jones( beam, az, za, opts.frequency, delays, amps, zenith_norm );

        // Write out the matrices to the respective files
        // First the ANALYTIC...
        for (i = 0; i < 4; i++)
        {
            fprintf( fA, "%lf %lf ", creal( JA[i] ), cimag( JA[i] ) );
        }
        fprintf( fA, "\n" );

        // .. and then the FEE2016
        for (i = 0; i < 8; i++)
        {
            fprintf( fF, "%lf ", JF[i] );
        }
        fprintf( fF, "\n" );

        // calc_jones allocates memory, so need to clean up with every loop
        free( JF );
    }

    // Close files for writing
    fclose(fF);
    fclose(fA);

    // Free up memory
    destroy_complex_weights( complex_weights_array, npointing, nstation, nchan );
    destroy_invJi( invJi, nstation, nchan, npol );

    destroy_metafits_info( &mi );

    free( opts.pointings    );
    free( opts.metafits     );

    // Clean up Hyperbeam
    free_fee_beam( beam );

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

}

