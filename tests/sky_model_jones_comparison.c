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
#include "../make_beam/mycomplex.h"

#define DD2R      1.74532925199433e-02 /* = PI/180 */
#define DPIBY2    1.570796326794897    /* = PI/2   */
#define NDELAYS  16

typedef struct sky_opts_t
{
    char     *pointings;       // File containing az_za pointings
    char     *metafits;        // Filename of the metafits file
    long int  frequency;       // Observing frequency (in Hz)
    char     *analytic_output; // Name of analytic output file
    char     *fee2016_output;  // Name of fee2016 output file
    int       print_header;    // boolean: print header in output files?
    int       swap_columns;    // boolean: swap the columns of the FEE2016 Jones matrix
    int       conjugate;       // boolean: conjugate each element of FEE2016 Jones matrix
    int       parallactic;     // boolean: apply parallactic angle correction to the FEE2016 Jones matrix
} sky_opts;

void sky_model_parse_cmdline( int , char **, sky_opts * );
void fprint_header( int, char **, FILE *, unsigned int *, double *,
        sky_opts *, struct metafits_info *, int );
void fprint_jones( FILE *, ComplexDouble * );
void swap_columns( ComplexDouble *, ComplexDouble *);
void conjugate( ComplexDouble *, ComplexDouble *);

int main(int argc, char **argv)
{
    // A place to hold the beamformer settings
    sky_opts opts;

    /* Set default beamformer settings */

    // Variables for required options
    opts.pointings       = NULL; // File containing az_za pointings
    opts.metafits        = NULL; // Filename of the metafits file
    opts.frequency       = 0;    // Observing frequency (in Hz)
    opts.analytic_output = NULL; // Name of analytic output file
    opts.fee2016_output  = NULL; // Name of fee2016 output file
    opts.print_header    = 0;    // boolean: print header in output files?
    opts.swap_columns    = 1;    // boolean: swap the columns of the FEE2016 Jones matrix
    opts.conjugate       = 0;    // boolean: conjugate the FEE2016 Jones matrix
    opts.parallactic     = 1;    // boolean: apply parallactic angle correction to the FEE2016 Jones matrix

    // Variables for MWA/VCS configuration
    int chan_width    = 10000;  // The bandwidth of an individual fine chanel (Hz)

    // Parse command line arguments
    sky_model_parse_cmdline( argc, argv, &opts );

    // Read in info from metafits file
    struct metafits_info mi;
    get_metafits_info( opts.metafits, &mi, chan_width );

    // Always use the first tile/pol's set of delays and assume no dead
    // dipoles
    unsigned int *delays = (unsigned int *)mi.delays[0];
    double amps[] = { 1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0 };

    // Load the FEE2016 beam model
    FEEBeam *beam = new_fee_beam( HYPERBEAM_HDF5 );

    // Open files for reading and writing
    FILE *fP = fopen( opts.pointings,       "r" );
    FILE *fF = fopen( opts.fee2016_output,  "w" );
    FILE *fA = fopen( opts.analytic_output, "w" );

    if (fF == NULL || fA == NULL)
    {
        fprintf( stderr, "Unable to open file(s) for writing: "
                "analytic = \"%s\", fee2016 = \"%s\"\n",
                opts.analytic_output, opts.fee2016_output );
        exit(EXIT_FAILURE);
    }

    ComplexDouble  JA[4]; // Beam model Jones -- ANALYTIC
    double         *JFtmp;
    ComplexDouble  JF[4];    // Beam model Jones -- FEE2016

    // Parallactic angle correction matrix, if needed
    double P[4];

    double az, za; // (az)imuth & (z)enith (a)ngle in radians
    int zenith_norm = 1; // Normalise FEE2016 beam to zenith

    if (opts.print_header)
    {
        fprint_header( argc, argv, fA, delays, amps, &opts, &mi, zenith_norm );
        fprint_header( argc, argv, fF, delays, amps, &opts, &mi, zenith_norm );
    }

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
        JFtmp = calc_jones( beam, az, za, opts.frequency, delays, amps, zenith_norm );
        JF[0] = JFtmp[0] + JFtmp[1]*I;
        JF[1] = JFtmp[2] + JFtmp[3]*I;
        JF[2] = JFtmp[4] + JFtmp[5]*I;
        JF[3] = JFtmp[6] + JFtmp[7]*I;

        // Make alterations to FEE beam as requested
        if (opts.parallactic)
        {
            parallactic_angle_correction_fee2016(
                    P,                  // output = rotation matrix
                    (MWA_LAT*DD2R),     // observing latitude (radians)
                    az, za );
            mult2x2d_RxC( P, JF, JF );
        }

        if (opts.swap_columns)
            swap_columns( JF, JF );

        if (opts.conjugate)
            conj2x2( JF, JF );

        // Write out the matrices to the respective files
        // First the ANALYTIC...
        fprint_jones( fA, JA );
        fprintf( fA, "\n" );

        // .. and then the FEE2016
        fprint_jones( fF, JF );
        fprintf( fF, "\n" );

        // calc_jones allocates memory, so need to clean up with every loop
        free( JFtmp );
    }

    // Close files for writing
    fclose(fF);
    fclose(fA);

    // Free up memory
    destroy_metafits_info( &mi );

    free( opts.pointings    );
    free( opts.metafits     );

    free( opts.analytic_output );
    free( opts.fee2016_output );

    // Clean up Hyperbeam
    free_fee_beam( beam );

    return 0;
}


void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: sky_model_jones_comparison [OPTIONS]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "DESCRIPTION\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tCalculates the Jones matrices for both the original "
                      "analytic beam and the FEE2016 beam for a tile\n"
                    "\tpointing given in the supplied metafits file and a "
                      "list of tied-array beam pointings supplied in\n"
                    "\tthe supplied two-column ASCII file. Two files are "
                      "generated: sky_model_comparison_ANALYTIC.dat\n"
                    "\tand sky_model_comparison_FEE2016.dat.\n" );

    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-p, --pointing_file=FILE  ");
    fprintf(stderr, "Two-column ascii file containing az and za (in rad)\n");
    fprintf(stderr, "\t-m, --metafits-file=FILE  ");
    fprintf(stderr, "FILE is the metafits file\n");
    fprintf(stderr, "\t-f, --frequency=FREQ      ");
    fprintf(stderr, "FREQ is the frequency in Hz");
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OTHER OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-A, --analytic_output     ");
    fprintf(stderr, "Name of analytic output file (default: sky_model_comparison_ANALYTIC.dat\n");
    fprintf(stderr, "\t-c, --conjugate          ");
    fprintf(stderr, "Conjugate the elements of the FEE2016 Jones matrix\n");
    fprintf(stderr, "\t-F, --fee2016_output      ");
    fprintf(stderr, "Name of fee2016 output file (default: sky_model_comparison_FEE2016.dat\n");
    fprintf(stderr, "\t-H, --header              ");
    fprintf(stderr, "Include header in output files (default: no header)\n");
    fprintf(stderr, "\t-h, --help                ");
    fprintf(stderr, "Print this help and exit\n");
    fprintf(stderr, "\t-P, --no_parallactic      ");
    fprintf(stderr, "Turn OFF the parallactic angle correction for the FEE2016 beam\n");
    fprintf(stderr, "\t-s, --no_swap_columns     ");
    fprintf(stderr, "Turn OFF column swap for the FEE2016 Jones matrix\n");
    fprintf(stderr, "\t-V, --version             ");
    fprintf(stderr, "Print version number and exit\n");
    fprintf(stderr, "\n");
}



void sky_model_parse_cmdline( int argc, char **argv, sky_opts *opts )
{
    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"analytic_output", required_argument, 0, 'A'},
                {"conjugate",       no_argument,       0, 'c'},
                {"frequency",       required_argument, 0, 'f'},
                {"fee2016_output",  required_argument, 0, 'F'},
                {"help",            no_argument,       0, 'h'},
                {"header",          no_argument,       0, 'H'},
                {"metafits-file",   required_argument, 0, 'm'},
                {"pointings_file",  required_argument, 0, 'p'},
                {"no_parallactic",  no_argument,       0, 'P'},
                {"no_swap_columns", no_argument,       0, 's'},
                {"version",         no_argument,       0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "A:cf:F:hHm:p:PsV",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {
                case 'A':
                    opts->analytic_output = strdup(optarg);
                    break;
                case 'c':
                    opts->conjugate = 1;
                    break;
                case 'f':
                    opts->frequency = atoi(optarg);
                    break;
                case 'F':
                    opts->fee2016_output = strdup(optarg);
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'H':
                    opts->print_header = 1;
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'p':
                    opts->pointings = strdup(optarg);
                    break;
                case 'P':
                    opts->parallactic = 0;
                    break;
                case 's':
                    opts->swap_columns = 0;
                    break;
                case 'V':
                    fprintf( stderr, "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(0);
                    break;
                default:
                    fprintf(stderr, "error: sky_model_parse_cmdline: "
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

    // Set other defaults
    if (opts->analytic_output == NULL)
        opts->analytic_output = strdup( "sky_model_comparison_ANALYTIC.dat" );
    if (opts->fee2016_output == NULL)
        opts->fee2016_output = strdup( "sky_model_comparison_FEE2016.dat" );
}


void fprint_header( int argc, char **argv, FILE *f, unsigned int *delays,
        double *amps, sky_opts *opts, struct metafits_info *mi,
        int zenith_norm )
{
    int i;
    fprintf( f, "# VCSTOOLS version %s\n", VERSION_BEAMFORMER );
    fprintf( f, "#\n" );
    fprintf( f, "# Generating command:\n" );
    fprintf( f, "#" );
    for (i = 0; i < argc; i++)
        fprintf( f, " %s", argv[i] );
    fprintf( f, "\n" );
    fprintf( f, "#\n" );
    fprintf( f, "# Frequency: %ld Hz\n", opts->frequency );
    fprintf( f, "# Delays:" );
    for (i = 0; i < NDELAYS; i++)
        fprintf( f, "%5u", delays[i] );
    fprintf( f, "\n" );
    fprintf( f, "# Amps:  " );
    for (i = 0; i < NDELAYS; i++)
        fprintf( f, "%5.1lf", amps[i] );
    fprintf( f, "\n" );
    fprintf( f, "# Tile pointing (az, za) (deg): %lf, %lf\n",
            mi->tile_pointing_az,
            90.0 - mi->tile_pointing_el );
    fprintf( f, "# FEE2016 beam zenith normalised? %s\n",
            (zenith_norm ? "yes" : "no") );
    fprintf( f, "#\n" );
    fprintf( f, "# Columns:\n" );
    fprintf( f, "# [ (1)+(2)i,  (3)+(4)i ]\n" );
    fprintf( f, "# [ (5)+(6)i,  (7)+(8)i ]\n" );
    fprintf( f, "#\n" );
}

void fprint_jones( FILE *f, ComplexDouble *J )
{
    int i;
    for (i = 0; i < 4; i++)
        fprintf( f, "%lf %lf ", CReald( J[i] ), CImagd( J[i] ) );
}
