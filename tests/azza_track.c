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
#include "slalib.h"
#include "slamac.h"
#include "../make_beam/beam_common.h"
#include "psrfits.h"
#include "../make_beam/mycomplex.h"
#include "../make_beam/form_beam.h"

typedef struct azza_opts_t
{
    char   *metafits;     // Filename of the metafits file
    char   *output;       // Name of the output file
    int     print_header; // boolean: print header in output files?
    int     skip_seconds; // Interval between outputs (in seconds)
    double  ra;           // Right ascension in radians
    double  dec;          // Declination in radians
} azza_opts;

void fprint_header( int, char **, FILE * );
void azza_parse_cmdline( int, char **, azza_opts * );

int main(int argc, char **argv)
{
    // Create structure for command line options and set defaults
    azza_opts opts;

    opts.metafits     = NULL; // Filename of the metafits file
    opts.output       = NULL; // Name of the output file
    opts.print_header = 0;    // boolean: print header in output files?
    opts.skip_seconds = 1;    // Interval between outputs (in seconds)
    opts.ra           = NAN;  // Right ascension in radians
    opts.dec          = NAN;  // Declination in radians

    // Parse command line arguments
    azza_parse_cmdline( argc, argv, &opts );

    // Get a file handle for writing
    FILE *fout = NULL;
    if (opts.output == NULL)
        fout = stdout;
    else
        fout = fopen( opts.output, "w" );

    // Read in info from metafits file
    struct metafits_info mi;
    int chan_width = 10000; // This isn't used, but is needed for the following call
    get_metafits_info( opts.metafits, &mi, chan_width );

    // Convert the date obs to mjd
    double intmjd, fracmjd, date_obs_mjd, mjd;
    utc2mjd( mi.date_obs, &intmjd, &fracmjd );
    date_obs_mjd = intmjd + fracmjd;

    // Prepare other needed variables for calculation
    double lst; // LST
    double ha;  // Hour angle
    double eq = 2000.0;
    double ra_ap, dec_ap;
    double app_ha_rad;
    double az, za, el;

    // Print header, if requested
    if (opts.print_header)
    {
        fprint_header( argc, argv, fout );
    }

    int s; // seconds offset from beginning of observation
    for (s = 0; s < mi.exposure; s += opts.skip_seconds)
    {
        // Get the MJD and LST of this second
        mjd = date_obs_mjd + (s + 0.5)/86400.0;
        mjd2lst( mjd, &lst );

        // Convert to hour angle
        slaMap(opts.ra, opts.dec, 0, 0, 0, 0, eq, mjd, &ra_ap, &dec_ap);
        ha = slaRanorm( lst - ra_ap )*DR2H;
        app_ha_rad = ha * DH2R;

        slaDe2h(app_ha_rad, dec_ap, MWA_LAT*DD2R, &az, &el);
        za = DPIBY2 - el;

        fprintf( fout, "%lf %lf\n", az, za );
    }

    // Free up memory and close files
    free( opts.metafits );
    if (opts.output != NULL)
    {
        fclose( fout );
        free( opts.output );
    }


    return 0;
}


void usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: az_za_track [OPTIONS]\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "DESCRIPTION\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\tOutputs the azimuth (az) and zenith angle (za) "
                      "pointings for a given observation.\n" );

    fprintf(stderr, "\n");
    fprintf(stderr, "REQUIRED OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-D, --dec=DD:MM:SS.S      ");
    fprintf(stderr, "The declination\n");
    fprintf(stderr, "\t-m, --metafits-file=FILE  ");
    fprintf(stderr, "FILE is the metafits file\n");
    fprintf(stderr, "\t-R, --ra=HH:MM:SS.S       ");
    fprintf(stderr, "The right ascension\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OTHER OPTIONS\n");
    fprintf(stderr, "\t-h, --help                ");
    fprintf(stderr, "Print this help and exit\n");
    fprintf(stderr, "\t-H, --header              ");
    fprintf(stderr, "Include header in output files (default: no header)\n");
    fprintf(stderr, "\t-o, --output-file=FILE    ");
    fprintf(stderr, "FILE is the metafits file (if not supplied, output "
                    "will be written to stdout).\n");
    fprintf(stderr, "\t-s, --skip-seconds        ");
    fprintf(stderr, "Interval between pointings in seconds (default: 1)\n");
    fprintf(stderr, "\t-V, --version             ");
    fprintf(stderr, "Print version number and exit\n");
    fprintf(stderr, "\n");
}



void azza_parse_cmdline( int argc, char **argv, azza_opts *opts )
{
    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"dec",             required_argument, 0, 'D'},
                {"help",            required_argument, 0, 'h'},
                {"header",          no_argument,       0, 'H'},
                {"metafits-file",   required_argument, 0, 'm'},
                {"output-file",     required_argument, 0, 'o'},
                {"ra",              required_argument, 0, 'R'},
                {"skip-seconds",    required_argument, 0, 's'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv,
                             "D:hHm:o:R:s:V",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'D':
                    opts->dec = parse_dec( optarg )*DD2R;
                    break;
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                    break;
                case 'H':
                    opts->print_header = 1;
                    break;
                case 'm':
                    opts->metafits = strdup(optarg);
                    break;
                case 'o':
                    opts->output = strdup(optarg);
                    break;
                case 'R':
                    opts->ra = parse_ra( optarg )*DH2R;
                    break;
                case 's':
                    opts->skip_seconds = atoi( optarg );
                    break;
                case 'V':
                    fprintf( stderr, "MWA Beamformer %s\n", VERSION_BEAMFORMER);
                    exit(EXIT_SUCCESS);
                    break;
                default:
                    fprintf(stderr, "error: azza_parse_cmdline: "
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
    if ( opts->metafits == NULL )
    {
        fprintf( stderr, "error: required option -m missing\n" );
        exit(EXIT_FAILURE);
    }
    if (isnan( opts->ra ) || isnan( opts->dec ))
    {
        fprintf( stderr, "error: required options -R, -D missing\n" );
        exit(EXIT_FAILURE);
    }

    // Check that something sensible was given for skip_seconds
    if (opts->skip_seconds <= 0)
    {
        fprintf( stderr, "error: -s option must be >1\n" );
        exit(EXIT_FAILURE);
    }
}

void fprint_header( int argc, char **argv, FILE *f )
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
    fprintf( f, "# Columns:\n" );
    fprintf( f, "#   az_rad  za_rad\n" );
    fprintf( f, "#\n" );
}

