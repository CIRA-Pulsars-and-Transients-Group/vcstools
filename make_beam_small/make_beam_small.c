#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <complex.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include "mwac_utils.h"
#include "slalib.h"
#include "slamac.h"
#include "ascii_header.h"
#include "mwa_header.h"
#include <omp.h>
//#include <mpi.h>
#include <glob.h>
#include <fcntl.h>
#include <assert.h>
#include "beam_common.h"
#include "beam_psrfits.h"
#include "make_beam_small.h"

// Are GPU available

#ifdef HAVE_CUDA
#include "gpu_utils.h"
#include <cuda_runtime.h>
#else
#define Complex float _Complex
#endif

//
// write out psrfits directly
#include "psrfits.h"
#include "antenna_mapping.h"
#include "beamer_version.h"

#define MAX_COMMAND_LENGTH 1024

//#define PROFILE

void usage();



void get_metafits_info( char *metafits, struct metafits_info *mi, unsigned int chan_width ) {
/* Read in the relevant information from the metafits file.
 * This function allocates dynamic memory. Destroy it with
 *   destroy_metafits_info(...)
 */

    fitsfile *fptr = NULL;
    int status     = 0;
    int anynull    = 0;

    // Open the metafits file
    fits_open_file(&fptr, metafits, READONLY, &status);
    if (fptr == NULL) {
        fprintf( stderr, "Failed to open metafits file \"%s\"\n", metafits );
        exit(EXIT_FAILURE);
    }

    // Read in the tile pointing information (and channel width)
    fits_read_key(fptr, TDOUBLE, "RA",       &(mi->tile_pointing_ra),  NULL, &status);
    fits_read_key(fptr, TDOUBLE, "DEC",      &(mi->tile_pointing_dec), NULL, &status);
    fits_read_key(fptr, TDOUBLE, "AZIMUTH",  &(mi->tile_pointing_az),  NULL, &status);
    fits_read_key(fptr, TDOUBLE, "ALTITUDE", &(mi->tile_pointing_el),  NULL, &status);
    mi->chan_width = chan_width;

    if (status != 0) {
        fprintf(stderr, "Fits status set: failed to read az/alt, ");
        fprintf(stderr, "ra/dec fits keys from the header\n");
        exit(EXIT_FAILURE);
    }

    // Move to the binary table
    fits_movnam_hdu(fptr, BINARY_TBL, "TILEDATA", 0, &status);
    if (status != 0) {
        fprintf(stderr, "Error: Failed to move to TILEDATA HDU\n");
        exit(EXIT_FAILURE);
    }

    // Read in the number of inputs (= nstation * npol)
    fits_read_key(fptr, TINT, "NAXIS2", &(mi->ninput), NULL, &status);
    if (status != 0) {
        fprintf(stderr, "Error: Failed to read size of binary table in TILEDATA\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory
    mi->N_array         =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->E_array         =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->H_array         =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->cable_array     =     (float *)malloc( mi->ninput*sizeof(float)    );
    mi->flag_array      =       (int *)malloc( mi->ninput*sizeof(int)      );
    mi->weights_array   =    (double *)malloc( mi->ninput*sizeof(double)   );
    mi->antenna_num     = (short int *)malloc( mi->ninput*sizeof(short int));
    mi->tilenames       =     (char **)malloc( mi->ninput*sizeof(char *)   );
    int i;
    for (i = 0; i < (int)(mi->ninput); i++) {
        mi->tilenames[i] = (char *)malloc( 32*sizeof(char) );
    }
    char *testval = (char *) malloc(1024);

    /* Read the columns */
    int colnum;

    // Cable lengths
    for (i=0; i < (int)(mi->ninput); i++) {

        fits_get_colnum(fptr, 1, "Length", &colnum, &status);
        if(fits_read_col_str(fptr, colnum, i+1, 1, 1, "0.0", &testval, &anynull, &status)) {
            fprintf(stderr, "Error: Failed to cable column  in metafile\n");
            exit(EXIT_FAILURE);
        }

        sscanf(testval, "EL_%f", &(mi->cable_array[i]));
    }

    // Tile names
    fits_get_colnum(fptr, 1, "TileName", &colnum, &status);
    if (status != 0) {
        status = 0;
        fits_get_colnum(fptr, 1, "Tile", &colnum, &status);
    }
    if (status != 0) {
        fprintf(stderr, "Could not find either column \"TileName\" or \"Tile\" in metafits file\n");
        exit(EXIT_FAILURE);
    }

    fits_read_col(fptr, TSTRING, colnum, 1, 1, mi->ninput, NULL, mi->tilenames, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read Tile(Name) in metafile\n");
        exit(EXIT_FAILURE);
    }

    // North coordinate
    fits_get_colnum(fptr, 1, "North", &colnum, &status);
    fits_read_col_flt(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->N_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read  N coord in metafile\n");
        exit(EXIT_FAILURE);
    }

    // East coordinate
    fits_get_colnum(fptr, 1, "East", &colnum, &status);
    fits_read_col_flt(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->E_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read E coord in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Height coordinate
    fits_get_colnum(fptr, 1, "Height", &colnum, &status);
    fits_read_col_flt(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->H_array, &anynull, &status);

    if (status != 0){
        fprintf(stderr, "Error: Failed to read H coord in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Antenna number
    fits_get_colnum(fptr, 1, "Antenna", &colnum, &status);
    fits_read_col_sht(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->antenna_num, &anynull, &status);

    if (status != 0){
        fprintf(stderr, "Error: Failed to read field \"Antenna\" in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Flags & weights
    float wgt_sum = 0.0;
    fits_get_colnum(fptr, 1, "Flag", &colnum, &status);
    fits_read_col_int(fptr, colnum, 1, 1, mi->ninput, 0.0, mi->flag_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read flags column in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Invert value (flag off = full weight; flag on = zero weight)
    for (i = 0; i < mi->ninput; i++) {
        mi->weights_array[i] = 1.0 - (double)mi->flag_array[i];
        wgt_sum += mi->weights_array[i];
        // This differs from Ord's orig code, which sums squares. However,
        // all values should be = 1, so end result should be the same
    }

    // Exit with error if there are no weights
    if (wgt_sum == 0.0) {
        fprintf(stderr, "Zero weight sum on read\n");
        exit(EXIT_FAILURE);
    }

    // Clean up
    free( testval );
    fits_close_file(fptr, &status);
}


void destroy_metafits_info( struct metafits_info *mi ) {
/* Frees the memory allocated in the metafits_info struct
 */
    free( mi->N_array       );
    free( mi->E_array       );
    free( mi->H_array       );
    free( mi->cable_array   );
    free( mi->flag_array    );
    free( mi->weights_array );
    free( mi->antenna_num   );
    int i;
    for (i = 0; i < mi->ninput; i++)
        free( mi->tilenames[i] );
    free( mi->tilenames     );
}


void int8_to_uint8(int n, int shift, char * to_convert) {
    int j;
    int scratch;
    int8_t with_sign;

    for (j = 0; j < n; j++) {
        with_sign = (int8_t) *to_convert;
        scratch = with_sign + shift;
        *to_convert = (uint8_t) scratch;
        to_convert++;
    }
}
void float2int8_trunc(float *f, int n, float min, float max, int8_t *i) /*includefile*/
{
    int j;
    for (j = 0; j < n; j++) {
        f[j] = (f[j] > max) ? (max) : f[j];
        f[j] = (f[j] < min) ? (min) : f[j];
        i[j] = (int8_t) rint(f[j]);

    }
}
void to_offset_binary(int8_t *i, int n){
    int j;
    for (j = 0; j < n; j++) {
        i[j] = i[j] ^ 0x80;
    }
}

void flatten_bandpass(int nstep, int nchan, int npol, void *data, float *scales, float *offsets, int new_var, int iscomplex, int normalise, int update, int clear, int shutdown) {
    // putpose is to generate a mean value for each channel/polaridation

    int i=0, j=0;
    int p=0;
    float *data_ptr = NULL;

    static float **band;

    static float **chan_min;

    static float **chan_max;


    static int setup = 0;

    if (setup == 0) {
        band = (float **) calloc (npol, sizeof(float *));
        chan_min = (float **) calloc (npol, sizeof(float *));
        chan_max = (float **) calloc (npol, sizeof(float *));
        for (i=0;i<npol;i++) {
            band[i] = (float *) calloc(nchan, sizeof(float));
            chan_min[i] = (float *) calloc(nchan, sizeof(float));
            chan_max[i] = (float *) calloc(nchan, sizeof(float));
        }
        setup = 1;
    }

    if (update) {
        for (p = 0;p<npol;p++) {
            for (j=0;j<nchan;j++){

                band[p][j] = 0.0;
            }
        }

        if (iscomplex == 0) {
            data_ptr = (float *) data;

            for (i=0;i<nstep;i++) {
                for (p = 0;p<npol;p++) {
                    for (j=0;j<nchan;j++){


                        if (i==0) {
                            chan_min[p][j] = *data_ptr;
                            chan_max[p][j] = *data_ptr;
                        }
                        band[p][j] += fabsf(*data_ptr);
                        if (*data_ptr < chan_min[p][j]) {
                            chan_min[p][j] = *data_ptr;
                        }
                        else if (*data_ptr > chan_max[p][j]) {
                            chan_max[p][j] = *data_ptr;
                        }
                        data_ptr++;
                    }
                }

            }
        }
        else {
            complex float  *data_ptr = (complex float *) data;
            for (i=0;i<nstep;i++) {
                for (p = 0;p<npol;p++) {
                    for (j=0;j<nchan;j++){

                        band[p][j] += cabsf(*data_ptr);
                        data_ptr++;
                    }
                }

            }

        }

    }
    // set the offsets and scales - even if we are not updating ....

    float *out=scales;
    float *off = offsets;
    for (p = 0;p<npol;p++) {
        for (j=0;j<nchan;j++){

            // current mean
            *out = ((band[p][j]/nstep))/new_var; // removed a divide by 32 here ....
            //fprintf(stderr, "Channel %d pol %d mean: %f normaliser %f (max-min) %f\n", j, p, (band[p][j]/nstep), *out, (chan_max[p][j] - chan_min[p][j]));
            out++;
            *off = 0.0;

            off++;

        }
    }
    // apply them to the data

    if (normalise) {

        data_ptr = (float *) data;

        for (i=0;i<nstep;i++) {
            float *normaliser = scales;
            float *off  = offsets;
            for (p = 0;p<npol;p++) {
                for (j=0;j<nchan;j++){

                    *data_ptr = ((*data_ptr) - (*off))/(*normaliser); // 0 mean normalised to 1
                    //fprintf(stderr, "%f %f %f\n", *data_ptr, *off, *normaliser);
                    off++;
                    data_ptr++;
                    normaliser++;
                }
            }

        }
    }

    // clear the weights if required

    if (clear) {

        float *out=scales;
        float *off = offsets;
        for (p = 0;p<npol;p++) {
            for (j=0;j<nchan;j++){

                // reset
                *out = 1.0;



                out++;
                *off = 0.0;

                off++;

            }
        }
    }

    // free the memory
    if (shutdown) {
        for (i=0;i<npol;i++) {
            free(band[i]);
            free(chan_min[i]);
            free(chan_max[i]);
        }


        free(band);
        free(chan_min);
        free(chan_max);
        setup = 0;
    }
}


void read_data( char *filename, uint8_t *data, int nbytes ) {

    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (f == NULL) {
        fprintf( stderr, "Error opening file '%s'\n", filename );
        exit(EXIT_FAILURE);
    }

    // Read in nbytes bytes
    int nread = fread( (void *)data, sizeof(uint8_t), nbytes, f );

    // Check that I got all nbytes
    // Any other number is disallowed
    if (nread != nbytes) {
        fprintf( stderr, "Error: file %s does not contain %d bytes\n",
                filename, nbytes );
        exit(EXIT_FAILURE);
    }

    fclose(f);

}

void correct_stt( struct psrfits *pf ) {
    /* now we have to correct the STT_SMJD/STT_OFFS as they will have been broken by the write_psrfits*/
    int    itmp    = 0;
    int    itmp2   = 0;
    double dtmp    = 0;
    int    status  = 0;

    //fits_open_file(&(pf.fptr), pf.filename, READWRITE, &status);

    fits_read_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp,  NULL, &status);
    fits_read_key(pf->fptr, TINT,    "STT_SMJD", &itmp,  NULL, &status);
    fits_read_key(pf->fptr, TINT,    "STT_IMJD", &itmp2, NULL, &status);

    if (dtmp > 0.5) {
        itmp = itmp+1;
        if (itmp == 86400) {
            itmp = 0;
            itmp2++;
        }
    }
    dtmp = 0.0;

    fits_update_key(pf->fptr, TINT, "STT_SMJD", &itmp, NULL, &status);
    fits_update_key(pf->fptr, TINT, "STT_IMJD", &itmp2, NULL, &status);
    fits_update_key(pf->fptr, TDOUBLE, "STT_OFFS", &dtmp, NULL, &status);

}

/*****************
 * MAIN FUNCTION *
 *****************/

int main(int argc, char **argv) {

    // Variables for required options
    char              *obsid       = NULL; // The observation ID
    unsigned long int  begin       = 0;    // GPS time -- when to start beamforming
    unsigned long int  end         = 0;    // GPS time -- when to stop beamforming
    char              *time_utc    = NULL; // utc time string "yyyy-mm-ddThh:mm:ss"
    char              *dec_ddmmss  = NULL; // "dd:mm:ss"
    char              *ra_hhmmss   = NULL; // "hh:mm:ss"
    char              *datadir     = NULL; // The path to where the recombined data live
    char              *metafits    = NULL; // filename of the metafits file
    char              *rec_channel = NULL; // 0 - 255 receiver 1.28MHz channel
    long int           frequency   = 0;    // = rec_channel expressed in Hz

    // Variables for MWA/VCS configuration
    int                nstation      = 128;    // The number of antennas
    int                nchan         = 128;    // The number of fine channels (per coarse channel)
    unsigned int       chan_width    = 10000;  // The bandwidth of an individual fine chanel (Hz)
    unsigned int       sample_rate   = 10000;  // The VCS sample rate (Hz)
    int                use_ant_flags = 0;      // Use flags in metafits file?
    const int          npol          = 2;      // X,Y
    const int          outpol        = 4;      // I,Q,U,V

    // Output options
    int                out_incoh     = 0;  // Default = incoherent output turned OFF
    int                out_coh       = 1;  // Default = coherent   output turned ON

    // Variables for calibration settings
    struct calibration cal;
    cal.filename          = NULL;
    cal.bandpass_filename = NULL;
    cal.chan_width        = 40000;
    cal.nchan             = 0;
    cal.cal_type          = NO_CALIBRATION;
    cal.offr_chan_num     = 0;

    // These are used to calculate how the input data are ordered
    const int npfb = 4;
    const int nrec = 16;
    const int ninc = 4;

    if (argc > 1) {

        int c;
        while (1) {

            static struct option long_options[] = {
                {"obsid",           required_argument, 0, 'o'},
                {"begin",           required_argument, 0, 'b'},
                {"end",             required_argument, 0, 'e'},
                {"incoh",           no_argument,       0, 'i'},
                {"utc-time",        required_argument, 0, 'z'},
                {"dec",             required_argument, 0, 'D'},
                {"ra",              required_argument, 0, 'R'},
                {"data-location",   required_argument, 0, 'd'},
                {"metafits-file",   required_argument, 0, 'm'},
                {"coarse-chan",     required_argument, 0, 'f'},
                {"antennas",        required_argument, 0, 'a'},
                {"num-fine-chans",  required_argument, 0, 'n'},
                {"fine-chan-width", required_argument, 0, 'w'},
                {"sample-rate",     required_argument, 0, 'r'},
                {"use-ant-flags",   no_argument,       0, 'F'},
                {"dijones-file",    required_argument, 0, 'J'},
                {"bandpass-file",   required_argument, 0, 'B'},
                {"rts-chan-width",  required_argument, 0, 'W'},
                {"offringa-file",   required_argument, 0, 'O'},
                {"offringa-chan",   required_argument, 0, 'C'},
                {"help",            required_argument, 0, 'h'},
                {"version",         required_argument, 0, 'V'}
            };

            int option_index = 0;
            c = getopt_long( argc, argv, "a:b:B:C:d:D:e:f:FhiJ:m:n:o:O:r:R:Vw:W:z:",
                             long_options, &option_index);
            if (c == -1)
                break;

            switch(c) {

                case 'a':
                    nstation = atoi(optarg);
                    break;
                case 'b':
                    begin = atol(optarg);
                    break;
                case 'B':
                    cal.bandpass_filename = strdup(optarg);
                    cal.cal_type = RTS_BANDPASS;
                    break;
                case 'C':
                    cal.offr_chan_num = atoi(optarg);
                    break;
                case 'd':
                    datadir = strdup(optarg);
                    break;
                case 'D':
                    dec_ddmmss = strdup(optarg);
                    break;
                case 'e':
                    end = atol(optarg);
                    break;
                case 'f':
                    rec_channel = strdup(optarg);
                    frequency = atoi(optarg) * 1.28e6 - 640e3; // The base frequency of the coarse channel in Hz
                    break;
                case 'F':
                    use_ant_flags = 1;
                    break;
                case 'h':
                    usage();
                    exit(0);
                    break;
                case 'i':
                    out_incoh = 1;
                    break;
                case 'J':
                    cal.filename = strdup(optarg);
                    if (cal.cal_type != RTS_BANDPASS)
                        cal.cal_type = RTS;
                    break;
                case 'm':
                    metafits = strdup(optarg);
                    break;
                case 'n':
                    nchan = atoi(optarg);
                    break;
                case 'o':
                    obsid = strdup(optarg);
                    break;
                case 'O':
                    cal.filename = strdup(optarg);
                    cal.cal_type = OFFRINGA;
                    break;
                case 'r':
                    sample_rate = atoi(optarg);
                    break;
                case 'R':
                    ra_hhmmss = strdup(optarg);
                    break;
                case 'V':
                    printf("%s\n", MAKE_BEAM_VERSION);
                    exit(0);
                    break;
                case 'w':
                    chan_width = atoi(optarg);
                    break;
                case 'W':
                    cal.chan_width = atoi(optarg);
                    break;
                case 'z':
                    time_utc = strdup(optarg);
                    break;
                default:
                    fprintf(stderr, "Error: unrecognised option '%s'\n", optarg);
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
    if (obsid == NULL || begin == 0 || end  == 0 || time_utc == NULL ||
        dec_ddmmss == NULL || ra_hhmmss == NULL || datadir == NULL ||
        metafits == NULL || rec_channel == NULL)
    {
        fprintf(stderr, "Error: missing required options\n");
        usage();
        exit(EXIT_FAILURE);
    }

    // Check that a calibration solution was supplied
    if (cal.cal_type == NO_CALIBRATION)
    {
        fprintf(stderr, "Error: no calibration solution supplied\n");
        usage();
        exit(EXIT_FAILURE);
    }

    // Check that at least one output option has been selected
    if (!out_incoh && !out_coh)
    {
        fprintf(stderr, "Error: no output format selected\n");
        usage();
        exit(EXIT_FAILURE);
    }

    // Start counting time from here (i.e. after parsing the command line)
    double begintime = omp_get_wtime();
    printf("[%f]  Starting make_beam\n", omp_get_wtime()-begintime);

    // Calculate the number of files
    int nfiles = end - begin + 1;
    if (nfiles <= 0) {
        fprintf(stderr, "Cannot beamform on %d files (between %lu and %lu)\n", nfiles, begin, end);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the file name list
    char **filenames = NULL;
    filenames = (char **)malloc( nfiles*sizeof(char *) );

    // Allocate memory and write filenames
    int second;
    unsigned long int timestamp;
    for (second = 0; second < nfiles; second++) {
        timestamp = second + begin;
        filenames[second] = (char *)malloc( MAX_COMMAND_LENGTH*sizeof(char) );
        sprintf( filenames[second], "%s/%s_%ld_ch%s.dat", datadir, obsid, timestamp, rec_channel );
    }

    // Allocate memory for complex weights matrices
    int ant, p, ch; // Loop variables

    complex double  ***complex_weights_array = NULL; // [ant][ch][pol]
    
    complex_weights_array = (complex double ***)malloc( nstation * sizeof(complex double **) );
    for (ant = 0; ant < nstation; ant++) {
        complex_weights_array[ant] = (complex double **)malloc( nchan * sizeof(complex double *) );
        for (ch = 0; ch < nchan; ch++) {
            complex_weights_array[ant][ch] = (complex double *)malloc( npol * sizeof(complex double) );
        }
    }

    // Allocate memory for (inverse) Jones matrices
    complex double ****invJi = NULL; // [ant][ch][pol][pol]

    invJi = (complex double ****)malloc( nstation * sizeof(complex double ***) );
    for (ant = 0; ant < nstation; ant++) {
        invJi[ant] =(complex double ***)malloc( nchan * sizeof(complex double **) );
        for (ch = 0; ch < nchan; ch++) {
            invJi[ant][ch] = (complex double **)malloc( npol * sizeof(complex double *) );
            for (p = 0; p < npol; p++) {
                invJi[ant][ch][p] = (complex double *)malloc( npol * sizeof(complex double) );
            }
        }
    }

    // Structures for holding psrfits header information
    struct psrfits pf;
    struct psrfits pf_incoh;

    // Read in info from metafits file
    printf("[%f]  Reading in metafits file information from %s\n", omp_get_wtime()-begintime, metafits);
    struct metafits_info mi;
    get_metafits_info( metafits, &mi, chan_width );

    // If using bandpass calibration solutions, calculate number of expected bandpass channels
    if (cal.cal_type == RTS_BANDPASS)
        cal.nchan = (nchan * chan_width) / cal.chan_width;

    int i;
    if (!use_ant_flags)
        for (i = 0; i < nstation*npol; i++)
            mi.weights_array[i] = 1.0;

    double wgt_sum = 0;
    for (i = 0; i < nstation*npol; i++)
        wgt_sum += mi.weights_array[i];
    double invw = 1.0/wgt_sum;

    // Run get_delays to populate the delay_vals struct
    printf("[%f]  Setting up output header information\n", omp_get_wtime()-begintime);
    struct delays delay_vals;
    get_delays(
            dec_ddmmss,    // dec as a string "dd:mm:ss"
            ra_hhmmss,     // ra  as a string "hh:mm:ss"
            frequency,     // middle of the first frequency channel in Hz
            &cal,          // struct holding info about calibration
            sample_rate,   // = 10000 samples per sec
            time_utc,      // utc time string
            0.0,           // seconds offset from time_utc at which to calculate delays
            &delay_vals,   // Populate psrfits header info
            &mi,           // Struct containing info from metafits file
            NULL,          // complex weights array (ignore this time)
            NULL           // invJi array           (ignore this time)
    );

    // now we need to create a fits file and populate its header
    populate_psrfits_header( &pf, metafits, obsid, time_utc, sample_rate,
            frequency, nchan, chan_width, outpol, rec_channel, &delay_vals );

    // If incoherent sum requested, populate incoherent psrfits header
    if (out_incoh)
    {
        int incoh_npol = 1;
        populate_psrfits_header( &pf_incoh, metafits, obsid, time_utc, sample_rate,
                frequency, nchan, chan_width, incoh_npol, rec_channel, &delay_vals );

        // The above sets the RA and Dec to be the beamforming pointing, but it ought
        // to have the tile pointing, not the beam pointing
        pf_incoh.hdr.ra2000  = mi.tile_pointing_ra;
        pf_incoh.hdr.dec2000 = mi.tile_pointing_dec;
        pf_incoh.hdr.summed_polns = 1;

        // Also, we need to give it a different base name for the output files
        sprintf(pf_incoh.basefilename, "%s_%s_ch%03d_incoh",
                pf_incoh.hdr.project_id, pf_incoh.hdr.source, atoi(rec_channel));
    }

    // Create array for holding the raw data
    int bytes_per_file = sample_rate * nstation * npol * nchan;
    uint8_t *data = (uint8_t *)malloc( bytes_per_file * sizeof(uint8_t) );
    assert(data);

    int8_t *out_buffer_8_psrfits = (int8_t *)malloc( outpol*nchan*pf.hdr.nsblk * sizeof(int8_t) );
    float  *data_buffer_psrfits  =  (float *)malloc( nchan*outpol*pf.hdr.nsblk * sizeof(float) );
    int8_t *out_buffer_8_incoh   = NULL;
    float  *data_buffer_incoh    = NULL;
    if (out_incoh)
    {
        out_buffer_8_incoh = (int8_t *)malloc( nchan*pf_incoh.hdr.nsblk * sizeof(int8_t) );
        data_buffer_incoh = (float *)malloc( nchan*pf_incoh.hdr.nsblk * sizeof(float) );
    }

    int offset_in_psrfits;
    int offset_in_incoh;

    int file_no = 0;
    int sample;

    printf("[%f]  **BEGINNING BEAMFORMING**\n", omp_get_wtime()-begintime);
    for (file_no = 0; file_no < nfiles; file_no++) {

        // Read in data from next file
        printf("[%f]  Reading in data from %s [%d/%d]\n", omp_get_wtime()-begintime,
                filenames[file_no], file_no+1, nfiles);
        read_data( filenames[file_no], data, bytes_per_file  );

        // Get the next second's worth of phases / jones matrices, if needed
        printf("[%f]  Calculating delays\n", omp_get_wtime()-begintime);
        get_delays(
                dec_ddmmss,    // dec as a string "dd:mm:ss"
                ra_hhmmss,     // ra  as a string "hh:mm:ss"
                frequency,     // middle of the first frequency channel in Hz
                &cal,          // struct holding info about calibration
                sample_rate,   // = 10000 samples per sec
                time_utc,      // utc time string
                (double)file_no, // seconds offset from time_utc at which to calculate delays
                NULL,          // Don't update delay_vals
                &mi,           // Struct containing info from metafits file
                complex_weights_array,  // complex weights array (answer will be output here)
                invJi );       // invJi array           (answer will be output here)

#ifdef PROFILE
        printf("[%f]  Calculating beam --       ", omp_get_wtime()-begintime);
#else
        printf("[%f]  Calculating beam\n", omp_get_wtime()-begintime);
#endif

        offset_in_psrfits  = 0;
        offset_in_incoh    = 0;

        for (i = 0; i < nchan*outpol*pf.hdr.nsblk; i++)
            data_buffer_psrfits[i] = 0.0;

        if (out_incoh)
            for (i = 0; i < nchan*pf_incoh.hdr.nsblk; i++)
                data_buffer_incoh[i] = 0.0;

#ifdef PROFILE
        double sec1=0.0, sec2=0.0, sec3=0.0, sec4=0.0, sec5=0.0;
        int progress = 0;
        int num_threads = omp_get_max_threads();
#endif
#ifdef PROFILE
#pragma omp parallel for reduction(+:sec1,sec2,sec3,sec4,sec5,progress)
#else
#pragma omp parallel for
#endif
        for (sample = 0; sample < (int)sample_rate; sample++ ) {

#ifdef PROFILE
            double chkpt1, chkpt2, chkpt3, chkpt4, chkpt5, chkpt6, chkpt7;
            printf("\b\b\b\b\b\b%5.1f%%", progress*num_threads/100.0 + 0.05);
            fflush(stdout);
            progress++;
#endif

            int ch, ant, pol, opol, opol1, opol2;
            uint8_t uD, uDr, uDi;
            int sDr, sDi;
            int pfb, rec, inc;
            int data_idx;

            complex float beam[nchan][nstation][npol];
            float         incoh_beam[nchan][nstation][npol];
            complex float detected_beam[nchan][npol];
            float         detected_incoh_beam[nchan];
            float         spectrum[nchan*outpol];
            complex float noise_floor[nchan][npol][npol];
            complex float e_true[npol], e_dash[npol];


            // Initialise noise floor to zero
            for (ch    = 0; ch    < nchan;  ch++   )
            for (opol1 = 0; opol1 < npol; opol1++)
            for (opol2 = 0; opol2 < npol; opol2++)
                noise_floor[ch][opol1][opol2] = 0.0;

            // Initialise detected beam to zero
            for (ch  = 0; ch  < nchan   ; ch++ )
            for (pol = 0; pol < npol    ; pol++)
                detected_beam[ch][pol] = 0.0 + 0.0*I;

            // Initialise incoherent beam arrays to zero, if necessary
            if (out_incoh)
                for (ch  = 0; ch  < nchan   ; ch++ )
                    detected_incoh_beam[ch] = 0.0;

            // Calculate the beam, noise floor
            for (ant = 0; ant < nstation; ant++) {

                // Get the index for the data that corresponds to this
                //   sample, channel, antenna, polarisation
                // Justification for the rather bizarre mapping is found in
                // the docs.
                // (rec depends on polarisation, so is calculating in the inner loop)
                pfb = ant / 32;
                inc = (ant / 8) % 4;

                for (ch = 0; ch < nchan; ch++ ) {

#ifdef PROFILE
                    chkpt1 = omp_get_wtime();
#endif
                    // Calculate quantities that depend only on "input" polarisation
                    for (pol = 0; pol < npol    ; pol++) {

                        rec = (2*ant+pol) % 16;

                        data_idx = sample * (ninc*nrec*npfb*nchan) +
                                   ch     * (ninc*nrec*npfb)       +
                                   pfb    * (ninc*nrec)            +
                                   rec    * (ninc)                 +
                                   inc;

                        uD  = data[data_idx];
                        uDr = uD & 0xf;        // Real part = least significant nibble
                        uDi = (uD >> 4) & 0xf; // Imag part = most  significant nibble

                        // Convert from unsigned to signed
                        sDr = (uDr >= 0x8 ? (signed int)uDr - 0x10 : (signed int) uDr);
                        sDi = (uDi >= 0x8 ? (signed int)uDi - 0x10 : (signed int) uDi);

                        // Form a single complex number
                        e_dash[pol]  = (float)sDr + (float)sDi * I;

                        // Detect the incoherent beam, if requested
                        if (out_incoh)
                            incoh_beam[ch][ant][pol] = creal(e_dash[pol] * conj(e_dash[pol]));

                        // Apply complex weights
                        e_dash[pol] *= complex_weights_array[ant][ch][pol];

                    }

#ifdef PROFILE
                    chkpt2 = omp_get_wtime();
                    sec1 += chkpt2-chkpt1;
#endif
                    // Calculate quantities that depend on output polarisation
                    // (i.e. apply inv(jones))
                    for (pol = 0; pol < npol; pol++) {

                        e_true[pol] = 0.0 + 0.0*I;

                        for (opol = 0; opol < npol; opol++)
                            e_true[pol] += invJi[ant][ch][pol][opol] * e_dash[opol];

                        for (opol = 0; opol < npol; opol++)
                            noise_floor[ch][pol][opol] += e_true[pol] * conj(e_true[opol]);

                        beam[ch][ant][pol] = e_true[pol];
                    }
#ifdef PROFILE
                    chkpt3 = omp_get_wtime();
                    sec2 += chkpt3-chkpt2;
#endif
                }
            }

#ifdef PROFILE
            chkpt4 = omp_get_wtime();
#endif
            // Detect the beam = sum over antennas
            for (ant = 0; ant < nstation; ant++)
            for (pol = 0; pol < npol    ; pol++)
            for (ch  = 0; ch  < nchan   ; ch++ )
            {
                detected_beam[ch][pol] += beam[ch][ant][pol];

                // ...and the incoherent beam, if requested
                if (out_incoh)
                    detected_incoh_beam[ch] += incoh_beam[ch][ant][pol];
            }

#ifdef PROFILE
            chkpt5 = omp_get_wtime();
            sec3 += chkpt5-chkpt4;
#endif

            // Calculate the Stokes parameters
            double beam00, beam11;
            double noise0, noise1, noise3;
            complex double beam01;
            unsigned int stokesIidx, stokesQidx, stokesUidx, stokesVidx;

            for (ch = 0; ch < nchan; ch++) {

                beam00 = (double)(detected_beam[ch][0] * conj(detected_beam[ch][0]));
                beam11 = (double)(detected_beam[ch][1] * conj(detected_beam[ch][1]));
                beam01 = detected_beam[ch][0] * conj(detected_beam[ch][1]);

                noise0 = noise_floor[ch][0][0];
                noise1 = noise_floor[ch][0][1];
                noise3 = noise_floor[ch][1][1];

                stokesIidx = 0*nchan + ch;
                stokesQidx = 1*nchan + ch;
                stokesUidx = 2*nchan + ch;
                stokesVidx = 3*nchan + ch;

                // Looking at the dspsr loader the expected order is <ntime><npol><nchan>
                // so for a single timestep we do not have to interleave - I could just stack these
                spectrum[stokesIidx] = (beam00 + beam11 - noise0 - noise3) * invw;
                spectrum[stokesQidx] = (beam00 - beam11 - noise0 + noise3) * invw;
                spectrum[stokesUidx] = 2.0 * (creal(beam01) - noise1)*invw;
                spectrum[stokesVidx] = -2.0 * cimag((beam01 - noise1)*invw);
            }

#ifdef PROFILE
            chkpt6 = omp_get_wtime();
            sec4 += chkpt6-chkpt5;
#endif

            offset_in_psrfits  = sizeof(float)*nchan*outpol * sample;
            offset_in_incoh    = sizeof(float)*nchan * sample;

            memcpy((void *)((char *)data_buffer_psrfits + offset_in_psrfits), spectrum, sizeof(float)*nchan*outpol);

            if (out_incoh)
                memcpy((void *)((char *)data_buffer_incoh + offset_in_incoh), detected_incoh_beam, sizeof(float)*nchan);

#ifdef PROFILE
            chkpt7 = omp_get_wtime();
            sec5 += chkpt7-chkpt6;
#endif

        }
#ifdef PROFILE
        printf("\n");
        double sectotal = sec1 + sec2 + sec3 + sec4 + sec5;
        printf("%.1f  %.1f  %.1f  %.1f  %.1f out of %.1f sec on %d threads\n",
                100.0*sec1/sectotal,
                100.0*sec2/sectotal,
                100.0*sec3/sectotal,
                100.0*sec4/sectotal,
                100.0*sec5/sectotal,
                sectotal/num_threads,
                num_threads);
#endif

        // We've arrived at the end of a second's worth of data...

        printf("[%f]  Flattening bandpass\n", omp_get_wtime()-begintime);
        flatten_bandpass(pf.hdr.nsblk, nchan, outpol,
                data_buffer_psrfits, pf.sub.dat_scales,
                pf.sub.dat_offsets, 32, 0, 1, 1, 1, 0);

        float2int8_trunc(data_buffer_psrfits, pf.hdr.nsblk*nchan*outpol,
                -126.0, 127.0, out_buffer_8_psrfits);

        int8_to_uint8(pf.hdr.nsblk*nchan*outpol, 128,
                (char *) out_buffer_8_psrfits);

        memcpy(pf.sub.data, out_buffer_8_psrfits, pf.sub.bytes_per_subint);

        if (out_incoh)
        {
            flatten_bandpass(pf_incoh.hdr.nsblk, nchan, 1,
                    data_buffer_incoh, pf_incoh.sub.dat_scales,
                    pf_incoh.sub.dat_offsets, 32, 0, 1, 1, 1, 0);

            float2int8_trunc(data_buffer_incoh, pf_incoh.hdr.nsblk*nchan,
                    -126.0, 127.0, out_buffer_8_incoh);

            int8_to_uint8(pf_incoh.hdr.nsblk*nchan, 128,
                    (char *) out_buffer_8_incoh);

            memcpy(pf_incoh.sub.data, out_buffer_8_incoh, pf_incoh.sub.bytes_per_subint);
        }

        printf("[%f]  Writing data to file\n", omp_get_wtime()-begintime);
        if (psrfits_write_subint(&pf) != 0) {
            fprintf(stderr, "Write subint failed file exists?\n");
            break; // Exit from while loop
        }

        pf.sub.offs = roundf(pf.tot_rows * pf.sub.tsubint) + 0.5*pf.sub.tsubint;
        pf.sub.lst += pf.sub.tsubint;

        if (out_incoh)
        {
            if (psrfits_write_subint(&pf_incoh) != 0) {
                fprintf(stderr, "Write incoherent subint failed file exists?\n");
                break; // Exit from while loop
            }

            pf_incoh.sub.offs = roundf(pf_incoh.tot_rows * pf_incoh.sub.tsubint) +
                                0.5*pf_incoh.sub.tsubint;
            pf_incoh.sub.lst += pf_incoh.sub.tsubint;
        }

    }

    printf("[%f]  **FINISHED BEAMFORMING**\n", omp_get_wtime()-begintime);
    printf("[%f]  Starting clean-up\n", omp_get_wtime()-begintime);

    if (pf.status == 0) {
        correct_stt( &pf );

        flatten_bandpass(pf.hdr.nsblk, nchan, outpol, data_buffer_psrfits,
                pf.sub.dat_scales, pf.sub.dat_offsets, 32, 0, 0, 0, 0, 1);
    }

    if (pf_incoh.status == 0) {
        correct_stt( &pf_incoh );

        flatten_bandpass(pf_incoh.hdr.nsblk, nchan, 1, data_buffer_incoh,
                pf_incoh.sub.dat_scales, pf_incoh.sub.dat_offsets, 32, 0, 0, 0, 0, 1);
    }

    // Free up memory for filenames
    if (datadir) {
        int second;
        for (second = 0; second < nfiles; second++)
            free( filenames[second] );
        free( filenames );
    }

    destroy_metafits_info( &mi );
    free( out_buffer_8_psrfits );
    free( out_buffer_8_incoh );
    free( data_buffer_psrfits  );
    free( data_buffer_incoh  );
    free( data );

    return 0;
}


void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "usage: make_beam_small [OPTIONS]\n");

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
    fprintf(stderr, "\t-D, --dec=dd:mm:ss.s      ");
    fprintf(stderr, "Declination of pointing direction\n");
    fprintf(stderr, "\t-R, --ra=hh:mm:ss.s       ");
    fprintf(stderr, "Right ascension of pointing direction\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-d, --data-location=PATH  ");
    fprintf(stderr, "PATH is the directory containing the recombined data\n");
    fprintf(stderr, "\t-m, --metafits-file=FILE  ");
    fprintf(stderr, "FILE is the metafits file pertaining to the OBSID given by the\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr,  "-o option\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-f, --coarse-chan=N       ");
    fprintf(stderr, "Absolute coarse channel number (0-255)\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "OUTPUT OPTIONS\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\t-i, --incoh               ");
    fprintf(stderr, "Turn on incoherent beam output (off by default)\n");
    fprintf(stderr, "\t--no-psrfits              ");
    fprintf(stderr, "Turn off PSRFITS output (on by default) (NOT YET IMPLEMENTED)\n");

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
    fprintf(stderr, "\t-F, --use-ant-flags       ");
    fprintf(stderr, "Only include those antennas in the beamformer that have not      ");
    fprintf(stderr, "[default: off]\n");
    fprintf(stderr, "\t                          ");
    fprintf(stderr, "been flagged in the metafits file given by the -m option.\n");

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
    fprintf(stderr, "\t-V, --version             ");
    fprintf(stderr, "Print version number and exit\n");
    fprintf(stderr, "\n");
}

