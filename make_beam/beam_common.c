/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "beam_common.h"
#include "psrfits.h"
#include "mycomplex.h"
#include "mwa_hyperbeam.h"

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
    mi->delays          =      (int **)malloc( mi->ninput*sizeof(int*)     );
    mi->amps            =   (double **)malloc( mi->ninput*sizeof(double*)  );
    int i;
    for (i = 0; i < (int)(mi->ninput); i++) {
        mi->tilenames[i] =   (char *)malloc( 32*sizeof(char)        );
        mi->delays[i]    =    (int *)malloc( NDELAYS*sizeof(int)    );
        mi->amps[i]      = (double *)malloc( NDELAYS*sizeof(double) );
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
    fits_read_col_int(fptr, colnum, 1, 1, mi->ninput, 0, mi->flag_array, &anynull, &status);
    if (status != 0){
        fprintf(stderr, "Error: Failed to read flags column in metafile\n");
        exit(EXIT_FAILURE);
    }

    // Beamforming delays (this only reads a single set of delays from the first tile,
    // and assumes that all tiles are the same)
    fits_get_colnum(fptr, 1, "Delays", &colnum, &status);
    for (i=0; i<mi->ninput; i++){
        fits_read_col_int(fptr, colnum, i+1, 1, NDELAYS, 0, mi->delays[i], &anynull, &status);
        if (status != 0){
            fprintf(stderr, "Error: Failed to read delays column in metafile\n");
            exit(EXIT_FAILURE);
        }

        // The amps should all be '1', except when the corresponding delay = '32'
        for (int j = 0; j < NDELAYS; j++){
            //fprintf(stderr, "%d ", mi->delays[i][j]);
            mi->amps[i][j] = (mi->delays[i][j] == 32 ? 0.0 : 1.0);
        }
        //fprintf(stderr, "\n");
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
    for (i = 0; i < mi->ninput; i++){
        free( mi->tilenames[i] );
        free( mi->delays[i]    );
        free( mi->amps[i]      );
    }
    free( mi->tilenames     );
    free( mi->delays        );
    free( mi->amps          );
}


void flatten_bandpass(int nstep, int nchan, int npol, void *data)
{

    // nstep -> ? number of time steps (ie 10,000 per 1 second of data)
    // nchan -> number of (fine) channels (128)
    // npol -> number of polarisations (1 for icoh or 4 for coh)

    // magical mystery normalisation constant
    int new_var = 32;

    // purpose is to generate a mean value for each channel/polaridation

    int i=0, j=0;
    int p=0;

    float *data_ptr = (float *) data;
    float **band;


    band = (float **) calloc (npol, sizeof(float *));
    for (i=0;i<npol;i++) {
      band[i] = (float *) calloc(nchan, sizeof(float));
    }

    // initialise the band array
    for (p = 0;p<npol;p++) {
        for (j=0;j<nchan;j++){
            band[p][j] = 0.0;
        }
    }


    // accumulate abs(data) over all time samples and save into band
    data_ptr = data;
    for (i=0;i<nstep;i++) { // time steps
        for (p = 0;p<npol;p++) { // pols
            for (j=0;j<nchan;j++){ // channels
                band[p][j] += fabsf(*data_ptr);
                data_ptr++;
            }
        }

    }

    // calculate and apply the normalisation to the data
    data_ptr = data;
    for (i=0;i<nstep;i++) {
        for (p = 0;p<npol;p++) {
            for (j=0;j<nchan;j++){
                *data_ptr = (*data_ptr)/( (band[p][j]/nstep)/new_var );
                data_ptr++;
            }
        }

    }

    // free the memory
    for (i=0;i<npol;i++) {
        free(band[i]);
    }
    free(band);
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


int read_rts_file(ComplexDouble **G, ComplexDouble *Jref,
                  double *amp, char *fname)
{
    FILE *fp = NULL;
    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open gain Jones matrix file: %s\n",
                fname);
        exit(EXIT_FAILURE);
    }

    char line[BUFSIZE];
    int index = 0;
    double re0, im0, re1, im1, re2, im2, re3, im3;

    while ((fgets(line, BUFSIZE - 1, fp)) != NULL) {

        if (line[0] == '\n' || line[0] == '#' || line[0] == '\0')
            continue; // skip blank/comment lines
        if (line[0] == '/' && line[1] == '/')
            continue; // also a comment (to match other input files using this style)

        if (index == 0) {

            // read the amplitude and the Alignment Line
            sscanf(line, "%lf", amp);
            fgets(line, BUFSIZE - 1, fp);
            sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);

            Jref[0] = CMaked( re0, im0 );
            Jref[1] = CMaked( re1, im1 );
            Jref[2] = CMaked( re2, im2 );
            Jref[3] = CMaked( re3, im3 );

        }
        if (index > 0) {
            sscanf(line, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf", &re0,
                           &im0, &re1, &im1, &re2, &im2, &re3, &im3);
            G[index - 1][0] = CMaked( re0, im0 );
            G[index - 1][1] = CMaked( re1, im1 );
            G[index - 1][2] = CMaked( re2, im2 );
            G[index - 1][3] = CMaked( re3, im3 );
        }

        index++;

    }

    fclose(fp);

    return 0;

}


int read_bandpass_file(
        ComplexDouble ***Jm, // Output: measured Jones matrices (Jm[ant][ch][pol,pol])
        ComplexDouble ***Jf, // Output: fitted Jones matrices   (Jf[ant][ch][pol,pol])
        int chan_width,       // Input:  channel width of one column in file (in Hz)
        int nchan,            // Input:  (max) number of channels in one file (=128/(chan_width/10000))
        int nant,             // Input:  (max) number of antennas in one file (=128)
        char *filename        // Input:  name of bandpass file
        )
{

    // Open the file for reading
    FILE *f = NULL;
    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error: cannot open Bandpass file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Read in top row = frequency offsets
    int max_len = 4096; // Overkill
    char freqline[max_len];
    if (fgets( freqline, max_len, f ) == NULL) {
        fprintf(stderr, "Error: could not read first line of %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Parse top row
    // Find out which channels are actually in the Bandpass file
    // (i.e. which channels have not been flagged)
    char *freqline_ptr = freqline;
    int pos;
    double freq_offset;
    int chan_count = 0;
    int chan_idxs[nchan];
    int chan_idx;
    while (sscanf(freqline_ptr, "%lf,%n", &freq_offset, &pos) == 1) {

        chan_count++;

        // Make sure we haven't exceeded the total
        if (chan_count > nchan) {
            fprintf(stderr, "Error: More than nchan = %d columns in Bandpass file %s\n", nchan, filename);
            exit(EXIT_FAILURE);
        }
        chan_idx = (int)roundf( freq_offset*1e6 / (double)chan_width );
        chan_idxs[chan_count-1] = chan_idx;

        freqline_ptr += pos;
    }

    // Read in all the values
    int ant, curr_ant = 0;     // The antenna number, read in from the file
    int ant_row       = 0;     // Number between 0 and 7. Each antenna has 8 rows.
    int ch;                    // A counter for channel numbers
    int ci;                    // A counter for channel number indices
    double amp,ph;             // For holding the read-in value pairs (real, imaginary)
    int pol;                   // Number between 0 and 3. Corresponds to position in Jm/Jf matrices: [0,1]
                               //                                                                    [2,3]
    ComplexDouble ***J;       // Either points to Jm or Jf, according to which row we're on

    while (1) {   // Will terminate when EOF is reached

        if (fscanf(f, "%d,", &ant) == EOF)     // Read in first number = antenna number
            break;

        if (ant > nant) {                      // Check that the antenna number is not bigger than expected
            fprintf(stderr, "Error: More than nant = %d antennas in Bandpass file %s\n", nant, filename);
            exit(EXIT_FAILURE);
        }

        ant--;                                 // Convert to 0-offset

        if (ant == curr_ant) {                 // Ensure that there is not an unusual (!=8) number of rows for this antenna
            ant_row++;
            if (ant_row > 8) {
                fprintf(stderr, "Error: More than 8 rows for antenna %d in Bandpass file %s\n",
                        ant, filename);
                exit(EXIT_FAILURE);
            }
        }
        else {
            if (ant_row < 7) {
                fprintf(stderr, "Error: Fewer than 8 rows for antenna %d in Bandpass file %s\n",
                        ant, filename);
                exit(EXIT_FAILURE);
            }
            curr_ant = ant;
            ant_row  = 1;
        }

        if ((ant_row-1) % 2 == 0)  J = Jm;      // Decide if the row corresponds to the Jm values (even rows)
        else                       J = Jf;      // or Jf values (odd rows)

        if (J == NULL) {                        // If the caller doesn't care about this row
            fgets( freqline, max_len, f );      // Skip the rest of this line (freqline isn't needed any more)
            continue;                           // And start afresh on the next line
        }

        pol = (ant_row-1) / 2;                  // Get the polarisation index

        for (ci = 0; ci < chan_count; ci++) {   // Loop over the row

            ch = chan_idxs[ci];                 // Get the channel number
            fscanf(f, "%lf,%lf,", &amp, &ph);   // Read in the re,im pairs in each row

            J[ant][ch][pol] = CScld( CExpd( CMaked(0.0, ph) ), amp );
                                                // Convert to complex number and store in output array
        }

        // (Assumes that the number of values in each row are correct)
    }

    return 1;
}


int read_offringa_gains_file( ComplexDouble **antenna_gain, int nant,
                              int coarse_chan, char *gains_file, int *order )
{
    // Assumes that memory for antenna has already been allocated

    // Open the calibration file for reading
    FILE *fp = NULL;
    fp = fopen(gains_file,"r");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open %s: quitting\n",gains_file);
        exit(EXIT_FAILURE);
    }

    // Read in the necessary information from the header

    uint32_t intervalCount, antennaCount, channelCount, polarizationCount;

    fseek(fp, 16, SEEK_SET);
    fread(&intervalCount,     sizeof(uint32_t), 1, fp);
    fread(&antennaCount,      sizeof(uint32_t), 1, fp);
    fread(&channelCount,      sizeof(uint32_t), 1, fp);
    fread(&polarizationCount, sizeof(uint32_t), 1, fp);

    // Error-checking the info extracted from the header
    if (intervalCount > 1) {
        fprintf(stderr, "Warning: Only the first interval in the calibration ");
        fprintf(stderr, "solution (%s) will be used\n", gains_file);
    }
    if ((int)antennaCount != nant) {
        fprintf(stderr, "Error: Calibration solution (%s) ", gains_file);
        fprintf(stderr, "contains a different number of antennas (%d) ", antennaCount);
        fprintf(stderr, "than specified (%d)\n", nant);
        exit(1);
    }
    if (channelCount != 24) {
        fprintf(stderr, "Warning: Calibration solution (%s) ", gains_file);
        fprintf(stderr, "contains a different number (%d) ", channelCount);
        fprintf(stderr, "than the expected (%d) channels. ", 24);
    }
    if ((int)channelCount <= coarse_chan) {
        fprintf(stderr, "Error: Requested channel number (%d) ", coarse_chan);
        fprintf(stderr, "is more than the number of channels (0-%d) ", channelCount-1);
        fprintf(stderr, "available in the calibration solution (%s)\n", gains_file);
        exit(1);
    }
    int npols = polarizationCount; // This will always = 4

    // Prepare to jump to the first solution to be read in
    int bytes_left_in_header = 16;
    int bytes_to_first_jones = bytes_left_in_header + (npols * coarse_chan * sizeof(ComplexDouble));
         //     (See Offringa's specs for details)
         //     Assumes coarse_chan is zero-offset
         //     sizeof(complex double) *must* be 64-bit x 2 = 16-byte
    int bytes_to_next_jones = npols * (channelCount-1) * sizeof(ComplexDouble);

    int ant, pol;           // Iterate through antennas and polarisations
    int pol_idx, ant_idx;   // Used for "re-ordering" the antennas and pols
    int count = 0;          // Keep track of how many solutions have actually been read in
    double re, im;          // Temporary placeholders for the real and imaginary doubles read in

    // Loop through antennas and read in calibration solution
    int first = 1;
    for (ant = 0; ant < nant; ant++) {

        // Get correct antenna index
        // To wit: The nth antenna in the Offringa binary file will get put into
        // position number order[n]. Default is no re-ordering.
        if (order) {
            ant_idx = order[ant];
        }
        else {
            ant_idx = ant;
        }

        // Jump to next Jones matrix position for this channel
        if (first) {
            fseek(fp, bytes_to_first_jones, SEEK_CUR);
            first = 0;
        }
        else {
            fseek(fp, bytes_to_next_jones, SEEK_CUR);
        }

        // Read in the data
        for (pol = 0; pol < npols; pol++) {

            pol_idx = 3 - pol; // Read them in "backwards", because RTS's "x" = Offringa's "y"

            fread(&re, sizeof(double), 1, fp);
            fread(&im, sizeof(double), 1, fp);

            // Check for NaNs
            if (isnan(re) | isnan(im)) {

                // If NaN, set to identity matrix
                if (pol_idx == 0 || pol_idx == 3)
                    antenna_gain[ant_idx][pol_idx] = CMaked( 1.0, 0.0 );
                else
                    antenna_gain[ant_idx][pol_idx] = CMaked( 0.0, 0.0 );

            }
            else {
                antenna_gain[ant_idx][pol_idx] = CMaked( re, im );
            }

            count++;

        }
    }

    // Close the file, print a summary, and return
    fclose(fp);
    fprintf(stdout, "Read %d inputs from %s\n", count, gains_file);

    return count/npols; // should equal the number of antennas
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


void float2int8_trunc(float *f, int n, float min, float max, int8_t *i)
{
    int j;
    for (j = 0; j < n; j++) {
        f[j] = (f[j] > max) ? (max) : f[j];
        f[j] = (f[j] < min) ? (min) : f[j];
        i[j] = (int8_t) rint(f[j]);

    }
}

void float_to_unit8(float * in, int n, int8_t *out)
{
    int j;
    float min = -128.0; // -126.0 and -128.0 give the same result on test data
    float max = 127.0;
    // use a temp var so we don't modify the input data
    float scratch;
    for (j = 0; j < n; j++) {
        // TODO: count the number of samples that were clipped, store that and put it in the psrfits header
        // the if branching and ternary updates seem to be equivalent execution time
        if (in[j]> max) {
            scratch = max;
        } else if (in[j] < min) {
            scratch = min;
        } else {
            scratch = in[j];
        }
//        scratch = (in[j] > max) ? (max) : in[j];
//        scratch = (in[j] < min) ? (min) : scratch;
        out[j] = (uint8_t)( (int8_t)rint(scratch) + 128);
    }

}

void cp2x2(ComplexDouble *Min, ComplexDouble *Mout)
{
    Mout[0] = Min[0];
    Mout[1] = Min[1];
    Mout[2] = Min[2];
    Mout[3] = Min[3];
}


void inv2x2(ComplexDouble *Min, ComplexDouble *Mout)
{
    ComplexDouble m00 = Min[0];
    ComplexDouble m01 = Min[1];
    ComplexDouble m10 = Min[2];
    ComplexDouble m11 = Min[3];

    ComplexDouble m1 = CMuld( m00, m11 );
    ComplexDouble m2 = CMuld( m01, m10 );

    ComplexDouble det = CSubd( m1, m2 );
    ComplexDouble inv_det = CRcpd( det );

    Mout[0] = CMuld(       inv_det,  m11 );
    Mout[1] = CMuld( CNegd(inv_det), m01 );
    Mout[2] = CMuld( CNegd(inv_det), m10 );
    Mout[3] = CMuld(       inv_det,  m00 );
}

void inv2x2d(double *Min, double *Mout)
{
    double m00 = Min[0];
    double m01 = Min[1];
    double m10 = Min[2];
    double m11 = Min[3];

    double m1 = m00 * m11;
    double m2 = m01 * m10;

    double det = m1 - m2;
    double inv_det = 1/det;

    Mout[0] =  inv_det * m11;
    Mout[1] = -inv_det * m01;
    Mout[2] = -inv_det * m10;
    Mout[3] =  inv_det * m00;
}


void inv2x2S(ComplexDouble *Min, ComplexDouble **Mout)
// Same as inv2x2(), but the output is a 2x2 2D array, instead of a 4-element
// 1D array
{
    ComplexDouble m1 = CMuld( Min[0], Min[3] );
    ComplexDouble m2 = CMuld( Min[1], Min[2] );
    ComplexDouble det = CSubd( m1, m2 );
    ComplexDouble inv_det = CRcpd( det );
    Mout[0][0] = CMuld(       inv_det,  Min[3] );
    Mout[0][1] = CMuld( CNegd(inv_det), Min[1] );
    Mout[1][0] = CMuld( CNegd(inv_det), Min[2] );
    Mout[1][1] = CMuld(       inv_det,  Min[0] );
}


void mult2x2d(ComplexDouble *M1, ComplexDouble *M2, ComplexDouble *Mout)
{
    ComplexDouble m00 = CMuld( M1[0], M2[0] );
    ComplexDouble m12 = CMuld( M1[1], M2[2] );
    ComplexDouble m01 = CMuld( M1[0], M2[1] );
    ComplexDouble m13 = CMuld( M1[1], M2[3] );
    ComplexDouble m20 = CMuld( M1[2], M2[0] );
    ComplexDouble m32 = CMuld( M1[3], M2[2] );
    ComplexDouble m21 = CMuld( M1[2], M2[1] );
    ComplexDouble m33 = CMuld( M1[3], M2[3] );
    Mout[0] = CAddd( m00, m12 );
    Mout[1] = CAddd( m01, m13 );
    Mout[2] = CAddd( m20, m32 );
    Mout[3] = CAddd( m21, m33 );
}

void mult2x2d_RxC(double *M1, ComplexDouble *M2, ComplexDouble *Mout)
/* Mout = M1 x M2
 */
{
    ComplexDouble m00 = CScld( M2[0], M1[0] );
    ComplexDouble m12 = CScld( M2[2], M1[1] );
    ComplexDouble m01 = CScld( M2[1], M1[0] );
    ComplexDouble m13 = CScld( M2[3], M1[1] );
    ComplexDouble m20 = CScld( M2[0], M1[2] );
    ComplexDouble m32 = CScld( M2[2], M1[3] );
    ComplexDouble m21 = CScld( M2[1], M1[2] );
    ComplexDouble m33 = CScld( M2[3], M1[3] );
    Mout[0] = CAddd( m00, m12 );
    Mout[1] = CAddd( m01, m13 );
    Mout[2] = CAddd( m20, m32 );
    Mout[3] = CAddd( m21, m33 );
}

void mult2x2d_CxR(ComplexDouble *M1, double *M2, ComplexDouble *Mout)
/* Mout = M1 x M2
 */
{
    ComplexDouble m00 = CScld( M1[0], M2[0] );
    ComplexDouble m21 = CScld( M1[2], M2[1] );
    ComplexDouble m01 = CScld( M1[0], M2[1] );
    ComplexDouble m13 = CScld( M1[1], M2[3] );
    ComplexDouble m20 = CScld( M1[2], M2[0] );
    ComplexDouble m32 = CScld( M1[3], M2[2] );
    ComplexDouble m12 = CScld( M1[1], M2[2] );
    ComplexDouble m33 = CScld( M1[3], M2[3] );
    Mout[0] = CAddd( m00, m12 );
    Mout[1] = CAddd( m01, m13 );
    Mout[2] = CAddd( m20, m32 );
    Mout[3] = CAddd( m21, m33 );
}

void conj2x2(ComplexDouble *M, ComplexDouble *Mout)
/* Calculate the conjugate of a matrix
 * It is safe for M and Mout to point to the same matrix
 */
{
    int i;
    for (i = 0; i < 4; i++)
        Mout[i] = CConjd(M[i]);
}


double norm2x2(ComplexDouble *M, ComplexDouble *Mout)
/* Normalise a 2x2 matrix via the Frobenius norm
 * It is safe for M and Mout to point to the same matrix.
 */
{
    // Calculate the normalising factor
    double Fnorm = 0.0;
    int i;
    for (i = 0; i < 4; i++)
        Fnorm += CReald( CMuld( M[i], CConjd(M[i]) ) );

    Fnorm = sqrt(Fnorm);

    // Divide each element through by the normalising factor.
    // If norm is 0, then output zeros everywhere
    for (i = 0; i < 4; i++) {
        if (Fnorm == 0.0)
            Mout[i] = CMaked( 0.0, 0.0 );
        else
            Mout[i] = CScld( M[i], 1.0 / Fnorm );
    }

    return Fnorm;
}


void dec2hms( char *out, double in, int sflag )
{
    int sign  = 1;
    char *ptr = out;
    int h, m;
    double s;

    if (in < 0.0)
    {
        sign = -1;
        in = fabs(in);
    }

    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    if (m >= 60)
    {
        // if minutes is 60 convert that to 1 hour
        h += 1;
        m -= 60;
    }
    s = in;
    if (s >= 59.995)
    {
        // if seconds is 60 convert that to 1 minute
        m += 1;
        s = 00.00;
    }
    if (sign==1 && sflag)
    {
        *ptr='+';
        ptr++;
    }
    else if (sign==-1)
    {
        *ptr='-';
        ptr++;
    }
    // Limiting the output's pointings' smallest significant figure to
    // 0.01 arc seconds
    sprintf( ptr, "%2.2d:%2.2d:%05.2f", h, m, s );
}
