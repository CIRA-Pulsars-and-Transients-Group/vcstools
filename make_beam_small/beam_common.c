#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mycomplex.h>
#include "beam_common.h"
#include "psrfits.h"

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
            int di = 0; // Index for data[]
            for (i = 0;i<nstep;i++) {
                for (p = 0;p<npol;p++) {
                    for (j = 0;j<nchan;j++){

                        band[p][j] += CAbsf(data[di++]);
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


void cp2x2(ComplexDouble *Min, ComplexDouble *Mout)
{
    Mout[0] = Min[0];
    Mout[1] = Min[1];
    Mout[2] = Min[2];
    Mout[3] = Min[3];
}


void inv2x2(ComplexDouble *Min, ComplexDouble *Mout)
{
    ComplexDouble m1 = CMuld( Min[0], Min[3] );
    ComplexDouble m2 = CMuld( Min[1], Min[2] );
    ComplexDouble det = CSubd( m1, m2 );
    ComplexDouble inv_det = CRcpd( det );
    Mout[0] = CMuld(       inv_det,  Min[3] );
    Mout[1] = CMuld( CNegd(inv_det), Min[1] );
    Mout[2] = CMuld( CNegd(inv_det), Min[2] );
    Mout[3] = CMuld(       inv_det,  Min[0] );
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
    Mout[0] = CSumd( m00, m12 );
    Mout[1] = CSumd( m01, m13 );
    Mout[2] = CSumd( m20, m32 );
    Mout[3] = CSumd( m21, m33 );
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
        Fnorm += CMuld( M[i], CConjd(M[i]) );

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

