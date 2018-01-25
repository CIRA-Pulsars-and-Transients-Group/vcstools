#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "filter.h"

void create_filter( filter *fil, int size, int ntaps )
{
    fil->size  = size;
    fil->ntaps = ntaps;
    fil->coeffs = (complex double *)malloc( size * sizeof(complex double) );
}


void destroy_filter( filter *fil )
{
    free( fil->coeffs );
}

void load_filter( char *filename, int dtype, int ntaps, filter *fil )
/* This function allocates memory for, and loads into memory, the filter
 * coefficients contained in [filename]. The file is expected to be white-
 * space separated ASCII text.
 *
 * dtype can be REAL_COEFFS (=0) or CPLX_COEFFS (=1). If dtype is REAL_COEFFS,
 * the file will be assumed to formatted as one float per entry. If dtype is
 * CPLX_COEFFS, two floats will be read in per filter coefficient and will be
 * assumed to represent the real and imaginary parts respectively. In either
 * case, the output is cast to complex doubles.
 *
 * The memory allocated here should be freed using destroy_filter().
 */
{
    // Check that a valid dtype was passed in
    if (dtype != REAL_COEFFS && dtype != CPLX_COEFFS)
    {
        fprintf( stderr, "error: load_filter: unrecognised dtype (%d). Valid "
                         "values are REAL_COEFFS(=0) or CPLX_COEFFS(=1).\n",
                         dtype );
        exit(EXIT_FAILURE);
    }

    // Variables to read the numbers in [filename] into
    double re, im;

    // Open the file for reading
    FILE *f = fopen( filename, "r" );
    if (!f)
    {
        fprintf( stderr, "error: load_filter: Couldn't open file %s\n",
                 filename );
        exit(EXIT_FAILURE);
    }

    // Run through file and count how many coefficients there are
    int ncoeffs = 0;
    while (fscanf( f, "%lf", &re) != EOF)
        ncoeffs++;

    if (dtype == CPLX_COEFFS)
        ncoeffs /= 2;

    // Allocate memory for that many coefficients
    create_filter( fil, ncoeffs, ntaps );

    // Rewind to the start of the file and read in the coefficients
    rewind(f);
    int i;
    for (i = 0; i < ncoeffs; i++)
    {
        // Read in a real value
        fscanf( f, "%lf", &re );

        // Read in or assign an imag value
        if (dtype == REAL_COEFFS)
            im = 0.0;
        else // dtype == CPLX_COEFFS
            fscanf( f, "%lf", &im );

        fil->coeffs[i] = re + im*I;
    }

    // Close the file
    fclose(f);
}

void upsample( complex double *x, int xsize, int ufact, complex double *y )
/* Upsample signal x by inserting ufact-1 zeros between each element of x.
 * 
 * Inputs:
 *   complex double *x     = the input array
 *   int             xsize = the size of array x
 *   int             ufact = the upsampling factor
 *
 * Outputs:
 *   complex double *y = the output array. Must have size at least xsize*ufact.
 */
{
    int i, j;
    for (i = 0; i < xsize; i++)
    {
        y[i*ufact] = x[i];
        for (j = 1; j < ufact; j++)
            y[i*ufact + j] = 0.0;
    }
}

void apply_phase_ramp( filter *in, double slope, filter *out )
/* Applies a phase ramp across the input filter. If in->size differs from
 * out->size, then only the first [size] elements are calculated, where
 * size is the smaller of in->size and out->size.
 *
 * Inputs:
 *   filter *in   = an arbitrary input array
 *   double slope = the phase ramp slope, in revolutions,
 *                  i.e. unique in interval [0,1), so
 *                  slope = 0 is the same as slope = 1
 *
 * Outputs:
 *   filter *out  = the output filter
 */
{
    // Get the size of the smallest filter
    int size = (in->size < out->size ? in->size : out->size);

    // Calculate the phase ramp and apply it
    int i;
    for (i = 0; i < size; i++)
        out->coeffs[i] = in->coeffs[i] * cexp( 2*M_PI*I*slope*i );
}

void apply_mult_phase_ramps( filter *in, int N, filter outs[] )
/* Applies multiple phase ramps to the array x. The slopes are chosen such
 * that the nth ramp has slope n/N (in revolutions, see apply_phase_ramp()).
 *
 * Inputs:
 *   filter *in     = an arbitrary input array
 *   int     N      = the number of ramps to apply
 *
 * Outputs:
 *   filter *outs[] = the output filters (there must be at least N).
 */
{
    double slope;
    int n;
    for (n = 0; n < N; n++)
    {
        slope = (double)n / (double)N;
        apply_phase_ramp( in, slope, &(outs[n]) );
    }
}

void fir_filter_1D( filter *fil, complex double *signal, int size,
                    complex double *res )
/* This implementation of a simple FIR filter is designed to
 * match scipy's "lfilter" function.
 *
 * Inputs:
 *   filter *fil            = the filter to be applied
 *   complex double *signal = the signal to be filtered
 *   int             size   = the size of the signal
 *
 * Outputs:
 *   complex double *res = the result of the filtering operation. It is
 *                         assumed that res points to a block of memory equal
 *                         to or bigger than signal
 */
{
    int n, i, m;

    for (n = 0; n < size; n++)
    {
        // Reset res element to zero
        res[n] = 0.0;

        // Sum of signal weighted by coefficients
        for (i = 0; i < fil->size; i++)
        {
            m = n - i;

            if (m < 0)
                continue;

            res[n] += signal[m] * fil->coeffs[i];
        }
    }
}

int test_fir_filter_1D()
{
    int       test_success = 1;
    double    tol          = 1e-8;
    int       size         = 20;

    // Create a custom filter
    filter fil;
    create_filter( &fil, 5, 12 );
    fil.coeffs[0] = -1.0+0.5*I;
    fil.coeffs[1] =  3.5+0.6*I;
    fil.coeffs[2] = -6.2-0.7*I;
    fil.coeffs[3] =  4.2-0.8*I;
    fil.coeffs[4] =  0.1+0.1*I;

    // Create a custom signal
    complex double signal[] = { -0.2792024707796282  + 0.14465327403336403*I,
                                -1.0731649697881844  + 1.6093704244328415 *I,
                                 1.7993616760985767  - 0.6111110177007358 *I,
                                -0.8477317486335060  - 1.2089562799954578 *I,
                                -0.4137688089010841  + 2.586665621944545  *I,
                                 0.3602252310877078  - 0.5419415562745188 *I,
                                 1.0689860399053286  - 0.9709662927314314 *I,
                                -0.16163105521320711 - 1.098072153391534  *I,
                                 1.1514768730460412  + 0.397194363766194  *I,
                                -0.6216168936876152  - 2.03873831556825   *I,
                                -1.213465847198062   - 2.64379586598601   *I,
                                 0.8886769156948763  - 0.36005340064800706*I,
                                -0.7902763726457631  - 1.16122584905192   *I,
                                -0.17489273465168498 - 0.2975180708506569 *I,
                                -1.6772004599784247  + 0.19889432497112922*I,
                                 0.517066561838601   + 0.3786269126784602 *I,
                                -1.0109700847494039  - 0.17189634510155302*I,
                                 1.242139345745306   + 0.24583995557764265*I,
                                -0.3428196213690234  - 0.5861638816699087 *I,
                                -0.7468335148086475  + 0.47996872640690363*I };

    // Create the known solution
    complex double known[] = { 0.20687583 - 0.28425451*I,
                              -0.79552085 - 1.80718793*I,
                              -4.38319321 + 5.79828079*I,
                              14.83989672 - 8.67015658*I,
                             -17.96725357 + 2.60023192*I,
                               8.12039947 +13.66352865*I,
                               1.09189180 -20.2034956 *I,
                               2.78956848 +12.35748807*I,
                              -7.78498003 - 0.83703135*I,
                               9.46948554 + 5.77876646*I,
                              -6.63148584 -13.21292639*I,
                               4.30500723 + 4.51940071*I,
                               6.20376833 + 9.36952709*I,
                             -14.57729539 -13.11733283*I,
                               8.81842831 + 2.96015829*I,
                              -9.94311433 - 2.65531193*I,
                              12.28187496 - 0.06207099*I,
                             -14.61395460 - 1.41254754*I,
                              13.27054399 + 4.82266309*I,
                             -12.24024464 - 5.32797534*I };

    // Create a results array
    complex double res[size];

    // Run the function to be tested
    fir_filter_1D( &fil, signal, size, res );

    // Compare the result with the known solution
    int i;
    for (i = 0; i < size; i++)
    {
        //fprintf( stderr, "(%15e, %15e)  (%15e, %15e)  ",
        //                     creal(res[i]),   cimag(res[i]),
        //                     creal(known[i]), cimag(known[i]) );
        if ((abs(creal(res[i]) - creal(known[i])) > tol) ||
            (abs(cimag(res[i]) - cimag(known[i])) > tol))
        {
            test_success = 0;
        //    fprintf( stderr, "test_fir_filter_1D() failed: diff > %e", tol );
            break;
        }
        //fprintf( stderr, "\n" );
    }

    // Free memory
    destroy_filter( &fil );

    return test_success;
}

int test_upsample()
{
    int test_success = 1;

    int ufact = 3;
    int xsize = 5;
    complex double x[] = { 1.0, 2.0*I, 3.0, 4.0*I, -5.0 };
    complex double y[xsize*ufact];
    upsample( x, xsize, ufact, y );

    complex double ans[] = { 1.0,   0.0, 0.0,
                             2.0*I, 0.0, 0.0,
                             3.0  , 0.0, 0.0,
                             4.0*I, 0.0, 0.0,
                            -5.0  , 0.0, 0.0 };

    int i;
    for (i = 0; i < xsize*ufact; i++)
        if (y[i] != ans[i])
            test_success = 0;

    return test_success;
}

void run_all_tests()
{
    int successful;

    // Test the test_fir_filter_1D() function
    successful = test_fir_filter_1D();
    if (successful)
        fprintf( stderr, "FIR filter test successful\n" );
    else
        fprintf( stderr, "FIR filter test failed\n" );

    // Test the upsample() function
    successful = test_upsample();
    if (successful)
        fprintf( stderr, "Upsample test successful\n" );
    else
        fprintf( stderr, "Upsample test failed\n" );
}

/*
void main()
{
    run_all_tests();

    filter fil;
    int ntaps = 12;
    load_filter( "filter_coeffs.txt", REAL_COEFFS, ntaps, &fil );

    int N = fil.size;
    filter fil_ramps[N];
    int n;
    for (n = 0; n < N; n++)
        create_filter( &fil_ramps[n], fil.size, ntaps );

    apply_mult_phase_ramps( &fil, N, fil_ramps );

    int i;
    for (n = 0; n < N; n++)
    {
        for (i = 0; i < fil.size; i++)
            printf( "%lf ", carg(fil_ramps[n].coeffs[i]), cimag(fil.coeffs[i]) );

        printf( "\n" );
    }

    destroy_filter( &fil );
    for (n = 0; n < N; n++)
        destroy_filter( &fil_ramps[n] );
}
*/
