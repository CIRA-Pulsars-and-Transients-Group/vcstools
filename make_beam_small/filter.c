#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

void apply_phase_ramp( complex double *x, int size, double slope,
                       complex double *y )
/* Applies a phase ramp across the array x.
 *
 * Inputs:
 *   complex double *x     = an arbitrary input array
 *   int             size  = the size of x
 *   double          slope = the phase ramp slope, in revolutions,
 *                           i.e. unique in interval [0,1), so
 *                           slope = 0 is the same as slope = 1
 *
 * Outputs:
 *   complex double *y = the output array. It is assumed to point to a
 *                       block of memory at least as big as x.
 */
{
    int i;
    for (i = 0; i < size; i++)
        y[i] = x[i] * cexp( 2*M_PI*I*slope*i );
}

void apply_mult_phase_ramps( complex double *x, int xsize, int N,
                             complex double **y )
/* Applies multiple phase ramps to the array x. The slopes are chosen such
 * that the nth ramp has slope n/N (in revolutions, see apply_phase_ramp()).
 *
 * Inputs:
 *   complex double *x     = an arbitrary input array
 *   int             xsize = the size of x
 *   int             N     = the number of ramps to apply
 *
 * Outputs:
 *   complex double **y = the output arrays. It is assumed that y points to
 *                        N valid arrays, each of size xsize (or greater).
 */
{
    double slope;
    int n;
    for (n = 0; n < N; n++)
    {
        slope = (double)n / (double)N;
        apply_phase_ramp( x, xsize, slope, y[n] );
    }
}

void fir_filter_1D( complex double *coeff, complex double *signal, int *size,
                    complex double *res )
/* This implementation of a simple FIR filter is designed to
 * match scipy's "lfilter" function.
 *
 * Inputs:
 *   complex double *coeff  = the filter coefficients
 *   complex double *signal = the signal to be filtered
 *   int            *size   = an array of 2 ints:
 *                            { size_of_coeff, size_of_signal }
 *
 * Outputs:
 *   complex double *res = the result of the filtering operation. It is
 *                         assumed that res points to a block of memory equal
 *                         to or bigger than signal
 */
{
    int n, i, m;

    for (n = 0; n < size[1]; n++)
    {
        // Reset res element to zero
        res[n] = 0.0;

        // Sum of signal weighted by coefficients
        for (i = 0; i < size[0]; i++)
        {
            m = n - i;

            if (m < 0)
                continue;

            res[n] += signal[m] * coeff[i];
        }
    }
}

int test_fir_filter_1D()
{
    int       test_success = 1;
    double    tol          = 1e-8;
    int       size[2]      = { 5, 20 };

    // Create a custom filter
    complex double coeff[] = { -1.0+0.5*I,
                                3.5+0.6*I,
                               -6.2-0.7*I,
                                4.2-0.8*I,
                                0.1+0.1*I };

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
    complex double res[size[1]];

    // Run the function to be tested
    fir_filter_1D( coeff, signal, size, res );

    // Compare the result with the known solution
    int i;
    for (i = 0; i < size[1]; i++)
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

    return test_success;
}

void run_all_tests()
{
    // Test the test_fir_filter_1D() function
    int successful = test_fir_filter_1D();
    if (successful)
        fprintf( stderr, "FIR filter test successful\n" );
    else
        fprintf( stderr, "FIR filter test failed\n" );
}

void main()
{
    //run_all_tests();

    int xsize = 128;
    int N     = 128;
    complex double x[xsize];
    complex double *y[N];
    int i;
    for (i = 0; i < N; i++)
        y[i] = (complex double *)malloc( xsize * sizeof(complex double) );

    for (i = 0; i < xsize; i++)
        x[i] = 1.0 + 0.0*I;

    apply_mult_phase_ramps(x, xsize, N, (complex double **)y);

    for (i = 0; i < xsize; i++)
        printf( "%d %f %f\n", i, creal(y[2][i]), cimag(y[2][i]) );

    for (i = 0; i < N; i++)
        free(y[i]);
}
