/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "filter.h"
#include "mycomplex.h"


void load_filter( char *filename, int dtype, ComplexDouble *fil )
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
 * The memory allocated here should be freed by the caller.
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
    if (fil) // != NULL
    {
        fprintf( stderr, "load_filter: warning: non-NULL value of filter to "
                         "be initialised\n" );
    }
    fil = (ComplexDouble *)malloc( ncoeffs * sizeof(ComplexDouble) );
    if (!fil) // i.e. fil == NULL
    {
        fprintf( stderr, "load_filter: error: unable to allocate memory for "
                         "filter\n" );
        exit(EXIT_FAILURE);
    }

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

        fil[i] = CMaked( re, im );
    }

    // Close the file
    fclose(f);
}


void apply_phase_ramp( ComplexDouble *in, int size, double slope,
                       ComplexDouble *out )
/* Applies a phase ramp across the input filter. This assumes that both "in"
 * and "out" are at least as big as "size".
 *
 * Inputs:
 *   ComplexDouble  *in   = an arbitrary input array
 *   int             size = the size of arrays in and out
 *   double slope         = the phase ramp slope, in revolutions,
 *                          i.e. unique in interval [0,1), so
 *                          slope = 0 is the same as slope = 1
 *
 * Outputs:
 *   ComplexDouble *out  = the output filter
 */
{
    // Calculate the phase ramp and apply it
    int i;
    for (i = 0; i < size; i++)
        out[i] = CMuld( in[i], CExpd( CMaked( 0.0, 2*M_PI*slope*i ) ) );
}


ComplexDouble **apply_mult_phase_ramps( ComplexDouble *in, int size, int N )
/* Applies multiple phase ramps to the array x. The slopes are chosen such
 * that the nth ramp has slope (n-c)/N, where c=N/2 is the central channel
 * and where the slope is given in revolutions (see apply_phase_ramp()).
 *
 * This function allocates memory for "out", such that it has dimensions
 *   out[N][size]
 * This should be freed by the caller.
 *
 * Inputs:
 *   ComplexDouble *in   = an arbitrary input array
 *   int             size = the size of "in"
 *   int             N    = the number of ramps to apply
 *
 * Outputs:
 *   ComplexDouble **out = the output filters (there must be at least N).
 */
{
    // Allocate (2D) memory for out
    ComplexDouble **out = (ComplexDouble **)malloc( N * sizeof(ComplexDouble *) );
    int n;
    for (n = 0; n < N; n++)
    {
        out[n] = (ComplexDouble *)malloc( size * sizeof(ComplexDouble) );
    }

    // Make sure the memory allocation worked
    if (!out) // i.e. out == NULL
    {
        fprintf( stderr, "apply_mult_phase_ramps: error: unable to allocate "
                         "memory for \"out\" variable\n" );
        exit(EXIT_FAILURE);
    }

    // Apply a phase ramp to each row
    double slope;
    for (n = 0; n < N; n++)
    {
        slope = (double)(n-N/2) / (double)N;
        apply_phase_ramp( in, size, slope, out[n] );
    }
    return out;
}


void fir_filter_1D( ComplexDouble *fil, int fil_size, ComplexDouble *signal,
                    int size, ComplexDouble *res )
/* This implementation of a simple FIR filter is designed to
 * match scipy's "lfilter" function.
 *
 * Inputs:
 *   ComplexDouble  *fil       = the filter to be applied
 *   int             fil_size  = the size of the filter
 *   ComplexDouble  *signal    = the signal to be filtered
 *   int             size      = the size of the signal
 *
 * Outputs:
 *   ComplexDouble  *res       = the result of the filtering operation. It is
 *                               assumed that res points to a block of memory
 *                               equal to or bigger than signal
 */
{
    int n, i, m;

    for (n = 0; n < size; n++)
    {
        // Reset res element to zero
        res[n] = CMaked( 0.0, 0.0 );

        // Sum of signal weighted by coefficients
        for (i = 0; i < fil_size; i++)
        {
            m = n - i;

            if (m < 0)
                continue;

            res[n] = CAddd( res[n], CMuld( signal[m], fil[i] ) );
        }
    }
}



