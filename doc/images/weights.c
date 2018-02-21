#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "filter.h"

int main()
{
    // Load filter coefficients
    filter fil;
    load_filter( "filter_coeffs.txt", REAL_COEFFS, 12, &fil );
    int i;
    for (i = 0; i < fil.size; i++)
        fil.coeffs[i] /= 120000.0;

    // Set up bank of filters with phase ramps applied
    int nchan = 128;
    filter fil_ramps[nchan];
    int ch;
    for (ch = 0; ch < nchan; ch++)
        create_filter( &fil_ramps[ch], fil.size, fil.ntaps );

    apply_mult_phase_ramps( &fil, nchan, fil_ramps );

    // Print them out, row by row
    FILE *f_abs = fopen( "weights_abs.txt", "w" );
    FILE *f_arg = fopen( "weights_arg.txt", "w" );
    for (ch = 0; ch < nchan; ch++)
    {
        for (i = 0; i < fil.size; i++)
        {
            fprintf( f_abs, "%f ", cabs(fil_ramps[ch].coeffs[i]) );
            fprintf( f_arg, "%f ", carg(fil_ramps[ch].coeffs[i]) );
        }

        fprintf( f_abs, "\n" );
        fprintf( f_arg, "\n" );
    }

    // Clean up
    fclose( f_abs );
    fclose( f_arg );

    for (ch = 0; ch < nchan; ch++)
        destroy_filter( &fil_ramps[ch] );

    return EXIT_SUCCESS;
}
