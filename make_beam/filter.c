/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "filter.h"
#include "mycomplex.h"


ComplexDouble *roots_of_unity( int N )
/* Creates a complex-valued array containing the N roots of unity.
 * The caller must free this memory (by passing the returned pointer to
 * free()).
 */
{
    // Allocate memory
    ComplexDouble *roots = (ComplexDouble *)malloc( N*sizeof(ComplexDouble) );

    // Make sure memory was allocated correctly
    if (!roots)
    {
        fprintf( stderr, "error: roots_of_unity: could not allocate "
                         "memory\n" );
        exit(EXIT_FAILURE);
    }

    // Calculate the roots and store them in the roots array
    int n;
    for (n = 0; n < N; n++)
    {
        // e^{2πin/N}
        double th = (double)n/(double)N;
        roots[n] = CExpd( CMaked( 0.0, 2*M_PI*th ) );
    }

    return roots;
}


