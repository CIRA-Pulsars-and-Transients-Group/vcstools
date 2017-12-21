#ifndef BEAM_PSRFITS_H
#define BEAM_PSRFITS_H

#include <complex.h>
#include "psrfits.h"

void printf_psrfits( struct psrfits *pf );  /* Prints values in psrfits struct to stdout */

void populate_psrfits_header(
        struct psrfits *pf,
        char           *metafits,
        char           *obsid,
        char           *time_utc,
        unsigned int    sample_rate,
        long int        frequency,
        int             nchan,
        long int        chan_width,
        int             outpol,
        char           *rec_channel,
        struct delays  *delay_vals );

#endif
