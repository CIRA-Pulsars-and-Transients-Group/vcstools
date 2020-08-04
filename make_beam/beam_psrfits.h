/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef BEAM_PSRFITS_H
#define BEAM_PSRFITS_H

#include "psrfits.h"
#include "mycomplex.h"

void printf_psrfits( struct psrfits *pf );  /* Prints values in psrfits struct to stdout */

void populate_psrfits_header(
        struct psrfits  pf[],
        char           *metafits,
        char           *obsid,
        char           *time_utc,
        unsigned int    sample_rate,
        int             max_sec_per_file,
        long int        frequency,
        int             nchan,
        long int        chan_width,
        int             outpol,
        char           *rec_channel,
        struct delays  *delay_vals,
        struct metafits_info mi,
        int             npointing,
        int             is_coherent );

void correct_psrfits_stt( struct psrfits *pf );

void psrfits_write_second( struct psrfits *pf, float *data_buffer, int nchan,
        int outpol, int p);


#endif
