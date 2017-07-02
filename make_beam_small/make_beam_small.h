#ifndef MAKE_BEAM_H
#define MAKE_BEAM_H

#include <complex.h>
#include "psrfits.h"

/* A structure to read in all the relevant info from the observation metafits
 * file.
 */
/*
typedef struct metafits_info_t {
    float      *N_array;
    float      *E_array;
    float      *H_array;
    float      *cable_array;
    int        *flag_array;
    short int  *antenna_num;
    char      **tilenames;
    int         ninput;
} metafits_info;
*/

struct delays {
    double mean_ra;
    double mean_dec;
    double az;
    double el;
    double lmst;
    double fracmjd;
    double intmjd;
};

/* Running get_delays from within make_beam */
void get_delays(
        int coarse_chan,
        char *dec_ddmmss,
        char *ra_hhmmss,
        long int frequency,
        char *metafits,
        int get_offringa,
        int get_rts,
        char *DI_Jones_file,
        float samples_per_sec,
        long int chan_width,
        char *time_utc,
        double sec_offset,
        struct delays *delay_vals,
        complex double **complex_weights_array,  // output
        double *weights_array,
        complex double **invJi           // output
);


void printf_psrfits( struct psrfits *pf );  /* Prints values in psrfits struct to stdout */

#endif
