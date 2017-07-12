#ifndef MAKE_BEAM_H
#define MAKE_BEAM_H

#include <complex.h>
#include "psrfits.h"

// Calibration solution types
#define RTS           0
#define RTS_BANDPASS  1
#define OFFRINGA      2

// A structure to read in all the relevant info from the observation metafits
// file.
struct metafits_info {
    double      tile_pointing_ra;
    double      tile_pointing_dec;
    double      tile_pointing_az;
    double      tile_pointing_el;
    float      *N_array;
    float      *E_array;
    float      *H_array;
    float      *cable_array;
    int        *flag_array;
    double     *weights_array;
    short int  *antenna_num;
    char      **tilenames;
    int         ninput;
    int         chan_width;
} metafits_info;

struct delays {
    double mean_ra;
    double mean_dec;
    double az;
    double el;
    double lmst;
    double fracmjd;
    double intmjd;
};

struct calibration {
    char *filename;        // The file that houses the calibration solution
    int   cal_type;        // Either RTS or OFFRINGA
    int   offr_chan_num;   // The channel number in the Offringa calibration solution file
};

/* Running get_delays from within make_beam */
void get_delays(
        char                  *dec_ddmmss,
        char                  *ra_hhmmss,
        long int               frequency,
        struct                 calibration *cal,
        float                  samples_per_sec,
        char                  *time_utc,
        double                 sec_offset,
        struct delays         *delay_vals,
        struct metafits_info  *mi,
        complex double       **complex_weights_array,  // output
        complex double       **invJi                   // output
);


void printf_psrfits( struct psrfits *pf );  /* Prints values in psrfits struct to stdout */

#endif
