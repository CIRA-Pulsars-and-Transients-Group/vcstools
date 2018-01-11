#ifndef BEAM_COMMON_H
#define BEAM_COMMON_H

#include <inttypes.h>

// Calibration solution types
#define NO_CALIBRATION  0
#define RTS             1
#define RTS_BANDPASS    2
#define OFFRINGA        3

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
    char *filename;           // The file that houses the calibration solution
    char *bandpass_filename;  // The file that houses the RTS bandpass information
    int   chan_width;         // Channel width used in RTS bandpass solutions (in Hz)
    int   nchan;              // The number of channels in the RTS bandpass solutions
    int   cal_type;           // Either RTS or OFFRINGA
    int   offr_chan_num;      // The channel number in the Offringa calibration solution file
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
        complex double      ***complex_weights_array,  // output
        complex double     ****invJi                   // output
);


void get_metafits_info( char *metafits, struct metafits_info *mi, unsigned int chan_width );
void destroy_metafits_info( struct metafits_info *mi );

void int8_to_uint8(int n, int shift, char * to_convert);
void float2int8_trunc(float *f, int n, float min, float max, int8_t *i);
void flatten_bandpass(
        int nstep,
        int nchan,
        int npol,
        void *data,
        float *scales,
        float *offsets,
        int new_var,
        int iscomplex,
        int normalise,
        int update,
        int clear,
        int shutdown );

void read_data( char *filename, uint8_t *data, int nbytes );

#endif
