#ifndef BEAM_COMMON_H
#define BEAM_COMMON_H

#include <inttypes.h>
#include "mycomplex.h"

// Calibration solution types
#define NO_CALIBRATION  0
#define RTS             1
#define RTS_BANDPASS    2
#define OFFRINGA        3

#define MAX_POLS   4

#define BUFSIZE    4096
#define VEL_LIGHT  299792458.0
#define N_COPOL    2
#define R2C_SIGN   -1.0

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

struct make_beam_opts {
    // Variables for required options
    char              *obsid;         // The observation ID
    unsigned long int  begin;         // GPS time -- when to start beamforming
    unsigned long int  end;           // GPS time -- when to stop beamforming
    char              *time_utc;      // utc time string "yyyy-mm-ddThh:mm:ss"
    char              *dec_ddmmss;    // "dd:mm:ss"
    char              *ra_hhmmss;     // "hh:mm:ss"
    char              *datadir;       // The path to where the recombined data live
    char              *metafits;      // filename of the metafits file
    char              *rec_channel;   // 0 - 255 receiver 1.28MHz channel
    long int           frequency;     // = rec_channel expressed in Hz

    // Variables for MWA/VCS configuration
    int                nstation;      // The number of antennas
    int                nchan;         // The number of fine channels (per coarse channel)
    unsigned int       chan_width;    // The bandwidth of an individual fine chanel (Hz)
    unsigned int       sample_rate;   // The VCS sample rate (Hz)
    int                use_ant_flags; // Use flags in metafits file?

    // Output options
    int                out_incoh;     // Default = PSRFITS (incoherent) output turned OFF
    int                out_coh;       // Default = PSRFITS (coherent)   output turned OFF
    int                out_vdif;      // Default = VDIF                 output turned OFF
    int                out_uvdif;     // Default = upsampled VDIF       output turned OFF

    struct calibration cal;           // Variables for calibration settings
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
        ComplexDouble       ***complex_weights_array,  // output
        ComplexDouble      ****invJi                   // output
);


void get_metafits_info( char *metafits, struct metafits_info *mi, unsigned int chan_width );
void destroy_metafits_info( struct metafits_info *mi );

void int8_to_uint8(int n, int shift, char * to_convert);
void float2int8_trunc(float *f, int n, float min, float max, int8_t *i);

void flatten_bandpass(
        int nstep,
        int nchan,
        int npol,
        void *data);

void read_data( char *filename, uint8_t *data, int nbytes );
int read_rts_file(ComplexDouble **G, ComplexDouble *Jref,
                  double *amp, char *fname);
int read_bandpass_file( ComplexDouble ***Jm, ComplexDouble ***Jf,
                        int chan_width, int nchan, int nant, char *filename );
int read_offringa_gains_file( ComplexDouble **antenna_gain, int nant,
                              int coarse_chan, char *gains_file, int *order );


void dec2hms( char *out, double in, int sflag );

/**** MATRIX OPERATIONS ****/

void cp2x2(ComplexDouble *Min, ComplexDouble *Mout);
void inv2x2(ComplexDouble *Min, ComplexDouble *Mout);
void inv2x2S(ComplexDouble *Min, ComplexDouble **Mout);
void mult2x2d(ComplexDouble *M1, ComplexDouble *M2, ComplexDouble *Mout);
void conj2x2(ComplexDouble *M, ComplexDouble *Mout);
double norm2x2(ComplexDouble *M, ComplexDouble *Mout);

#endif
