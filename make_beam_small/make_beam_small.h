#ifndef MAKE_BEAM_SMALL_H
#define MAKE_BEAM_SMALL_H

#include <stdlib.h>
#include "beam_common.h"
#include "mycomplex.h"


#define MAX_COMMAND_LENGTH 1024

#define REAL_NIBBLE_TO_UINT8(X)  ((X) & 0xf)
#define IMAG_NIBBLE_TO_UINT8(X)  (((X) >> 4) & 0xf)
#define UINT8_TO_INT(X)          ((X) >= 0x8 ? (signed int)(X) - 0x10 : (signed int)(X))
#define UCMPLX4_TO_CMPLX_FLT(X)  (CMakef((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))), \
                                         (float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X)))))
#define DETECT(X)                (CRealf((X)*CConjf(X)))

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

void usage();
void make_beam_parse_cmdline( int argc, char **argv, struct make_beam_opts *opts );

char **create_filenames( struct make_beam_opts *opts );
void  destroy_filenames( char **filenames, struct make_beam_opts *opts );

ComplexDouble ***create_complex_weights( int nstation, int nchan, int npol );
void             destroy_complex_weights( ComplexDouble ***array, int nstation, int nchan );

ComplexDouble ****create_invJi( int nstation, int nchan, int pol );
void              destroy_invJi( ComplexDouble ****array, int nstation, int nchan, int npol );

ComplexDouble ***create_detected_beam( int nsamples, int nchan, int npol );
void            destroy_detected_beam( ComplexDouble ***array, int nsamples, int nchan );

float *create_data_buffer_psrfits( size_t size );
float *create_data_buffer_vdif( struct vdifinfo *vf );
//float *create_data_buffer_uvdif( size_t );

#endif
