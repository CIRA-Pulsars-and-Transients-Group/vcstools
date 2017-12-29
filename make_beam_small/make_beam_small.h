#ifndef MAKE_BEAM_SMALL_H
#define MAKE_BEAM_SMALL_H

#include "beam_common.h"

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
    int                out_coh;       // Default = PSRFITS (coherent)   output turned ON
    int                out_vdif;      // Default = VDIF                 output turned OFF

    struct calibration cal;           // Variables for calibration settings
};

void usage();
void make_beam_parse_cmdline( int argc, char **argv, struct make_beam_opts *opts );

#endif
