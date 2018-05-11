#ifndef BEAM_VDIF_H
#define BEAM_VDIF_H

#include <fftw3.h>
#include "beam_common.h"
#include "vdifio.h"
#include "filter.h"
#include "mycomplex.h"

#define  VDIF_HEADER_SIZE  32

/* convenience type - this just collects all the vdif info together */
struct vdifinfo {

    int frame_length; //length of the vdif frame
    int frame_rate; // frames per second
    size_t samples_per_frame; // number of time samples per vdif frame
    int bits; // bits per sample
    int nchan; // channels per framme
    int chan_width; // channel width in hertz
    int sample_rate;  // sample rate in hertz
    int iscomplex; // complex sample flag
    int threadid; // which thread are we
    char stationid[3]; // which station are we.
    char exp_name[17]; // Experiment name
    char scan_name[17]; // scan_name
    size_t block_size; // size of one second of output data including headers
    size_t sizeof_buffer; // size of one second of 32bit complex beam data (no headers)
    size_t sizeof_beam; // size of 1 sample of 32bit complex beam data (no headers)
    float *b_scales; // bandpass mean
    float *b_offsets; // bandpass offset

    // observation info
    char telescope[24];
    char source[24];
    char obs_mode[8];
    double fctr;

    char ra_str[16];
    char dec_str[16];

    double BW;

    long double MJD_epoch;

    char date_obs[24];
    char basefilename[1024];

    int got_scales;

};

void vdif_write_data( struct vdifinfo *vf, int8_t *output );
void vdif_write_second( struct vdifinfo *vf, vdif_header *vhdr,
        float *data_buffer_vdif, float *gain );

void populate_vdif_header(
        struct vdifinfo *vf,
        vdif_header     *vhdr,
        char            *metafits,
        char            *obsid,
        char            *time_utc,
        int              sample_rate,
        long int         frequency,
        int              nchan, 
        long int         chan_width,
        char            *rec_channel,
        struct delays   *delay_vals );

ComplexFloat get_std_dev_complex( ComplexFloat *input, int nsamples );

void set_level_occupancy( ComplexFloat *input, int nsamples,
                          float *new_gain );

void get_mean_complex( ComplexFloat *input, int nsamples, float *rmean,
                       float *imean, ComplexFloat *cmean );

void normalise_complex( ComplexFloat *input, int nsamples, float scale );

void to_offset_binary( int8_t *i, int n );

#ifndef HAVE_CUDA
void invert_pfb_ifft( ComplexDouble ***detected_beam, int file_no,
                      int nsamples, int nchan, int npol,
                      float *data_buffer_vdif );

void invert_pfb_ord( ComplexDouble ***detected_beam, int file_no,
                      int nsamples, int nchan, int npol,
                      ComplexDouble **fils, int fil_size,
                      float *data_buffer_uvdif )
#endif

#endif
