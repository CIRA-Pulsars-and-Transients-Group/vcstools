#ifndef VDIFBEAM_H
#define VDIFBEAM_H

/* convenience type - this just collects all the vdif info together */
typedef struct vdifinfo_t {

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

} vdifinfo;

void vdif_write_second (vdifinfo *vf,int8_t *output);

#endif
