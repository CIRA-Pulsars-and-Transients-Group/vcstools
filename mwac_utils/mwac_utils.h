/*
 * mwac_utils.h
 *
 *  Created on: Jun 13, 2012
 *      Author: sord
 */

#ifndef MWAC_UTILS_H_
#define MWAC_UTILS_H_



#include <complex.h>
#include <stdlib.h>
// general stuff
#define MAX_POLS 4
#define STRSIZE 1024
#define VEL_LIGHT 299792458.0
#define BUFSIZE 4096
#define SHORTBUF 64
#define R2C_SIGN -1.0             // sign of the lm->uv Fourier transform
#define C2R_SIGN +1.0             // sign of the uv->lm Fourier transform
#define N_COPOL  2                // the number of polarised receptors on an antenna
#define MAX_FILTER_SIZE 32768     // For the beamformer

#define CASA_GAINS_FILE      1
#define MIRIAD_GAINS_FILE    2

struct filter_context {
   int ntaps;
   int nsamples;
   float filter[MAX_FILTER_SIZE];
};    


#ifdef __cplusplus
extern "C" {
#endif
    
    extern int nstation;
    extern int npol;
    extern int nfrequency;

    void mult2x2d(complex double *M1, complex double *M2, complex double *Mout);
    void mult2x2t(complex double *M1, complex double *M2, complex double *Mout);
    void mult2x2h(complex double *M1, complex double *M2, complex double *Mout);
    void multaccum2x2dd(complex double *M1, complex double *M2, complex double *Mout);
    void multaccum2x2dt(complex double *M1, complex double *M2, complex double *Mout);
    void multaccum2x2dh(complex double *M1, complex double *M2, complex double *Mout);
    void multaccum2x2hd(complex double *M1, complex double *M2, complex double *Mout);
    void cp2x2(complex double *Min, complex double *Mout);
    void inv2x2(complex double *Min, complex double *Mout);
    void mult2x2tlum(complex double *M1, complex double *M2, complex double *M3, complex double *Mout);
    double norm2x2(complex double *M, complex double *Mout);
    void conj2x2(complex double *M, complex double *Mout);
    
    int read_cal_file(complex double **G, int ninp, double *amp);
    int read_rts_file(complex double **G, complex double *M, int nant, double *amp, char *fname);
    int read_offringa_gains_file(complex double **antenna_gain, int nant, int coarse_chan, char *gains_file, int *order);
    int read_miriad_gains_file(char *fname, complex double **gains);
    int read_casa_gains_file(char *fname, complex double **gains,int nant, int chan_to_get);
    int gain_file_id(char *fname);
//    int calcEjones(complex response[MAX_POLS], const float freq, const float lat, const float az0, const float za0, const float az, const float za);
    
    
    void dec2hms(char *out, double in, int sflag);
    
    void fill_mapping_matrix();
    
    void get_baseline( int st1, int st2,int pol1, int pol2, complex float * data,
                      complex float * baseline);
    
    void get_baseline_lu( int st1, int st2,int pol1, int pol2, complex float * data,
                         complex float * baseline);
    
    
    void get_baseline_r(int st1, int st2, int pol1, int pol2, complex float * data,
                        complex float * baseline, int npol,int nstation, int nfrequency,int true_st1,int true_st2,int true_pol1,int true_pol2,int conjugate);
    
    void extractMatrix(complex float *matrix, complex float * packed);
    void extractMatrix_slow(complex float *matrix, complex float * packed);
    void full_reorder(complex float *full_matrix_h, complex float *reordered);
    
    // for the beamformer
    void invert_pfb(complex float * input, complex float *output, int nchan_in, int nchan_out, int npol_in, int npol_out,int mode,int last,void *context);
    int default_read_pfb_call(int in_fd, int out_fd, char *heap);

#ifdef __cplusplus
}
#endif

#endif /* MWAC_UTILS_H_ */
