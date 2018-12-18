/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef FORM_BEAM_H
#define FORM_BEAM_H

#include "mycomplex.h"
#include "beam_common.h"

#define NANT  128
#define NPFB  4
#define NREC  16
#define NINC  4

#define ANT2PFB(ant)    ((ant)>>5)                 /* pfb = ant/32           */
#define ANT2INC(ant)    (((ant)>>3)&0x03)          /* inc = (ant/8) % 4      */
#define AP2REC(ant,pol) ((((ant)<<1)+(pol))&0x0F)  /* rec = (2*ant+pol) % 16 */

/* structure for managing data arrays to be allocated on both host and device */
struct gpu_formbeam_arrays
{
    size_t coh_size;
    size_t incoh_size;
    size_t data_size;
    size_t Bd_size;
    size_t W_size;
    size_t J_size;
    ComplexDouble *W, *d_W;
    ComplexDouble *J, *d_J;
    ComplexDouble *Bd, *d_Bd;
    uint8_t *d_data;
    float   *d_coh;
    float   *d_incoh;
};


void malloc_formbeam( struct gpu_formbeam_arrays **g, unsigned int sample_rate,
        int nstation, int nchan, int npol, int outpol_coh, int outpol_incoh, int npointing );
void free_formbeam( struct gpu_formbeam_arrays **g );

/* Calculating array indices for GPU inputs and outputs */

#define D_IDX(s,c,a,p,nc)  ((s)         * (NINC*NREC*NPFB*(nc)) + \
                            (c)         * (NINC*NREC*NPFB)      + \
                            ANT2PFB(a)  * (NINC*NREC)           + \
                            AP2REC(a,p) * (NINC)                + \
                            ANT2INC(a))

#define W_IDX(p,a,c,pol,nc)   ((p) * (NPOL*(nc)*NANT)  + \
                               (a) * (NPOL*(nc))       + \
                               (c) * (NPOL)            + \
                               (pol))

#define J_IDX(p,a,c,p1,p2,nc) ((p) *  (NPOL*NPOL*(nc)*NANT) + \
                               (a)  * (NPOL*NPOL*(nc))      + \
                               (c)  * (NPOL*NPOL)           + \
                               (p1) * (NPOL)                + \
                               (p2))

#define B_IDX(p,s,c,pol,ns,nc) ((p)  * (NPOL*(nc)*(ns))   + \
                                (s)  * (NPOL*(nc))        + \
                                (c)  * (NPOL)             + \
                                (pol))
 
#define C_IDX(p,s,st,c,ns,nc)  ((p)  * ((nc)*NSTOKES*(ns)) + \
                                (s)  * ((nc)*NSTOKES)      + \
                                (st) *  (nc)               + \
                                (c))

#define I_IDX(s,c,nc)  ((s)*(nc) + (c))




#define W_FLAT_IDX(ch,ant,pol)   (((ch) <<8)         + \
                                 (((ant)<<1) & 0xC0) + \
                                 (((ant)<<2) & 0x38) + \
                                 (((ant)>>3) & 0x03) + \
                                 (((pol)<<2)))
#define J_FLAT_IDX(ch,ant,p1,p2) (((ch) <<9)          + \
                                 (((ant)<<2) & 0x180) + \
                                 (((ant)<<3) & 0x70)  + \
                                 (((ant)>>2) & 0x06)  + \
                                 (((p1 )<<3)))



/* Converting from 4+4 complex to full-blown complex doubles */

#define REAL_NIBBLE_TO_UINT8(X)  ((X) & 0xf)
#define IMAG_NIBBLE_TO_UINT8(X)  (((X) >> 4) & 0xf)
#define UINT8_TO_INT(X)          ((X) >= 0x8 ? (signed int)(X) - 0x10 : (signed int)(X))
#define RE_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))))
#define IM_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X))))
#define UCMPLX4_TO_CMPLX_FLT(X)  (CMaked((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))), \
                                         (float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X)))))
#define DETECT(X)                (CReald(CMuld(X,CConjd(X))))



#ifdef HAVE_CUDA

void cu_form_beam( uint8_t *data, struct make_beam_opts *opts, ComplexDouble ****W,
                   ComplexDouble *****J, int file_no, 
                   int npointing, int nstation, int nchan,
                   int npol, int outpol_coh, double invw, struct gpu_formbeam_arrays **g,
                   ComplexDouble ****detected_beam, float *coh, float *incoh );

#else

void form_beam( uint8_t *data, struct make_beam_opts *opts, ComplexDouble ***W,
                ComplexDouble ****J, int file_no, int nstation, int nchan,
                int npol, int outpol_coh, int outpol_incoh, double invw,
                ComplexDouble ***detected_beam, float *coh, float *incoh );

void form_stokes( ComplexDouble **detected_beam,
                  ComplexDouble noise_floor[][2][2],
                  int nchan, double invw, float *spectrum );

#endif

#endif
