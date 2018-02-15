#ifndef FORM_BEAM_H
#define FORM_BEAM_H

#include "mycomplex.h"
#include "beam_common.h"

#define NPFB  4
#define NREC  16
#define NINC  4

#define REAL_NIBBLE_TO_UINT8(X)  ((X) & 0xf)
#define IMAG_NIBBLE_TO_UINT8(X)  (((X) >> 4) & 0xf)
#define UINT8_TO_INT(X)          ((X) >= 0x8 ? (signed int)(X) - 0x10 : (signed int)(X))
#define RE_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))))
#define IM_UCMPLX4_TO_FLT(X)  ((float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X))))
#define UCMPLX4_TO_CMPLX_FLT(X)  (CMaked((float)(UINT8_TO_INT(REAL_NIBBLE_TO_UINT8(X))), \
                                         (float)(UINT8_TO_INT(IMAG_NIBBLE_TO_UINT8(X)))))
#define DETECT(X)                (CReald(CMuld(X,CConjd(X))))



#ifdef HAVE_CUDA

void cu_form_beam( uint8_t *data, struct make_beam_opts *opts, ComplexDouble ***W,
                   ComplexDouble ****J, int file_no, int nstation, int nchan,
                   int npol, int outpol_coh, int outpol_incoh, double invw,
                   ComplexDouble ***detected_beam, float *coh, float *incoh );

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
