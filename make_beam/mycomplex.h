/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#ifndef MYCOMPLEX_H
#define MYCOMPLEX_H

/* Define how to handle complex numbers */
#ifdef HAVE_CUDA

#include <cuComplex.h>

#define  ComplexDouble  cuDoubleComplex
#define  ComplexFloat   cuFloatComplex

#define  CMakef(X,Y)    (make_cuFloatComplex(X,Y))
#define  CMaked(X,Y)    (make_cuDoubleComplex(X,Y))

#define  CAddf(X,Y)     (cuCaddf(X,Y))
#define  CSubf(X,Y)     (cuCsubf(X,Y))
#define  CMulf(X,Y)     (cuCmulf(X,Y))
#define  CDivf(X,Y)     (cuCdivf(X,Y))

#define  CRealf(X)      (cuCrealf(X))
#define  CImagf(X)      (cuCimagf(X))

#define  CAddd(X,Y)     (cuCadd(X,Y))
#define  CSubd(X,Y)     (cuCsub(X,Y))
#define  CMuld(X,Y)     (cuCmul(X,Y))
#define  CDivd(X,Y)     (cuCdiv(X,Y))

#define  CReald(X)      (cuCreal(X))
#define  CImagd(X)      (cuCimag(X))

#define  CConjf(X)      (cuConjf(X))
#define  CConjd(X)      (cuConj(X))

#define  CAbsf(X)       (cuCabsf(X))
#define  CAbsd(X)       (cuCabs(X))

#define  CExpf(X)       (CMakef(expf(CRealf(X))*cos(CImagf(X)), \
                                expf(CRealf(X))*sin(CImagf(X))))
#define  CExpd(X)       (CMaked(expf(CReald(X))*cos(CImagd(X)), \
                                expf(CReald(X))*sin(CImagd(X))))

#define  CSclf(X,F)     (CMakef(F*CRealf(X),F*CImagf(X)))
#define  CScld(X,F)     (CMaked(F*CReald(X),F*CImagd(X)))

#define  CRcpf(X)       (CSclf(CConjf(X),1.0/(CRealf(X)*CRealf(X) + CImagf(X)*CImagf(X))))
#define  CRcpd(X)       (CScld(CConjd(X),1.0/(CReald(X)*CReald(X) + CImagd(X)*CImagd(X))))

#define  CNegf(X)       (CSclf((X),-1.0))
#define  CNegd(X)       (CScld((X),-1.0))

#define  CD2F(X)        (CMakef((float)CReald(X),(float)CImagd(X)))
#define  CF2D(X)        (CMaked((double)CRealf(X),(double)CImagf(X)))

#else

#include <complex.h>
#define ComplexDouble  complex double
#define ComplexFloat   complex float

#define  CMakef(X,Y)    ((X)+(Y)*I)
#define  CMaked(X,Y)    ((X)+(Y)*I)

#define  CAddf(X,Y)     ((X)+(Y))
#define  CSubf(X,Y)     ((X)-(Y))
#define  CMulf(X,Y)     ((X)*(Y))
#define  CDivf(X,Y)     ((X)/(Y))

#define  CRealf(X)      (crealf(X))
#define  CImagf(X)      (cimagf(X))

#define  CAddd(X,Y)     ((X)+(Y))
#define  CSubd(X,Y)     ((X)-(Y))
#define  CMuld(X,Y)     ((X)*(Y))
#define  CDivd(X,Y)     ((X)/(Y))

#define  CReald(X)      (creal(X))
#define  CImagd(X)      (cimag(X))

#define  CConjf(X)      (conjf(X))
#define  CConjd(X)      (conj(X))

#define  CAbsf(X)       (cabsf(X))
#define  CAbsd(X)       (cabs(X))

#define  CExpf(X)       (cexpf(X))
#define  CExpd(X)       (cexp(X))

#define  CSclf(X,F)     ((F)*(X))
#define  CScld(X,F)     ((F)*(X))

#define  CRcpf(X)       (1.0/(X))
#define  CRcpd(X)       (1.0/(X))

#define  CNegf(X)       (-(X))
#define  CNegd(X)       (-(X))

#define  CD2F(X)        ((ComplexFloat)(X))
#define  CF2D(X)        ((ComplexDouble)(X))

#endif



#endif
