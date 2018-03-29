#ifndef MAKE_BEAM_SMALL_H
#define MAKE_BEAM_SMALL_H

#include <stdlib.h>
#include "beam_common.h"
#include "mycomplex.h"


#define MAX_COMMAND_LENGTH 1024

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

#endif
