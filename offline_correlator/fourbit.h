#ifndef __FOURBIT_H
#define __FOURBIT_H

#include <stdint.h>

	
#ifdef __cplusplus
extern "C" {
#endif


void build_eight_bit_lookup();

void process_MWA_buffer(void *input,void *output);
void process_MWA_atomic(void *input,void *output, const size_t stride, size_t offset);
void process_MWA_atomic_raw(void *input,void *output, const size_t stride, size_t offset);

void float_to_4bit(int nsamps, float *in, int8_t *out);
void expand_4bit(uint16_t *input, int8_t *output);
int8_t split_8bit(int8_t *input);

#ifdef __cplusplus
}
#endif


#endif
