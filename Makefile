# C compiler
CC=gcc

# CUDA compiler
NVCC=nvcc

# CUDA architecture (i.e. compute capability) to compile for
NVCC_ARCH=sm_35

# Additional non-standard header directories
INC_DIRS=

# Additional non-standard library directories
LIB_DIRS=

# C flags
CFLAGS=-O3 -lm -lsla -lcfitsio

# CUDA flag
CUDAFLAGS=-lcuda



all: pabeam_gpu

pabeam_gpu: pabeam.cu
	$(NVCC) -o $@ $^ -arch=$(NVCC_ARCH) $(LIB_DIRS) $(INC_DIRS) $(CFLAGS) $(CUDAFLAGS)

profile: pabeam_gpu
	nvprof ./$^ -f 184.96e6 -r "05:34:31.97" -d "+22:00:52.06" -t "2014-11-07T16:53:20" -m 1099414416_metafits_ppds.fits -b flagged_tiles.txt -x 0.01 -y 0.01


pabeam_cpu: pabeam.c
	$(CC) -Wall -o $@ $^ $(LIB_DIRS) $(INC_DIRS) $(CFLAGS)



