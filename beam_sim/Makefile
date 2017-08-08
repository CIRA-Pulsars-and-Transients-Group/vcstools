# C compiler
CC=cc

# CUDA compiler
NVCC=nvcc

# CUDA architecture (i.e. compute capability) to compile for
NVCC_ARCH=sm_35

# Additional non-standard header/library directories
# cfitsio inc/lib dirs
CFITSIO_INC=${CFITSIO_ROOT}/include
CFITSIO_LIB=${CFITSIO_ROOT}/lib

# CUDA inc/lib dirs
CUDA_INC=${CUDATOOLKIT_HOME}/include
CUDA_LIB=${CUDATOOLKIT_HOME}/lib

# SLALIB inc/lib dirs
SLA_INC=${PSRSOFT}/include
SLA_LIB=${PSRSOFT}/lib

# Additional C flags
CFLAGS=-O3



###########################
# DON'T EDIT BELOW HERE * #
###########################
# * unless you know what you're doing...

INC_DIRS=$(addprefix -I,$(CFITSIO_INC) $(SLA_INC))
LIB_DIRS=$(addprefix -L,$(CFITSIO_LIB) $(SLA_LIB))

FLAGS=$(CFLAGS) -lm -lsla -lcfitsio
#CUDAFLAGS=-lcuda

# nvcc compile flags
NVCC_FLAGS=-arch=$(NVCC_ARCH)



all: pabeam_gpu

pabeam_gpu: pabeam.cu
	$(NVCC) -o $@ $^ $(NVCC_FLAGS) $(LIB_DIRS) $(INC_DIRS) $(FLAGS)

profile: pabeam_gpu
	nvprof ./$^ -f 184.96e6 -r "05:34:31.97" -d "+22:00:52.06" -t "2014-11-07T16:53:20" -m 1099414416_metafits_ppds.fits -b flagged_tiles.txt -x 0.01 -y 0.01

pabeam_cpu: pabeam.c
	$(CC) -Wall -o $@ $^ $(LIB_DIRS) $(INC_DIRS) $(FLAGS)

clean:
	rm -f pabeam_gpu
	rm -f pabeam_cpu

