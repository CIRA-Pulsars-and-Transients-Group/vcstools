# C and C++ compiler
CC = cc
CXX = cc

# CUDA compiler
NVCC = nvcc

# CUDA architecture (i.e. compute capability) to compile for
NVCC_ARCH = sm_35

# cfitsio inc/lib dirs
CFITSIO_INC = ${CFITSIO_ROOT}/include
CFITSIO_LIB = ${CFITSIO_ROOT}/lib

# CUDA inc/lib dirs
CUDA_INC = ${CUDATOOLKIT_HOME}/include
CUDA_LIB = ${CUDATOOLKIT_HOME}/lib

# SLALIB inc/lib dirs
SLA_INC = ${PSRSOFT}/include
SLA_LIB = ${PSRSOFT}/lib

# HDF5 inc/lib dirs
HDF5_INC = ${HDF5_ROOT}/include
HDF5_LIB = ${HDF5_ROOT}/lib

# BOOST inc/lib dirs
BOOST_INC = ${BOOST_ROOT}/include
BOOST_LIB = ${BOOST_ROOT}/lib


# full lists of inc/lib flags
INC_DIRS_NOTILE = $(addprefix -I,$(CFITSIO_INC) $(SLA_INC))
INC_DIRS = $(addprefix -I,$(CFITSIO_INC) $(SLA_INC) $(HDF5_INC) $(BOOST_INC))

LIB_DIRS_NOTILE = $(addprefix -L,$(CFITSIO_LIB) $(SLA_LIB))
LIB_DIRS = $(addprefix -L,$(CFITSIO_LIB) $(SLA_LIB) $(HDF5_LIB) $(BOOST_LIB))

# compiler flags
FLAGS_NOTILE = -O3 -lm -lsla -lcfitsio #-lcuda
FLAGS = $(FLAGS_NOTILE) -lhdf5_cpp -lhdf5 -lboost_system -lboost_filesystem #-g -D_GET_JONES_2016_TEST -std=gnu++11 

# nvcc compile flags
NVCC_FLAGS = -arch=$(NVCC_ARCH)

# targets
FEE2016_TARGETS = FEE2016/beam2016implementation.o FEE2016/mwa_beam_interface.o FEE2016/system.o
PABEAM_TARGETS = pabeam_kernal.o pabeam.o

TARGETS = $(PABEAM_TARGETS) #$(FEE2016_TARGETS)

all: pabeam_gpu_notile pabeam_gpu

pabeam_gpu_notile: pabeam_notile.cu
	$(NVCC) -o $@ $^ $(NVCC_FLAGS) $(LIB_DIRS_NOTILE) $(INC_DIRS_NOTILE) $(FLAGS_NOTILE)

pabeam_gpu: pabeam_kernal.cu pabeam.cu 
	$(NVCC) -c $(NVCC_FLAGS) $(LIB_DIRS) $(INC_DIRS) $(FLAGS) pabeam_kernal.cu
	$(NVCC) -c $(NVCC_FLAGS) $(LIB_DIRS) $(INC_DIRS) $(FLAGS) pabeam.cu
	$(NVCC) -o $@ $(TARGETS) $(NVCC_FLAGS) $(LIB_DIRS) $(INC_DIRS) $(FLAGS)

FEE2016/%.o:
	@$(MAKE) --no-print-directory -C $(@D) $(@F)

clean:
	$(RM) *.o
	$(RM) pabeam_gpu_notile
	$(RM) pabeam_gpu
	$(RM) pabeam_cpu
