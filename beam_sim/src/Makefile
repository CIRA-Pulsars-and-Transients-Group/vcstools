# CUDA compiler
CC = nvcc

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
FLAGS = -O3 # -D_GET_JONES_2016_TEST

# nvcc compiler flags
NVCC_FLAGS = -arch=$(NVCC_ARCH) -Xcompiler

# All compile flags
CFLAGS = $(NVCC_FLAGS) $(FLAGS)
INCFLAGS = $(INC_DIRS)
LDFLAGS = $(LIB_DIRS) -lm -lsla -lcfitsio -lhdf5_cpp -lhdf5 -lboost_system -lboost_filesystem # -lcuda

# targets
FEE2016_OBJECTS = FEE2016/beam2016implementation.o FEE2016/mwa_beam_interface.o FEE2016/system.o
PABEAM_OBJECTS = utils.o pabeam_kernal.o pabeam.o

OBJECTS = $(FEE2016_OBJECTS) $(PABEAM_OBJECTS)

#TARGETS = pabeam_gpu_notile pabeam_gpu
TARGETS = pabeam_gpu

BIN_DIR = ../bin

all: $(TARGETS)


pabeam_gpu_notile: pabeam_notile.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)


pabeam_gpu: $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)


%.o: %.cu
	$(CC) $(CFLAGS) -c $@ $< $(INCFLAGS)


FEE2016/%.o:
	@$(MAKE) -C $(@D) $(@F)


clean:
	$(RM) *.o
	$(RM) pabeam_gpu_notile
	$(RM) pabeam_gpu
	$(RM) pabeam_cpu


# clean this directory and the FEE2016 model directory
dist-clean:
	$(RM) *.o
	$(RM) pabeam_gpu_notile
	$(RM) pabeam_gpu
	$(RM) pabeam_cpu
	@$(MAKE) -C FEE2016/ clean

install:
	mkdir -p $(BIN_DIR)
	cp $(TARGETS) $(BIN_DIR)
