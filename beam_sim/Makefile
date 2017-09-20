# C and C++ compiler
#CC = cc
#CXX = cc

# CUDA compiler
NVCC = nvcc

# CUDA architecture (i.e. compute capability) to compile for
NVCC_ARCH = sm_35

# Additional non-standard header/library directories
# cfitsio inc/lib dirs
CFITSIO_INC = ${CFITSIO_ROOT}/include
CFITSIO_LIB = ${CFITSIO_ROOT}/lib

# CUDA inc/lib dirs
CUDA_INC = ${CUDATOOLKIT_HOME}/include
CUDA_LIB = ${CUDATOOLKIT_HOME}/lib

# SLALIB inc/lib dirs
SLA_INC = ${PSRSOFT}/include
SLA_LIB = ${PSRSOFT}/lib

# HDF5 directories
HDF5_INC = ${HDF5_ROOT}/include
HDF5_LIB = ${HDF5_ROOT}/lib

# BOOST directories
# -> need to do: module load boost/1.57.0
BOOST_INC = ${BOOST_ROOT}/include
BOOST_LIB = ${BOOST_ROOT}/lib




###########################
# DON'T EDIT BELOW HERE * #
###########################
# * unless you know what you're doing...

INC_DIRS = $(addprefix -I,$(CFITSIO_INC) $(SLA_INC) $(HDF5_INC) $(BOOST_INC))
INC_DIRS_NOTILE = $(addprefix -I,$(CFITSIO_INC) $(SLA_INC))
LIB_DIRS = $(addprefix -L,$(CFITSIO_LIB) $(SLA_LIB) $(HDF5_LIB) $(BOOST_LIB))
LIB_DIRS_NOTILE = $(addprefix -L,$(CFITSIO_LIB) $(SLA_LIB))

# compiler flags
CFLAGS = -O3 -lm -lsla -lcfitsio #-lcuda
CXXFLAGS = -lhdf5_cpp -lhdf5 -lboost_system -lboost_filesystem -g -D_GET_JONES_2016_TEST #-std=gnu++11 
FLAGS = $(CFLAGS) $(CXXFLAGS)
FLAGS_NOTILE = $(CFLAGS)

# nvcc compile flags
NVCC_FLAGS = -arch=$(NVCC_ARCH)

# targets
FEE2016_TARGETS = FEE2016/beam2016implementation.o  FEE2016/mwa_beam_interface.o  FEE2016/system.o



all: pabeam_gpu_notile pabeam_gpu

pabeam_gpu_notile: pabeam_notile.cu
	$(NVCC) -o $@ $^ $(NVCC_FLAGS) $(LIB_DIRS_NOTILE) $(INC_DIRS_NOTILE) $(FLAGS_NOTILE)

pabeam_gpu: pabeam.cu $(FEE2016_TARGETS)
	$(NVCC) -o $@ $^ $(NVCC_FLAGS) $(LIB_DIRS) $(INC_DIRS) $(FLAGS)

#pabeam_cpu: pabeam.c
#	$(CC) -Wall -o $@ $^ $(LIB_DIRS) $(INC_DIRS) $(FLAGS)

FEE2016/%.o:
	@$(MAKE) --no-print-directory -C $(@D) $(@F)

#profile: pabeam_gpu
#	nvprof ./$^ -f 184.96e6 -r "05:34:31.97" -d "+22:00:52.06" -t "2014-11-07T16:53:20" -m 1099414416_metafits_ppds.fits -b flagged_tiles.txt -x 0.01 -y 0.01

clean:
	$(RM) pabeam_gpu_notile
	$(RM) pabeam_gpu
	$(RM) pabeam_cpu
