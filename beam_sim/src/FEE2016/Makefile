# C++ compiler
CXX=cc

# HDF5 directories
HDF5_INC=${HDF5_ROOT}/include

# BOOST directories
# -> need to do: module load boost/1.57.0
BOOST_INC=${BOOST_ROOT}/include


CXXFLAGS=-O3 -D_MWA_2016_BEAM_MODEL_ -std=gnu++11 $(addprefix -I,$(HDF5_INC) $(BOOST_INC))

TARGETS = beam2016implementation.o mwa_beam_interface.o system.o

all: $(TARGETS)

clean: 
	$(RM) *.o
