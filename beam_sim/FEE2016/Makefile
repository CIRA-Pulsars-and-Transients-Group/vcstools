# C compiler
CC=cc
CXX=cc

# HDF5 directories
HDF5_INC=${HDF5_ROOT}/include

# BOOST directories
# -> need to do: module load boost/1.57.0
BOOST_INC=${BOOST_ROOT}/include


CXXFLAGS=-O3 -g -D_GET_JONES_2016_TEST -std=gnu++11 -lhdf5_cpp -lhdf5 -lboost_system -lboost_filesystem $(addprefix -I,$(HDF5_INC) $(BOOST_INC))
LDFLAGS=$(addprefix -L,$(HDF5_LIB) $(BOOST_LIB))

TARGETS = beam2016implementation.o mwa_beam_interface.o system.o

all: $(TARGETS)

clean: 
	$(RM) *.o
