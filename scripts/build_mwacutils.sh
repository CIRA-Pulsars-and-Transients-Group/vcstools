#!/bin/bash -l

# This script will build and install the version of MWACUtils specified in
# the first argument ($1) into the directory specified in the second
# argument ($2). The locations of the dependencies must be entered into
# this script "by hand".

# Set (and create, if necessary) an install folder for this remote branch
src_dir=`readlink -f "$1"`
cd "$src_dir"

install_dir=`readlink -f "$2"`
mkdir -p "${install_dir}"

# Set up the necessary dependencies for cmake
root_dir=/group/mwaops/PULSAR

cmake_install_prefix=${install_dir}

slalib_root_dir=${root_dir}/src/slalib_c

cfitsio_root_dir=/pawsey/cle52up04/devel/PrgEnv-gnu/5.2.82/gcc/4.8.2/sandybridge/cfitsio/3370
cfitsio_library=${cfitsio_root_dir}/lib/libcfitsio.a 
cfitsio_include_dir=${cfitsio_root_dir}/include 

fftw3_root_dir=/opt/cray/fftw/3.3.4.3/sandybridge
mpi_include_path=/opt/cray/mpt/7.0.0/gni/mpich2-gnu/48/include 
psrfits_utils_library=${root_dir}/lib/libpsrfits_utils.a
psrfits_utils_include_dir=${root_dir}/include/psrfits_utils

# Load the necessary modules
module swap PrgEnv-cray PrgEnv-gnu 
module swap craype-ivybridge craype-sandybridge
module load cmake cfitsio cudatoolkit
#module swap fftw/3.3.4.0 fftw/3.3.4.3

export CRAYPE_LINK_TYPE=dynamic

rm -r build
mkdir build
cd build

# Build
CC=cc CXX=CC cmake -DCMAKE_INSTALL_PREFIX=${cmake_install_prefix} \
                   -DCFITSIO_LIBRARY=${cfitsio_library} \
                   -DCFITSIO_INCLUDE_DIR=${cfitsio_include_dir} \
                   -DFFTW3_ROOT_DIR=${fftw3_root_dir} \
                   -DMPI_INCLUDE_PATH=${mpi_include_path} \
                   -DSLALIB_ROOT_DIR=${slalib_root_dir} \
                   -DPSRFITS_UTILS_LIBRARY=${psrfits_utils_library} \
                   -DPSRFITS_UTILS_INCLUDE_DIR=${psrfits_utils_include_dir} \
                   ..

make
make install
