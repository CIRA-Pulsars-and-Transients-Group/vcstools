FROM nvidia/cuda:10.1-devel-ubuntu18.04

# Install package dependencies
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y \
         autoconf \
         autotools-dev \
         build-essential \
         cmake \
         dirmngr \
         git \
         libboost-all-dev \
         libcfitsio-bin \
         libcfitsio-dev \
         libfftw3-dev \
         libhdf5-dev \
         libtool \
         wget \
         python3 \
         python3-pip \
         python3-setuptools \
         curl && \
    apt-get clean all && \
    rm -r /var/lib/apt/lists/*


### Build psrfits_utils ###
ARG PSR_CONFIG_OPTS="--prefix=/usr"
WORKDIR /tmp/psrfits-build
RUN git clone https://github.com/demorest/psrfits_utils.git && \
    cd psrfits_utils && \
    echo "m4_ifdef([AM_PROG_AR], [AM_PROG_AR])" >> configure.ac && \
    autoreconf --install --force && \
    ./configure ${PSR_CONFIG_OPTS} && \
    make && \
    make install


### Build pulsar ###
WORKDIR /tmp/pulsar-build
RUN wget https://github.com/CIRA-Pulsars-and-Transients-Group/vcstools/raw/master/docker/pulsar_src.tar.gz && \
    tar -zxvf pulsar_src.tar.gz && \
    cd pulsar_src && \
    make INSTALL_DIR=/usr


### Install xGPU ###
ARG XGPU_MAKE_OPTS="NSTATION=128 NFREQUENCY=128 NTIME=1000 NTIME_PIPE=100"
ARG XGPU_PREFIX="/usr"
WORKDIR /tmp/xgpu-build
RUN git clone https://github.com/GPU-correlators/xGPU.git && \
    cd xGPU/src && \
    make ${XGPU_MAKE_OPTS} && \
    make install prefix=${XGPU_PREFIX}


### Install rust and cargo
ENV PATH="/root/.cargo/bin:${PATH}"
ENV HDF5_INCLUDE_PATH = "/usr/include/hdf5/serial/"
RUN curl https://sh.rustup.rs -sSf |  bash -s -- -y


### Install hyperbeam
ENV LD_LIBRARY_PATH="/tmp/hyperbeam/mwa_hyperbeam/include:/tmp/hyperbeam/mwa_hyperbeam/lib:${LD_LIBRARY_PATH}"
ENV HYPERBEAM_ROOT="/tmp/hyperbeam/mwa_hyperbeam/"
ENV HYPERBEAM_INCLUDE_DIR="/tmp/hyperbeam/mwa_hyperbeam/include"
ENV HYPERBEAM_LIB="/tmp/hyperbeam/mwa_hyperbeam/lib"
WORKDIR /tmp/hyperbeam
RUN git clone https://github.com/MWATelescope/mwa_hyperbeam.git && \
    cd mwa_hyperbeam && \
    wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5 && \
    cargo test && \
    cargo build --release && \
    mkdir -p lib && \
    cp target/release/libmwa_hyperbeam.so target/release/libmwa_hyperbeam.a lib/

### Install vcstools ###
WORKDIR /tmp/vcstools-build
COPY . /tmp/vcstools-build
RUN pip3 install . && \
    mkdir -p build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
      -DSLALIB_ROOT_DIR=/usr \
      -DCMAKE_BEAMFORMER=1 \
      -DHYPERBEAM_ROOT=${HYPERBEAM_ROOT} \
      -DHYPERBEAM_HDF5_DIR=/tmp/hyperbeam/mwa_hyperbeam/ \
      -DCTESTS=1 \
      .. && \
    make install
# Once I make vcstools a proper package this PYTHONPATH addition won't be required
ENV PYTHONPATH="/usr/local/bin"
CMD ["/bin/bash"]
