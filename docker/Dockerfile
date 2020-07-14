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


### Install vcstools ###
WORKDIR /tmp/vcstools-build
COPY . /tmp/vcstools-build
RUN pip3 install . && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
      -DSLALIB_ROOT_DIR=/usr \
      -DCMAKE_BEAMFORMER=1 \
      .. && \
    make install
# Once I make vcstools a proper package this PYTHONPATH addition won't be required
ENV PYTHONPATH="/usr/local/bin"
ENV CMD_VCS_DB_FILE None
CMD ["/bin/bash"]
