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


RUN pip3 install astropy \
    astroplan \
    matplotlib \
    numpy \
    psrqpy \
    scipy \
    skyfield \
    h5py \
    setuptools


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


# Install mwa_pb
WORKDIR /tmp/mwa_pb
RUN git clone https://github.com/MWATelescope/mwa_pb.git && \
    cd mwa_pb && \
    python3 setup.py install && \
    cd mwa_pb/data && \
    #if cerberus is down use the backup server
    wget -t 1 http://cerberus.mwa128t.org/mwa_full_embedded_element_pattern.h5; if [ $? == "4" ]; then wget http://orthrus.mwa128t.org/mwa_full_embedded_element_pattern.h5; fi

# Install mwa_client
WORKDIR /tmp/mwa_client
RUN git clone https://github.com/ICRAR/mwa-voltage.git && \
    cd mwa-voltage/mwa_pulsar_client && \
    python3 setup.py install

# Download psrcat database
WORKDIR /tmp/psrcat
RUN wget https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz  && \
    gunzip psrcat_pkg.tar.gz && \
    tar -xvf psrcat_pkg.tar
ENV PSRCAT_FILE /tmp/psrcat/psrcat_tar/psrcat.db

### Install vcstools ###
WORKDIR /tmp/vcstools-build
COPY . /tmp/vcstools-build
RUN mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
      -DSLALIB_ROOT_DIR=/usr \
      -DCMAKE_BEAMFORMER=1 \
      .. && \
    make install && \
    cd .. && \
    python3 setup.py build --build-scripts=/usr/local/bin
ENV PYTHONPATH="/usr/local/bin"
ENV CMD_VCS_DB_FILE None
ENV KNOWN_RFRB_CSV /tmp/vcstools-build/database/known_repeating_FRBs.csv
CMD ["/bin/bash"]
