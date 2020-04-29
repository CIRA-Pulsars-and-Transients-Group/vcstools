FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

# Add aarnet mirror to speed up package update
RUN perl -p -i.orig -e \
      's/archive.ubuntu.com/mirror.aarnet.edu.au\/pub\/ubuntu\/archive/' /etc/apt/sources.list \
      && sed -i '0,/# deb-src/{s/# deb-src/deb-src/}' /etc/apt/sources.list

# Install package dependencies
RUN apt-get update \
      && apt-get install -y \
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
      && apt-get clean all \
      && rm -r /var/lib/apt/lists/*


### Build psrfits_utils ###
ARG PSR_CONFIG_OPTS="--prefix=/usr"
WORKDIR /tmp/psrfits-build
RUN git clone https://github.com/demorest/psrfits_utils.git \
      && cd psrfits_utils \
      && echo "m4_ifdef([AM_PROG_AR], [AM_PROG_AR])" >> configure.ac \
      && autoreconf --install --force \
      && ./configure ${PSR_CONFIG_OPTS} \
      && make \
      && make install


### Build pulsar ###
WORKDIR /tmp/pulsar-build
ADD docker/pulsar_src.tar.gz .
RUN cd pulsar_src \
      && make INSTALL_DIR=/usr


### Install CUDA ###
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
WORKDIR /tmp/cuda-build
RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin && \
    mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600 && \
    wget http://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda-repo-ubuntu1804-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb && \
    dpkg -i cuda-repo-ubuntu1804-10-2-local-10.2.89-440.33.01_1.0-1_amd64.deb && \
    apt-key add /var/cuda-repo-10-2-local-10.2.89-440.33.01/7fa2af80.pub && \
    apt-get update && \
    apt-get -y install cuda



### Install xGPU ###
ARG XGPU_MAKE_OPTS="NSTATION=128 NFREQUENCY=128 NTIME=1000 NTIME_PIPE=100"
ARG XGPU_PREFIX="/usr"
WORKDIR /tmp/xgpu-build
RUN git clone https://github.com/GPU-correlators/xGPU.git \
      && cd xGPU/src \
      && make ${XGPU_MAKE_OPTS} \
      && make install prefix=${XGPU_PREFIX}


### Install vcstools ###
WORKDIR /tmp/vcstools-build
COPY . /tmp/vcstools-build
RUN ls
RUN mkdir build \
      && cd build \
      && cmake -DCMAKE_INSTALL_PREFIX=/usr \
      -DCMAKE_CUDA_COMPILER=/usr/local/cuda/bin/nvcc \
      -DSLALIB_ROOT_DIR=/usr \
      -DCMAKE_BEAMFORMER=1 \
      .. \
      && make install

WORKDIR /
RUN rm -rf /tmp/*

CMD ["/bin/bash"]
