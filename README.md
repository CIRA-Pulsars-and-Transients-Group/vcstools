VCStools software tools
======


Installation
======
The installation is done in two steps. The first involves installing all the python scripts which is done with the command:
```bash
python setup.py build --build-scripts=<build_location>
```
Replacing <build_location> with where you would like to store the scripts.

The second step is to compile the beamformer which is much more difficult. All of the beamformer's dependancies must be taken into account as seen in this example cmake command:
```bash
cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_CUDA_COMPILER=$CUDA_COMPILER \
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} \
    -DCMAKE_CUDA_FLAGS=${CUDA_FLAGS} \
    -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} \
    -DFFTW3_ROOT_DIR=${FFTW3_ROOT_DIR} \
    -DFFTW3_INCLUDE_DIR=${FFTW_INCLUDE_DIR} \
    -DSLALIB_ROOT_DIR=${LIBSLA_LIB} \
    -DPSRFITS_UTILS_ROOT_DIR=${PSRFITS_UTILS_ROOT} \
    -DXGPU_ROOT_DIR=${XGPU_ROOT} \
    ..
```

For this reason we have created a docker image which is much easier to install and can be found [here](https://cloud.docker.com/u/cirapulsarsandtransients/repository/docker/cirapulsarsandtransients/vcstools)

You will have to make your own entry in config.py for your supercomputer which we are happy to help with.

Help
======
Documentation on how to process MWA VCS data can be found [here](https://wiki.mwatelescope.org/display/MP/Documentation)

Credit
======
If you use the MWA beamformer please give credit by citing:
[Ord et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...30O/abstract)

Status
=====
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/c5935a38ea9d4ec7a5c78fe107766cda)](https://www.codacy.com/manual/CIRA-Pulsars-and-Transients-Group/vcstools?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=CIRA-Pulsars-and-Transients-Group/vcstools&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/CIRA-Pulsars-and-Transients-Group/vcstools/badge.svg)](https://coveralls.io/github/CIRA-Pulsars-and-Transients-Group/vcstools)
[![Build Status](https://travis-ci.org/CIRA-Pulsars-and-Transients-Group/vcstools.svg?branch=master)](https://travis-ci.org/CIRA-Pulsars-and-Transients-Group/vcstools)
