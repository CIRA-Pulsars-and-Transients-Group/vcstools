VCStools software tools
======
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/77cdb25072364ee48dd2f1e3ca078af5)](https://www.codacy.com/gh/CIRA-Pulsars-and-Transients-Group/vcstools/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=CIRA-Pulsars-and-Transients-Group/vcstools&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/77cdb25072364ee48dd2f1e3ca078af5)](https://www.codacy.com/gh/CIRA-Pulsars-and-Transients-Group/vcstools/dashboard?utm_source=github.com&utm_medium=referral&utm_content=CIRA-Pulsars-and-Transients-Group/vcstools&utm_campaign=Badge_Coverage)
[![Build Status](https://travis-ci.org/CIRA-Pulsars-and-Transients-Group/vcstools.svg?branch=master)](https://travis-ci.org/CIRA-Pulsars-and-Transients-Group/vcstools)
[![Documentation Status](https://readthedocs.org/projects/mwa-vcstools/badge/?version=latest)](https://mwa-vcstools.readthedocs.io/en/latest/?badge=latest)

Installation
------
The installation is done in two steps. The first involves installing all the python scripts, which is done with the command:
```bash
pip install mwa-vcstools
```
Or git clone into the repo and run
```bash
python setup.py install
```
or
```bash
python3 setup.py install --prefix="<install_dir>" --single-version-externally-managed --record=record.txt
```

The second step is to compile the beamformer which is much more difficult. All of the beamformer's dependancies must be taken into account as seen in this example cmake command:
```bash
cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_CUDA_COMPILER=$CUDA_COMPILER \
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} \
    -DCMAKE_CUDA_FLAGS=${CUDA_FLAGS} \
    -DCFITSIO_ROOT_DIR=${MAALI_CFITSIO_HOME} \
    -DFFTW3_ROOT_DIR=${FFTW3_ROOT_DIR} \
    -DFFTW3_INCLUDE_DIR=${FFTW_INCLUDE_DIR} \
    -DPAL_ROOT_DIR=${PAL_ROOT} \
    -DPSRFITS_UTILS_ROOT_DIR=${PSRFITS_UTILS_ROOT} \
    ..
```

For this reason, we have created a docker image which is much easier to install and can be found [here](https://cloud.docker.com/u/cirapulsarsandtransients/repository/docker/cirapulsarsandtransients/vcstools)

You will have to make your own entry in vcstools/config.py for your supercomputer which we are happy to help with.

Help
------
Documentation on how to process MWA VCS data can be found [here](https://wiki.mwatelescope.org/display/MP/Documentation) and the developer documentation can be found [here](https://mwa-vcstools.readthedocs.io)

Credit
------
You can reference this repository using: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3762792.svg)](https://doi.org/10.5281/zenodo.3762792)

If you use the MWA beamformer, please give credit by citing:
[Ord et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...30O/abstract)

If you used polarimetry, please give credit by citing: 
[Xue et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019PASA...36...25X/abstract)

If you used the inverse PFB, please give credit by citing:
[McSweeney et al. (2020)](http://dx.doi.org/10.1017/pasa.2020.24)
