==============================================
Murchison WideField Array Correlator Utilities
==============================================

Contents:

1) - build_lfiles: A utility that will read the output format of the MWA correlator and generate the LFILE format as used by the MWA32T prototype
2) - check_corr: A utility that operates on a raw correlator dump to plot antenna correlation products (REQUIRES TWOPIP LIBRARY)


Build Instructions

Uses cmake (version 2.8). It will try and find the requirements for the utilities. They are currently:

build_lfiles --- NO DEPENDENCIES - will build from this tarball
check_corr   --- TWOPIP library from the MWA git repository (the TWOPIP library itself has many dependencies)

CMAKE instructions:

mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=<OPTIONAL> <path to source>
make -j4 && make install


