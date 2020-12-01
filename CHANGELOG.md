# VCS Tools

----
## v1.5.0 (2019-06-28)

Upgraded the beamforming to be multi-pixel so hundreads of tied-array beams can be created at once. This increases the speed of beamforming by a factor of 5.

----
## v1.4.3 (2019-06-20)

Upgraded all scripts to python 3

----
## v1.4.2 (2019-04-18)

*   Removed all functions that use MWA_Tools. Since we now use astropy to more accurately calculate altitude and azimuth, the results differ ~0.1 degrees from the MWA_Tools calculations. This leads to a ~1% difference in tile beam power calculations.

----
## v1.4.1 (2018-03-28)

*   Renamed `make_beam_small` to `make_beam`

----
## v1.4.0 (2018-03-26)

-   Remove deprecated code for the beamformer
-   Make `splice_psrfits` overwrite existing output file instead of aborting if output file already exists
-   Include the (abbreviated) git hash of the current build in the version number
-   Add a -V option (version) to all programs and scripts within VCS Tools

----
## v1.3.2 (2018-02-19)

-   Imported the offline correlator into this project

----
## v1.3.1 (2018-02-15)

### Added

-   GPU accelerated beamforming
-   GPU accelerated PFB inversion
-   Updated documentation describing the beamforming operations

### Removed

-   Support for **-v** when compiled with GPU support

----
## v1.3.0 (2018-02-03)

-   Added support for full PFB inversion
-   Added options for toggling output formats:
  - **-i** incoherently summed PSRFITS
  - **-p** coherently summed PSRFITS
  - **-v** PFB inverted (IFFT method) VDIF
  - **-u** PFB inverted (full, 'Ord' method) VDIF