# make\_beam\_small

----
## Unreleased

----
## v1.3.1

### Added

* GPU accelerated beamforming
* GPU accelerated PFB inversion
* Updated documentation describing the beamforming operations

### Removed

* Support for **-v** when compiled with GPU support

----
## v1.3.0

* Added support for full PFB inversion
* Added options for toggling output formats:
  - **-i** incoherently summed PSRFITS
  - **-p** coherently summed PSRFITS
  - **-v** PFB inverted (IFFT method) VDIF
  - **-u** PFB inverted (full, 'Ord' method) VDIF

