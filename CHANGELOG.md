# VCS Tools

----
## Unreleased

* Remove deprecated code for the beamformer
* Rename `make_beam_small` to `make_beam` 

----
## v1.3.2 (2018-02-19)

### Added

* Imported the offline correlator into this project

----
## v1.3.1 (2018-02-15)

### Added

* GPU accelerated beamforming
* GPU accelerated PFB inversion
* Updated documentation describing the beamforming operations

### Removed

* Support for **-v** when compiled with GPU support

----
## v1.3.0 (2018-02-03)

* Added support for full PFB inversion
* Added options for toggling output formats:
  - **-i** incoherently summed PSRFITS
  - **-p** coherently summed PSRFITS
  - **-v** PFB inverted (IFFT method) VDIF
  - **-u** PFB inverted (full, 'Ord' method) VDIF

