# Small-Grid Dithers

This is my scripts used to generate PSF, add noise to them, apply the LOCI
algorithm, and calculate the azimuthal average and contrast curves.  The first
three steps were originally written in Mathematica code and I have translated
(and made them run faster) into Python.

The overall procedure is summarized in the Figure below and the Python scripts
below are used sequentially to arrive at the final average 5-sigma contrast curves.

Note that some input data to the scripts below such as filter transmission, QE
curves, and/or background noise levels as well as Brian York's JWST Simulator
python modules (for AddNoisetoPSF.py) will only be provided upon request.

![alt text](https://github.com/Skyhawk172/SmallGridDithers/blob/master/SGDcartoon.png "SGD cartoon")


## 1. Generate_PSFs.py

Generate PSF grids using WebbPSF.

The idea here is that we generate one science target PSF using the expected
target acquisition accuracy, then simulate the acquision of another reference
PSF that we dither behind the coronagraphic mask, following a small square grid
pattern. Typically, the square grids have 9 points and the grid step is 20 mas,
although both can be modified on the command line (see below).  The PSF are all
generated using [WebbPSF](http://www.stsci.edu/jwst/software/webbpsf).

Command-line arguments (-- for optional):

| Argument     | Name      | Description                                      |
|--------------|-----------|--------------------------------------------------|
|-h, --help    |           | Show this help message and exit                  |
|-I            |INSTRUMENT | Instrument to use (MIRI or NIRCam)               |
|-f            |FILTER     | Filter to use                                    |
|-mask         |MASK_CORON | Mask coron. to use                               |
|-stop         |PUPIL_STOP | Pupil stop to use                                |
|--rms         |RMS        | Select OPD rms to use, if available (default: 136 nm)|
|--noopd       |           | Do not use any OPD (default: False)              |
|--jitter      |JITTER     | Sigma Jitter (default: 0 mas)                    |
|--fov         |FOV        | Field of view (diam.; default=7.04 arcsecond)    |
|--nruns       |NRUNS      | Number of SGD grids to generate (default: 1)     |
|--gstep       |GSTEP      | SGD grid steps (default: 20 mas)                 |
|--gnpts       |GNPTS      | SGD square grid points (default: 9)              |


## 2. AddNoisetoPSF.py

Add simple noise sources to pre-generated PSFs. Included are approximations for
stellar photon noise, background photon noise, detector noise (readnoise, dark
noise), as well as OTE throughput and detector quantum efficiency. The script
scales the original PSF (output of generate_PSFs.py) according to the user's
inputs regarding the stellar spectrum, distance, and exposure time and saves the
output images in a separate folder.

Note that the image size for MIRI and NIRCam (Short wave & Long wave) are
hardcoded in the function "prepareSpectrum". The default values are 64, 109, and
222 pixels across for MIRI, NIRCam short wave, and NIRCam long wave,
respectively. The user may need to change these values accordingly.

Command-line arguments (-- for optional):

| Argument     | Name      | Description                                      |
|--------------|-----------|--------------------------------------------------|
|-h, --help    |           | show this help message and exit                  |
|-I            |Instrument | JWST instrument                                  | 
|-f            |FILTER     | filter to use                                    |
|-step         |STEP       | dither grid step size (default: 20 mas)          |
|-jitter       |JITTER     | sigma Jitter (default: 0 mas)                    |
|--webbpsf     |           | use if images were created with WebbPSF          |
|--t           |T          | effective temperature (default: Solar T=5800 K)  |
|--z           |z          | metallicity (default: Solar z=0)                 |
|--g           |G          | surface gravity (default: Solar log_g=4.44)      |
|--R           |R          | radius of star (default: 1 Rsun )                |
|--dist        |DIST       | distance to star parsec (default: 10 pc)         |
|--exptime     |EXPTIME    | exposure time (default: 1 sec)                   |
|--run         |RUN        | run number (default: all)                        |
|--vega        |           | use Vega spectrum (default: False)               |
|--clobber     |           | overwrite output files (default: False)          |
|--nonoise     |           | to turn off noise sources (default: add noise)   |
|--rms         |RMS        | select which OPD rms to use, if available (default: 400 nm)|
|--planets     |           | Inject HR8799 planets (default: False)           |
|--roll        |           | Roll the planets clockwise in degrees (default: 0 deg.)|

Note: the --webbpsf flag is used for MIRI PSF generated with WebbPSF so that the
Quantum efficiency as well as the Germanium and OTE transmission profiles can be
applied. For MIRI PSF create with Mathematica or NIRCam PSF generated with
WebbPSF, this is not necessary as these contributions are taken into account in
the filter profiles included.

## 3. (ApplyLOCI.py)

See https://github.com/Skyhawk172/PythonLOCI for now...

## 4. AzimuthalAverage_SGD.py

Ingests the all outputs "Maps" of ApplyLOCI.py and calculates an average 5-sigma
contrast curve along with a +/- 1-sigma range.

Contrast is simply defined as the... 
