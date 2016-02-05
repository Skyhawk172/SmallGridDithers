# Small-Grid Dithers

This is my scripts used to generate PSF, add noise to them, apply the LOCI
algorithm, and calculate the azimuthal average and contrast curves.  The first
three steps were originally written in Mathematica code and I have translated
(and made them run faster) into Python.

The overall procedure is summarized in the Figure below:

![alt text](https://github.com/Skyhawk172/SmallGridDithers/blob/master/SGDcartoon.png "SGD cartoon")


## 1. Generate_PSFs.py

usage: generate_PSFs.py [-h] -I INSTRUMENT -f FILTER -mask MASK_CORON -stop
                        PUPIL_STOP [-gridstep GRIDSTEP] [-gridNpts GRIDNPTS]
                        [-jitter JITTER] [-rms RMS] [-nruns NRUNS] [-fov FOV]
                        [--noopd]

Generate PSF grids for LOCI using WebbPSF

optional arguments:
  -h, --help          show this help message and exit
  -I INSTRUMENT       instrument to use (MIRI or NIRCam)
  -f FILTER           Filter to use
  -mask MASK_CORON    Mask coron. to use
  -stop PUPIL_STOP    Pupil stop to use
  -gridstep GRIDSTEP  SGD grid steps (default: 20)
  -gridNpts GRIDNPTS  SGD square grid points (default: 9)
  -jitter JITTER      sigma Jitter (default: 0 mas)
  -rms RMS            select OPD rms to use, if available (default: 136 nm)
  -nruns NRUNS        number of SGD grids to generate (default: 1)
  -fov FOV            Field of view (diameter) in arcseconds (default=7.04)
  --noopd             don't use any OPD (default: False)



## 2. AddNoisetoPSF.py

usage: AddNoisetoPSF.py [-h] -I Instrument -f FILTER -step STEP -jitter JITTER
                        [--t T] [--z Z] [--g G] [--R R] [--dist DIST]
                        [--exptime EXPTIME] [--run RUN] [--vega] [--clobber]
                        [--nonoise] [--rms RMS]

Add simple noise sources to pre-generated PSFs

optional arguments:
  -h, --help         show this help message and exit
  -I Instrument      JWST instrument
  -f FILTER          Filter to use
  -step STEP         Dither grid step size (default: 20 mas)
  -jitter JITTER     sigma Jitter (default: 0 mas)
  --t T              effective temperature (default: solar T_eff=5800 K)
  --z Z              metallicity (default: solar z=0)
  --g G              surface gravity (default: solar log_g=4.44)
  --R R              radius of star (default: 1 Rsun )
  --dist DIST        distance to star parsec (default: 10 pc)
  --exptime EXPTIME  exposure time (default: 1 sec)
  --run RUN          run number (default: all)
  --vega             use Vega spectrum (default: False)
  --clobber          overwrite output files (default: False)
  --nonoise          to turn off noise sources (default: add noise)
  --rms RMS          select which OPD rms to use, if available (default: 400
                     nm)

## 3. (ApplyLOCI.py)


## 4. AzimuthalAverage_SGD.py


