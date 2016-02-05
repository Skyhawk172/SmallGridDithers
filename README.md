# Small-Grid Dithers

This is my scripts used to generate PSF, add noise to them, apply the LOCI
algorithm, and calculate the azimuthal average and contrast curves.  The first
three steps were originally written in Mathematica code and I have translated
(and made them run faster) into Python.

The overall procedure is summarized in the Figure below:

Reference-style: 
![alt text][logo]

[logo]: https://github.com/Skyhawk172/SmallGridDithers/blob/master/SGDcartoon.pdf "SGD cartoon"


1. Generate_PSFs.py



2. AddNoisetoPSF.py


3. (ApplyLOCI.py)


4. AzimuthalAverage_SGD.py


