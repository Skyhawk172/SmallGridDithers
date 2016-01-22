import numpy as np
import scipy.interpolate as scint
import os, pyfits, sys
import scipy.optimize as scopt


def interpolate(c,dim,targ,interp1,interp2):
    ref =np.zeros( [dim,dim] )
    mask=np.zeros( [dim,dim] )
    for i in xrange(dim):
        for j in xrange(dim):
            ref[i,j]  = c[2]*(interp1(i-c[0],j-c[1]))
            mask[i,j] = interp2(i-c[0],j-c[1])

    return np.log10( np.sum( ((ref-targ)*mask)**2 ) )


def align_images(run,dim,reference,target):
    radius=25
    mask=np.ones( [dim,dim] )
    x = np.arange(0,dim)
    y = np.arange(0,dim)
    xmid=dim/2
    ymid=dim/2

    mask=np.array( [ [0 if (np.sqrt( (i-xmid)**2 + (j-ymid)**2) < radius) else 1 for i in xrange(dim)] for j in xrange(dim)] )


    alpha=0.
    beta =0.
    nu   =1.
    guess=[alpha,beta,nu]

    interp_ref =scint.RectBivariateSpline(y,x,reference)
    interp_mask=scint.RectBivariateSpline(y,x,mask,kx=1,ky=1)

    results= scopt.minimize(interpolate,guess,args=(dim,target,interp_ref,interp_mask),tol=1e-5)
    print "  Run %d: alpha=%10.8f  beta=%10.8f  nu= %10.8f" %(run,results.x[0],results.x[1],results.x[2])

    
#    ref_aligned=np.array( [ [results.x[2]*(interp_ref(i-results.x[0],j-results.x[1])[0,0]) for j in xrange(dim)] for i in xrange(dim)] )
    ref_aligned=np.zeros( [dim,dim] )
    for i in xrange(dim):
        for j in xrange(dim):
            ref_aligned[i,j] = results.x[2]*(interp_ref(i-results.x[0],j-results.x[1]))

    return ref_aligned



dim=64
runNumber=11

os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/F1140C/DitherStep10/Jitter0/Originals')
hdu   = pyfits.open('PSF_run'+str(runNumber)+'_ScienceTarget.fits')
target=hdu[0].data
hdu.close()

hdu   = pyfits.open('PSF_run'+str(runNumber)+'_ReferenceTarget.fits')
reference=hdu[0].data
hdu.close()


aligned_ref=align_images(runNumber,dim,reference,target)


pyfits.writeto("TEST_ALIGN_run"+str(runNumber)+".fits",aligned_ref,clobber=True)
