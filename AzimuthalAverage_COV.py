##########################################################################
# WHAT: AzimuthalAverage_COV.py
#      
#       This is a modified and lighter version of AzimuthalAverage_SGD.py
#
#       Calculates the noise matrix from all the contrast maps in the 
#       input folder according to the JWST ETC (Kyle van Gorkom). The 
#       contrast is then plotted for all of the different methods 
#       (e.g. Classical, LOCI, LOCI5pt, etc) available.
#
# HOW: python AzimuthalAverage_COV.py --help
#       
# WHO:  C-P LAJOIE
# WHEN: April 2016
###########################################################################

import numpy as np
import matplotlib.pyplot as P
import argparse, pyfits
import os, glob, sys



# FROM KYLE:
def get_aperture(aperture_image):
    ''' Get the aperture at each pixel position. Used in calculating noise and radial profile
    (for the contrast normalization).

    Parameters:
        aperture_image : nd array
            An nxn image of the aperture centered at the central pixel. The dimensions
            of this image should match

    Returns:
        ap_array : np array
            An n^2xn^2 array representing the aperture centered at each pixel.
    '''

    kernel = aperture_image
    dim = kernel.shape

    #Embed the kernel in a large array
    desired = (dim[0]*2 + 1)/2
    add = ((desired - dim[0]/2,desired - dim[0]/2),(desired - dim[1]/2,desired - dim[1]/2))

    kernel_large = np.pad(kernel,add,mode='constant', constant_values=(0))

    #Make the array that represents the flattened aperture centered at each pixel
    ap_array = []
    size = dim[0]
    for x in np.arange(size):
        for v in np.arange(size):
            ap_array.append(kernel_large[size-x:2*size-x,
                size-v:2*size-v].flatten())

    return np.array(ap_array)



def azimuthalAverage(image, center=None, stddev=False, returnradii=False, return_nr=False, 
        binsize=0.5, weights=None, steps=False, interpnan=False, left=None, right=None):
    """
    see  https://code.google.com/p/agpy/source/browse/trunk/agpy/radialprofile.py?r=317

    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
    stddev - if specified, return the azimuthal standard deviation instead of the average
    returnradii - if specified, return (radii_array,radial_profile)
    return_nr   - if specified, return number of pixels per radius *and* radius
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    weights - can do a weighted average instead of a simple average if this keyword parameter
        is set.  weights.shape must = image.shape.  weighted stddev is undefined, so don't
        set weights and stddev.
    steps - if specified, will return a double-length bin array and radial
        profile so you can plot a step-form radial profile (which more accurately
        represents what's going on)
    interpnan - Interpolate over NAN values, i.e. bins where there is no data?
        left,right - passed to interpnan; they set the extrapolated values

    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    nbins = int((np.round((x-center[0]).max() / binsize)+1)) #CPL: changed the max bin and added int()
    maxbin = nbins * binsize
    bins = np.linspace(0,maxbin,nbins+1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:]+bins[:-1])/2.0

    # Find out which radial bin each point in the map belongs to
    whichbin = np.digitize(r.flat,bins)

    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    nr = np.bincount(whichbin)[1:]

    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or xrange(1,nbins+1) )
    # radial_prof.shape = bin_centers.shape
    if stddev:
        radial_prof = np.array([image.flat[whichbin==b].std() for b in xrange(1,nbins+1)])
    else:
        radial_prof = np.array([(image*weights).flat[whichbin==b].sum() / weights.flat[whichbin==b].sum() for b in xrange(1,nbins+1)])

    #import pdb; pdb.set_trace()

    if interpnan:
        radial_prof = np.interp(bin_centers,bin_centers[radial_prof==radial_prof],radial_prof[radial_prof==radial_prof],left=left,right=right)

    if steps:
        xarr = np.array(zip(bins[:-1],bins[1:])).ravel() 
        yarr = np.array(zip(radial_prof,radial_prof)).ravel() 
        return xarr,yarr
    elif returnradii: 
        return bin_centers,radial_prof
    elif return_nr:
        return nr,bin_centers,radial_prof
    else:
        return radial_prof


def plot_contrast(instr, filt, lambda0, x, av_prof, nsig, kw):
    print instr
    fig ,(ax,ax2)  = P.subplots(2,figsize=( (11,8) ) )

    # size of input images (7.04 arcseconds)
    if   instr=="MIRI":   pixelSize = 0.11
    elif instr=="NIRCam": pixelSize = 0.032 if lambda0<=2.35 else 0.065

    JWSTdiam = 6.61099137008
    lambdaD=lambda0*1e-6/JWSTdiam *180/np.pi*3600.  # in arcseconds
    xLoverD = x / (lambdaD/pixelSize)
    xmax=np.max(xLoverD)

    # PLOT PROFILES AND RELATIVE GAIN:
    arcseconds = x*pixelSize
    for i in xrange(len(kw)):
        ax.plot(x*pixelSize, np.log10(nsig*av_prof[i]),'-', label=kw[i])
        if i>0: ax2.plot(arcseconds, av_prof[0]/av_prof[i])
        print kw[i],"\n",np.log10(av_prof[i])


    # PLOT AXIS + LEGEND:
    ax.legend(loc=1)
    ax.set_xlabel("Separation (arcsec)")
    ax.set_ylabel('%d$\sigma$ contrast (Log)' % nsig)
    ax.set_ylim(-8., -2.) 
    ax.set_xlim(0., np.max(arcseconds) )

    ax2.text(0.32,-0.25, 'Gain over Classical Subtraction',horizontalalignment='left',
             verticalalignment='top',transform=ax.transAxes,fontsize=13, family='serif',fontweight='light') 
    ax2.set_xlim(0.,np.max(arcseconds))
    ax2.set_ylabel('Gain in contrast')
    ax2.set_xlabel("Separation (arcsec)")

    ax.grid(True)
    ax2.grid(True)

    axtop = ax.twiny()
    axtop.set_xticks(np.arange(0., xmax, int(np.ceil(xmax)/10.)))    
    axtop.set_xlabel('Arcsec')
    axtop.set_xlim(0.,xmax)
    axtop.set_xlabel(r"$\lambda/D$ (at %5.2f $\mu$m)" %lambda0)

    # SAVE FIGURE:
    filename="AzimuthalAverage_"+filt+".pdf"
    raw_input("\nPress enter to save PDF to '%s'" %filename)
    print '   --> Saving to file',filename
    P.savefig(filename)



    return 





#######################################
# MAIN PROGRAM
#######################################
if __name__ == "__main__":


    #############################################
    #INPUT PARAMETERS AND GLOB FILES:
    #############################################
    parser = argparse.ArgumentParser(description="Calculate noise and plot radial contrast profiles")
    parser.add_argument("dir" , type=str, help="Input directory")
    parser.add_argument("-s"  , type=float, default= 5., help="Nsigmas to plot (default=5)")
    args = parser.parse_args()


    indir= args.dir
    filt = indir.split('/')[1]
    nsigmas = args.s
    

    if "Planets" in indir: 
        print "\n These maps have planets, continue? (press enter)"
        raw_input()
    if "F210M"  in indir: lambda0= 2.10
    if "F430M"  in indir: lambda0= 4.30
    if "F460M"  in indir: lambda0= 4.60
    if "F1065C" in indir: lambda0=10.65
    if "F1140C" in indir: lambda0=11.40
    if "F1550C" in indir: lambda0=15.50

    # size of input images (7.04 arcseconds)
    if "MIRI" in indir:     instr = 'MIRI'
    elif "NIRCam" in indir: instr = 'NIRCAM'
        

    os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/Dither-LOCI/Results/'+indir)
    #os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/CWG/SGD/'+indir)
    directory=os.getcwd()


    input_Unocculted=glob.glob('*run1_Unocculted*.fits')
    hdu=pyfits.open(input_Unocculted[0])
    unocculted=hdu[0].data
    hdu.close()
    dim = unocculted.shape[0]
    npix = dim**2

    # Binary aperture centered on (xmid, ymid): 0 within "radius", 1 outside
    xmid, ymid= dim/2-1, dim/2-1
    radius = 4
    aperture = np.array([ [1. if (np.sqrt( (i-xmid)**2 + (j-ymid)**2) < radius) else 0. for i in xrange(dim)] for j in xrange(dim)])
    psf_aper = np.convolve( unocculted.flatten(), aperture.flatten() )

    binsize= 1.0 
    y, x = np.indices( (dim,dim) )
    center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    nbins = int((np.round((x-center[0]).max() / binsize)+1)) #CPL: changed the max bin and added int()

    keywords=['CLAS','LOCI','LOCI5pt','LOCI25pts']#,'BOCC'] #different maps already generated
    if 'LOCI' in sys.argv or 'CLAS' in sys.argv or 'BOCC' in sys.argv: keywords=[sys.argv[-1]]



    #############################################
    # CALCULATE NOISE AND AZIMUTHAL PROFILES:
    #############################################
    print '\n   Reading contrast files from %s\n' %indir

    ic = 0
    types=[]
    rad_prof = np.empty( (len(keywords), nbins) )
    for kw in keywords:

        files=glob.glob('*Map*'+kw+'_run*.fits')
        nfiles=len(files)
        
        if nfiles==0: 
            print '      Error: no %s contrast maps found\n' %kw

        else: 
            types.append(kw)
            print '      Processing %s files #:'%kw,

            images = np.empty( (nfiles, npix) )
            for i in xrange(nfiles):
                image=files[i]
                print i+1,
                hdu = pyfits.open(image)
                images[i] =  hdu[0].data.flatten()
                hdu.close()

            images        = images.T
            covmatrix     = images.dot(images.T)/(nfiles-1)
            A             = get_aperture(aperture)
            noise_matrix  = A.dot(covmatrix.dot(A.T))
            noise         = np.sqrt( np.diag(noise_matrix) ).reshape( (dim,dim) )
            x,rad_prof[ic]= azimuthalAverage(noise, binsize=binsize, stddev=False, interpnan=True, returnradii=True)
            rad_prof[ic]  = rad_prof[ic]/psf_aper.max()

            ic+=1
            print "\n"


    #############################################
    # PLOT CONTRAST CURVES:
    #############################################
    plot_contrast(instr, filt, lambda0, x, rad_prof, nsigmas, types)








