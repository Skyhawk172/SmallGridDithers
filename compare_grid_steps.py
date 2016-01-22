1065###########################################################################
# WHAT: AzimuthalAverage.py
#      
#      Calculates azimuthal average for every contrast image in the 
#      specified directory on command line 
#      (e.g. CDR for ~/CWG/Latency/Contrast_CDR). All the radial profiles
#      are averaged and only the average and the +/-1 sigma lines are 
#      plotted.
#
#      Also plots in inset the distribution of positions used to create
#      the images.
#
#      The contrast images are created by Mathematica Notebook 
#      JWST small-grid-dither reduction V2.0 LOOP.nb and 
#      the original PSF images live in /JWST/Simulations/.../Jitter**/.
#
#
#
#
# HOW: The call should include the directory where the contrast files live:
#      "python AzimuthalAverage.py 1550 CDR_PhotNoise"
#       
# WHO:  C-P LAJOIE
# WHEN: July 2013
###########################################################################

import numpy as np
from numpy import *
import matplotlib.pyplot as P
import pyfits, os, glob, re, sys
from matplotlib.ticker import MultipleLocator

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



#######################################
# MAIN PROGRAM
#######################################
fig ,ax  = P.subplots(figsize=( (11,8) ) )

indir=sys.argv[1:]
print indir

directory=os.getcwd()

binsize = 1
#nbins=int(dim/binsize/2 +1)
nsigmas=5
max_gain=0

xcutoff0= 1.5
xmin=0
ymin=-7.
ymax=-4.

colors = ['r', 'b','g','c','y','k','m','r']     #color code for loop
method="LOCI"



icount=0
for i in xrange(len(indir)):
    
###########################################################################
# GLOB FILES FOR ALL THE COMBINATIONS AVAILABLE IN THE CONTRAST FILES
###########################################################################
    if 'Originals' not in indir[i]: prefix='Scaled_'
    else: prefix=''

    if "F1065C" in indir[i]: 
        lambda0=10.65
        #hdu=pyfits.open(prefix+'PSF_1065_run1_Unocculted_204.fits')
    if "F1140C" in indir[i]: 
        lambda0=11.40
        #hdu=pyfits.open(prefix+'PSF_1140_run1_Unocculted_204.fits')
    if "F1550C" in indir[i]: 
        lambda0=15.50
        #hdu=pyfits.open(prefix+'PSF_1550_run1_Unocculted_204.fits')
    if "F210M"  in indir[i]: 
        lambda0=2.10
        #hdu=pyfits.open(prefix+'PSF_210_run1_Unocculted_136.fits')
    if "F430M"  in indir[i]: 
        lambda0=4.30
        #hdu=pyfits.open(prefix+'PSF_430_run1_Unocculted_136.fits')
    if "F460M"  in indir[i]: 
        lambda0=4.60
        #hdu=pyfits.open(prefix+'PSF_460_run1_Unocculted_136.fits')

    if "MIRI" in indir[i]: 
        dim=64
        pixelSize=0.11
    elif "NIRCam" in indir[i]:
        dim = 109 if lambda0>2.35 else 222
        pixelSize = 0.032 if lambda0<2.35 else 0.065

    y, x = np.indices( (dim,dim) )
    center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    nbins = int((np.round((x-center[0]).max() / binsize)+1)) #CPL: changed the max bin and added int()


    os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/'+indir[i])
    


#    unocculted=hdu[0].data
#    hdu.close()



    files=glob.glob('*Map*'+method+'_*run*.fits')
    print '\n   Reading contrast files from',indir[i]#,'\n'
    nfiles=len(files)
    print '     --> %d files\n' %nfiles
    ic=0
    rad_prof=np.empty( (nfiles,nbins) )
    for ifile in xrange(nfiles):
        hdu = pyfits.open(files[ifile])
        data=hdu[0].data
        hdu.close()
        #################################
        #NOW DONE IN MATHEMATICA NB
        #data=data/unocculted.max()
        #################################
        x,rad_prof[ic]=azimuthalAverage(data,binsize=binsize,stddev=True,interpnan=True,returnradii=True)
        ic+=1



    #get average of the stddev for each radial bin:
    av_prof =nsigmas* ( np.mean(rad_prof,axis=0) )
    std_prof=nsigmas* np.std( rad_prof,axis=0,dtype=np.float64)
    min_prof=nsigmas* np.abs( np.min( rad_prof,axis=0) )
    max_prof=nsigmas* np.abs( np.max( rad_prof,axis=0) )
    plusone =av_prof + std_prof
    minusone=av_prof - std_prof

###########################################################################
# PLOT STUFF
###########################################################################
    lambdaD=lambda0*1e-6/6.5 *180/pi*3600 #in arcseconds
    
    maxLambdaD=max(x)/(lambdaD/pixelSize)
    xLoverD = x / (lambdaD/pixelSize)
    xmax=xLoverD.max()-1.

    #for i in xrange(nfiles):
    #    ax.plot(xLoverD,np.log10(np.abs(rad_prof[i])),color=colors[icount])

    idx = (np.abs(xLoverD-xcutoff0)).argmin() +1
    xcutoff=0. #xLoverD[idx]

    tmpx  = xLoverD[xLoverD>=xcutoff]
    tmpy  = av_prof[xLoverD>=xcutoff]
    tmpp1 = plusone[xLoverD>=xcutoff]
    tmpm1 =minusone[xLoverD>=xcutoff]
    tmpmax=max_prof[xLoverD>=xcutoff]
    

    if 'DitherStep10' in indir[i]: label="Dither step 10 mas"
    if 'DitherStep20' in indir[i]: label="Dither step 20 mas"
    if 'DitherStep30' in indir[i]: label="Dither step 30 mas"
    if 'DitherStep50' in indir[i]: label="Dither step 50 mas"

#    ax.plot(xLoverD[:-1],np.log10(av_prof[:-1]),'--',label=label,color=colors[icount])
    ax.plot(tmpx[:-1],np.log10(tmpy[:-1]),'--',label=label,color=colors[icount])
    ax.grid(True)

#    ax.plot(xLoverD,np.log10(min_prof),'r')#,label='+/-1 $\sigma$')
#    ax.plot(xLoverD,np.log10(max_prof),'r')
#    ax.fill_between(xLoverD[:-1],np.log10(plusone[:-1]),np.log10(minusone[:-1]),color=colors[icount], alpha=0.25)
    ax.fill_between(tmpx[:-1],np.log10(tmpp1[:-1]),np.log10(tmpm1[:-1]),color=colors[icount], alpha=0.25)
    if icount==0: ax.add_patch(P.Rectangle( (0.,ymin), xcutoff,ymax-ymin, color='gray',alpha=0.25,hatch='.',zorder=10))
    
    icount+=1

######################
# X-AXIS ON TOP:
######################
axtop=ax.twiny()
#step=5
#xpixels=np.arange(0,max(x)+1,step)/(lambdaD/pixelSize)
#axtop.set_xticks(xpixels)
#axtop.set_xticklabels(np.arange(0,max(x)+1,step,dtype=int))
#axtop.set_xlabel('MIRI pixels')


step=0.5/pixelSize
xpixels=np.arange(0,max(x)+1,step)*pixelSize/lambdaD
axtop.set_xticks(xpixels)
axtop.set_xticklabels(xpixels*lambdaD)
axtop.set_xlabel('Arcsec')
axtop.set_xlim(0.,xmax)
               


ax.set_xticks(np.arange(0,maxLambdaD))

ax.set_ylabel('%d$\sigma$ contrast (log)' % nsigmas)
#ax.set_ylim(-6.,-3.5)#ymin,ymax) #np.min(np.log10(av_prof))-1,np.max(np.log10(av_prof))+0.5)
ax.set_ylim(-6.5,-4.5)#ymin,ymax) #np.min(np.log10(av_prof))-1,np.max(np.log10(av_prof))+0.5)
ax.set_xlim(0.,xmax)
ax.set_xticks(np.arange(0., xmax, 1.0))
ax.legend(loc=1)#bbox_to_anchor = (0.35, 0.27))#loc=3)
note='             '+method+"\nBased on %d pointings" %nfiles
P.text(0.75*xmax,ymax-0.5*np.abs(ymax-ymin),note)

ax.set_xlabel(r"$\lambda/D$ (at %5.2f $\mu$m)" %lambda0)






filename=directory+"/Dithered_10vs30mas"+str(lambda0)+".pdf"
raw_input("Press enter to save PDF")
print '   --> Saving to file',filename
P.savefig(filename)







