##########################################################################
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


def plot_TApositions():
    ###########################################################################
    # PARSE POSITIONS OUT OF STRING AND CALCULATE MEAN OFFSET AND STD DEV.
    # AND PLOT INSET
    ###########################################################################
    os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/'+indir[:19])
    print "    Reading position files from ",indir[:19]
    pos_files=glob.glob('run*.txt')

    xpos=[]
    ypos=[]
    for i in xrange(len(pos_files)):
        pos =np.loadtxt(pos_files[i],unpack=True)
        pos2=np.reshape(pos,(11,2))
        xpos.append(pos2[0,0])
        ypos.append(pos2[0,1])
        xpos.append(pos2[1,0])
        ypos.append(pos2[1,1])

    meanx=np.average(xpos)
    meany=np.average(ypos)
    radius=np.sqrt( (xpos-meanx)**2 + (ypos-meany)**2)
    stdr=np.sqrt(sum(radius*radius)/len(radius))#std[j]=np.std(radius)

    inset=P.axes([0.62,0.57,0.3,0.3],aspect='equal') #This defines the inset by its origin and size! and this is the only thing you need to do!
    P.plot(xpos,ypos,'b.')
    P.axhline(color='k',linestyle="--")
    P.axvline(color='k',linestyle="--")
    #P.title('FQPM'+sys.argv[1]+'\n'+sys.argv[2])
    #P.xticks(np.arange(-50,50+1,10))
    #P.yticks(np.arange(-50,50+1,10))
    theta=np.arange(0,2*pi,0.05)
    xx=meanx+stdr*cos(theta)
    yy=meany+stdr*sin(theta)
    P.plot(xx,yy,color='r')

    xx=meanx+2*stdr*cos(theta)
    yy=meany+2*stdr*sin(theta)
    P.plot(xx,yy,color='r')
    text="$\sigma$= "+str("%.3f"%stdr)
    P.text(-45,40,text)
    offset=np.sqrt(meanx*meanx + meany*meany)
    text="Offset= "+str("%.3f"%offset)
    P.text(-45,30,text)

    P.xlim(-50,50)
    P.ylim(-50,50)
    P.xlabel('x mas')
    P.ylabel('y mas')
    #P.setp(inset,xticks=[0.0005,0.0015,0.0025])

    return


def plot_labels(max_gain,indir):
    #ax2.add_patch(P.Rectangle( (0.,0), xcutoff,max_gain, color='gray',alpha=0.25,hatch='.',zorder=10))
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
               

    split=re.split('/',sys.argv[1])
    try: 
        note='Dither step='+re.split(r'(\d+)',split[2])[1]+'   Jitter='+re.split(r'(\d+)',split[3])[1]
        P.text(xcutoff+35.5,ymax-0.5*np.abs(ymax-ymin),note)
    except: pass

    ax.set_xticks(np.arange(0,maxLambdaD))
    ax.set_ylabel('%d$\sigma$ contrast (log)' % nsigmas)
    ax.set_ylim(ymin,ymax) #np.min(np.log10(av_prof))-1,np.max(np.log10(av_prof))+0.5)
    ax.set_xlim(0.,xmax)
    ax.set_xticks(np.arange(0., xmax, int(np.ceil(xmax)/10.)))

    note="Based on %d pointings" %nfiles
    P.text(xcutoff+35.5,ymax-0.6*np.abs(ymax-ymin),note)

    if len(sys.argv)<3:
        ax.legend(loc=1)#bbox_to_anchor = (0.35, 0.27))#loc=3)

        ax2.legend(loc=1)#bbox_to_anchor = (0.35, 0.27))#loc=3)
        ax2.set_xticks(np.arange(0., xmax+1, int(np.ceil(xmax)/10.)))
        ax2.set_ylim(0.,max_gain)
        ax2.set_xlim(0.,xmax)
        ax2.set_ylabel('Gain in contrast')
        ax2.set_xlabel(r"$\lambda/D$ (at %5.2f $\mu$m)" %lambda0)
        ax2.grid(True)
    else: ax.set_xlabel(r"$\lambda/D$ (at %5.2f $\mu$m)" %lambda0)

    ax2.set_title(indir)
        
    return 





def plot_contrast_curves(x,av_prof,xLoverD,plusone,minusone,max_prof,max_gain):

    tmpx  = xLoverD[xLoverD>=xcutoff]
    tmpy  = av_prof[xLoverD>=xcutoff]
    tmpp1 = plusone[xLoverD>=xcutoff]
    tmpm1 = minusone[xLoverD>=xcutoff]
    tmpmax= max_prof[xLoverD>=xcutoff]
    
    if kw=='LOCI'     : label='LOCI 9 pts'
    if kw=='LOCI5pt'  : label='LOCI 5 pts'
    if kw=='LOCI25pts': label='LOCI 25 pts'
    if kw=='CLAS'     : label='Classical subtraction'
    if kw=='BOCC'     : label='Classical sub. (Boccaletti @ 5mas)'

    ax.plot(tmpx[:-1],np.log10(tmpy[:-1]),'-',label=label,color=colors[icount])
    ax.grid(True)

    print xLoverD
    print av_prof

    #ax.plot(xLoverD,np.log10(min_prof),'r')#,label='+/-1 $\sigma$')
    if len(sys.argv)>2: ax.plot(tmpx[:-1],np.log10(tmpmax[:-1]),'r')
    ax.fill_between(tmpx[:-1],np.log10(tmpp1[:-1]),np.log10(tmpm1[:-1]),color='gray', alpha=0.25)
    if icount==1: ax.add_patch(P.Rectangle( (0.,ymin), xcutoff,ymax-ymin, color='gray',alpha=0.25,hatch='.',zorder=10))

    if icount>0: 
        gain=av_prof_CLAS/av_prof
        gain[gain==inf]=0.
        tmpgain=gain[xLoverD>=xcutoff]
        if tmpgain[:-1].max()>max_gain: max_gain=tmpgain[:-1].max()
        if kw!='CLAS': 
            ax2.text(0.32,-0.25, 'Gain over Classical Subtraction',horizontalalignment='left',
                    verticalalignment='top',transform=ax.transAxes,fontsize=13, family='serif',fontweight='light') 
            ax2.plot(tmpx[:-1],tmpgain[:-1],'-',color=colors[icount],label=label)

        #PLOT INDIVIDUAL CURVES
#        if kw=='LOCI':
#            for i in xrange(nfiles):
#                ax.plot(xLoverD,np.log10(np.abs(rad_prof[i])),color='gray')

    return max_gain









#######################################
# MAIN PROGRAM
#######################################
if len(sys.argv)<3: fig ,(ax,ax2)  = P.subplots(2,figsize=( (11,8) ) )
else: fig ,ax  = P.subplots(figsize=( (11,8) ) )

indir=sys.argv[1]
filt = indir.split('/')[1]

if "Planets" in indir: 
    print "\nThese maps have planets, continue? (press enter)\n"
    raw_input()
if "F210M" in indir : lambda0=2.10
if "F212N" in indir : lambda0=2.12
if "F430M" in indir : lambda0=4.30
if "F460M" in indir : lambda0=4.60
if "F1065C" in indir: lambda0=10.65
if "F1140C" in indir: lambda0=11.40
if "F1550C" in indir: lambda0=15.50


# size of input images (7.04 arcseconds)
if "MIRI" in sys.argv[1]: 
    instr     = 'MIRI'
    pixelSize = 0.11
    dim       = 64
elif "NIRCam" in sys.argv[1]: 
    instr    = 'NIRCAM'
    pixelSize= 0.032 if lambda0<=2.35 else 0.065
    dim      =   222 if lambda0<=2.35 else 109

    
binsize=1.1 #2.5
y, x = np.indices( (dim,dim) )
center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
nbins = int((np.round((x-center[0]).max() / binsize)+1)) #CPL: changed the max bin and added int()

nsigmas=5

xcutoff0=1.5 #where to start drawing lines and add the grey rectangle
xmin=0
ymin=-8.
ymax= -2.
keywords=['CLAS','LOCI5pt','LOCI','LOCI25pts']#,'BOCC'] #different maps already generated
colors = [  'b',   'g',   'r',  'm', 'y']     #color code for loop




#############################################
# GLOB FILES FOR ALL THE COMBINATIONS 
# AVAILABLE IN THE CONTRAST FILES
#############################################
os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/'+indir)
#os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/CWG/SGD/'+indir)
directory=os.getcwd()

if 'Originals' not in indir: prefix='Scaled_'
else: prefix=''


#input_Unocculted=    files=glob.glob('*run1_Unocculted*.fits')
#hdu=pyfits.open(input_Unocculted[0])
#unocculted=hdu[0].data
#hdu.close()


JWSTdiam = 6.61099137008
lambdaD=lambda0*1e-6/JWSTdiam *180/pi*3600 #in arcseconds


if 'LOCI' in sys.argv or 'CLAS' in sys.argv or 'BOCC' in sys.argv: keywords=[sys.argv[-1]]

icount=0
for kw in reversed(keywords):
    files=glob.glob('*Map*'+kw+'_run*.fits')
    if len(files)==0: 
        print ' Error: no Contrast maps fits files found for keyword %s\n' %kw
        keywords.remove(kw)

if len(keywords)==0: sys.exit()


###############################
# CALCULATE AZIMUTHAL PROFILES
################################
max_gain=0
for kw in keywords:
    files=glob.glob('*Map*'+kw+'_run*.fits')
    print '\n   Reading contrast files from',directory[93:],'\n'
    nfiles=len(files)

    print '      Processing %s files #:'%kw,
    rad_prof=np.empty( (nfiles,nbins) )
    for i in xrange(nfiles):
        image=files[i]
        print i+1,
        hdu = pyfits.open(image)
        data=hdu[0].data
        hdu.close()
        ##################################################
        # NOW DONE IN MATHEMATICA REDUCTION NOTEBOOK
        #data=data/unocculted.max()
        ##################################################
        x,rad_prof[i]=azimuthalAverage(data,binsize=binsize,stddev=True,interpnan=True,returnradii=True)

        #x,rad_prof[i]=azimuthalAverage(data,binsize=binsize,stddev=False,interpnan=True,returnradii=True)

    print '\n'

    #get average of the stddev for each radial bin:
    #av_prof = nsigmas* ( np.mean(rad_prof,axis=0) )
    av_prof = nsigmas* ( np.mean(rad_prof,axis=0) )
    std_prof= nsigmas* np.std( rad_prof,axis=0,dtype=np.float64)
    min_prof= nsigmas* np.abs( np.min( rad_prof,axis=0) )
    max_prof= nsigmas* np.abs( np.max( rad_prof,axis=0) )

    if kw=='CLAS': av_prof_CLAS=av_prof #save classical subtraction profile for normalization

    plusone =av_prof + std_prof
    minusone=av_prof - std_prof
    #minusone[minusone<=0]=1e-12
    #plusone[plusone<=0]=1e-12



    ###########################################
    # PLOT CURVES
    ###########################################

    maxLambdaD=max(x)/(lambdaD/pixelSize)
    xLoverD = x / (lambdaD/pixelSize)
    xmax=xLoverD.max()-1.
    idx = (np.abs(xLoverD-xcutoff0)).argmin() +1

    xcutoff= 0 #xLoverD[idx]

    max_gain=plot_contrast_curves(x,av_prof,xLoverD,plusone,minusone,max_prof,max_gain)
    icount+=1



########################
# PLOT & SAVE STUFF
########################
plot_labels(max_gain,indir)

if len(sys.argv)<3:
    #filename="AzimuthalAverage_F"+str(int(np.round(100*lambda0)))+"C.pdf"
    filename="AzimuthalAverage_"+filt+".pdf"
else: 
    #filename="AzimuthalAverage_F"+str(int(np.round(100*lambda0)))+"C_"+kw+".pdf"
    filename="AzimuthalAverage_"+filt+"_"+kw+".pdf"
   


if len(sys.argv)>2: plot_TApositions()

os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/'+indir)
#os.chdir('/Users/lajoie/Documents/Work/Projects/JWST/CWG/SGD/'+indir)
print indir
raw_input("Press enter to save PDF to '%s'" %filename)
print '   --> Saving to file',filename
P.savefig(filename)




