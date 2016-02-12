##############################################################
# WHAT: generate PSF with WebbPSF on a small grid based on 
#       various command line arguments provided by the users.
#       Used for small-grid dither simulations. 
#       Can be used for MIRI or NIRCAM.
#       
# output: fits files for all the psf generated as well as 
#         txt files (runxx.txt) for the positions of the
#         sources on the grid. 
#
#         These PSFs can be used by AddNoisetoPSF.py to
#         add noise to them, or they can be used directly 
#         in SGD-LOCI code to calculate PSF subtraction.
#
# HOW: python generate_PSFs.py --help
#      python generate_PSFs.py -I MIRI   -f f1065c -mask FQPM1065 -stop maskfqpm -gridstep 20 -jitter 7 -rms 204 -nruns 10
#      python generate_PSFs.py -I Nircam -f f210m  -mask mask210r -stop circlyot -gridstep 20 -jitter 7 -rms 136 -nruns 10
#
# WHO: Charles-Philippe Lajoie
#
# WHEN: April 2014
##############################################################

import webbpsf
import numpy as np
import pysynphot
import pyfits
import astropy
import math
import argparse, sys, os, logging, pkg_resources



def set_sim_params(args):
    #logging.basicConfig(level=logging.INFO, format='%(name)-12s: %(levelname)-8s %(message)s',)
    print "WebbPSF",pkg_resources.get_distribution("webbpsf").version," Poppy",pkg_resources.get_distribution("poppy").version
    #webbpsf.setup_logging(level='ERROR')

    nruns= args.nruns            #
    filt = args.f.upper()        #e.g. F430M
    mask_coron=args.mask.upper() #e.g. MASK430R
    pupil_stop=args.stop.upper() #e.g. CIRCLYOT
    instr=args.I.upper()

    if instr=='NIRCAM': lambda0=float(filt[1:4])/100
    elif instr=='MIRI': lambda0=float(filt[1:5])/100
    jitter= args.jitter 
    grid_step=args.gstep
    gridpoints=args.gnpts
    side=int(np.sqrt(gridpoints))
    max_step=np.floor(side/2)
    rms= args.rms 

    fovarcsec = args.fov #default 7.04 
    #fovpixels = 64

    sigmaTA = 4.7
    sigmaFSM= 2.0


    if instr=='NIRCAM': img=webbpsf.NIRCam()
    elif instr=='MIRI': img=webbpsf.MIRI()

    if filt not in img.filter_list:
        if instr == 'MIRI'  : print "\nFilter "+mask_coron+" not available. Select from: ",[f for f in img.filter_list if "C" in f[-1]],'\n'
        if instr == 'NIRCAM': print "\nFilter "+mask_coron+" not available. Select from: ",img.filter_list,'\n'
        sys.exit()

    if mask_coron not in img.image_mask_list:
        print "\nCoronagraphic mask "+mask_coron+" not available. Select from: ",img.image_mask_list,'\n'
        sys.exit()

    if pupil_stop not in img.pupil_mask_list[:2]:
        print "\nPupil mask "+pupil_stop+" not available. Select from: ",img.pupil_mask_list[:2],'\n' 
        sys.exit()

    opd_rms = [int(img.opd_list[i][-8:-5]) for i in xrange(len(img.opd_list))]
    if rms not in opd_rms:
        print "\nOPD rms ("+str(rms)+") not available. Select from: ",opd_rms,'\n'
        sys.exit()
    if instr=='NIRCAM' and args.noopd==False: opd='/Users/lajoie/WebbPSF/webbpsf-data/NIRCam/OPD/OPD_RevV_nircam_'+str(rms)+'.fits'
    if instr=='MIRI'   and args.noopd==False: opd='/Users/lajoie/WebbPSF/webbpsf-data/MIRI/OPD/OPD_RevV_miri_'+str(rms)+'.fits'

    img.filter=filt
    img._rotation=0.
    img.options["output_mode"]='detector sampled'
    if args.noopd==False: 
        img.pupilopd = (opd, 0) #select FITS extension for OPD
        print opd

    if instr=='NIRCAM': 
        outdir='./Results'+'/'+instr+'/'+filt+'_'+mask_coron+'/DitherStep'+str(int(grid_step))+'/Jitter'+str(int(jitter))+'/rms'+str(rms)+'nm/'
    else: outdir='./Results'+'/'+instr+'/'+filt+'/DitherStep'+str(int(grid_step))+'/Jitter'+str(int(jitter))+'/rms'+str(rms)+'nm/'

    if gridpoints != 9: 
        raw_input(" *** generating grids of %i points: ENTER to continue " %gridpoints)
        outdir = outdir[:-1]+"_Grid_%i_pts/" %gridpoints
        

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    x0,y0 = get_source_offset(lambda0, filt, mask_coron)
    if np.abs(x0) > fovarcsec/2: 
        print "\n optimal bar occulter offset out of FOV. Please fix FOV in code (fovarcsec).\n   x_offset=%f Total fov=%f\n" %(x0,fovarcsec)
        sys.exit()

    generate_PSF(img, outdir, nruns, x0, y0, grid_step, side, max_step, sigmaTA, sigmaFSM, mask_coron, pupil_stop, lambda0, jitter, rms, fovarcsec)

    return

    
def get_source_offset(lambda0, filt, mask_coron):
    ########################################################
    # Find optimal position for source behind bar occulters.
    # See document from John Stansberry (email 04/08/15) 
    # for equations below; modified for orientation of bars.
    #
    # Position (x,y)=(0,0) for both SWB and LWB are 4 lambda/D 
    # wide for 2.10 and 4.60 um respectively. Use scaling below  
    # to offset for other wavelengths. Results are arcseconds.
    #
    # Bars are 20" wide with flat edges. Linear part is 15" wide.
    # 
    # Width of bar occulters as a function of x, where x=0
    # is the center of the bars:
    # WSWB(x) (") = 0.2666" + 0.01777 * x
    # WLWB(x) (") = 0.5839" + 0.03893 * x
    ########################################################

    y0 = 0
    if mask_coron == 'MASKSWB':
        if filt[-1]=="W": x0 =-6.83*( np.sqrt(lambda0) - 2.196 )
        else:             x0 =-7.14*( lambda0 - 2.10 )
    elif mask_coron == 'MASKLWB':
        if filt[-1]=="W": x0 =-3.16*( np.sqrt(lambda0) - 4.747 )
        else:             x0 =-3.26*( lambda0 - 4.60 )
    else: x0,y0 = 0.,0.

    return x0,y0


def make_ref_grid(n, outdir, x0, y0, side, max_step, grid_step, sigmaTA, sigmaFSM):
    #####################################
    # GENERATE GRIDS + TARGET/REFERENCE
    # units are milli-arcseconds. 
    #####################################
    grid=[]
    Target_pos=([np.random.normal(x0*1000, sigmaTA),np.random.normal(y0*1000, sigmaTA)])
    Ref_pos=([np.random.normal(x0*1000, sigmaTA),np.random.normal(y0*1000, sigmaTA)])
    grid=[Target_pos,Ref_pos]

    # below: generating the pointings using Remi's method from Mathematica notebook.
    #for i in xrange(-1,2,1):
    #    for j in xrange(-1,2,1):
    #        grid.append( [np.random.normal(Ref_pos[0]+i*grid_step,sigmaFSM),np.random.normal(Ref_pos[1]+j*grid_step,sigmaFSM)] )


    # below: generating the pointings by walking from one poinitng to another, unlike just above where the pointings are 
    # generated from the reference target position. 
    grid.append( [np.random.normal(grid[-1][0]-max_step*grid_step,sigmaFSM),np.random.normal(grid[-1][1]-max_step*grid_step,sigmaFSM)] )
    yjump=0
    for j in xrange(side):
        for i in xrange(side):
            if i!=0 or j!=0:
                if i!=0:
                    if j%2 == 0: x=np.random.normal( grid[-1][1]+grid_step,sigmaFSM )
                    else       : x=np.random.normal( grid[-1][1]-grid_step,sigmaFSM )
                else: x=np.random.normal( grid[-1][1],sigmaFSM )

                y=np.random.normal(grid[-1][0]+yjump*grid_step,sigmaFSM) 
                grid.append( [y,x] )
            yjump=0
        yjump=1
              
    grid=np.array(grid) 
    outdir_runs=outdir.split("Jitter")[0]
    np.savetxt(outdir_runs+"run"+str(n)+".txt",grid.flatten())

    return grid


def generate_PSF(img, outdir, nruns, x0, y0, grid_step, side, max_step, sigmaTA, sigmaFSM, mask_coron, pupil_stop, lambda0, jitter, rms, fovarcsec):
    src=pysynphot.Icat("ck04models",5800, 0.0, 4.44) #temp, z, logg
    print src

    njitter=50 #number of PSF to use to create the final jittered PSF
    for n in xrange(1,nruns+1):
        grid = make_ref_grid(n, outdir, x0, y0, side, max_step, grid_step, sigmaTA, sigmaFSM)
        
        #########################################
        # GENERATE PSF FOR EACH POSITION IN GRID
        #########################################
        for isgd in xrange(len(grid)+1):
            if   isgd==0 : label='Unocculted'
            elif isgd==1 : label='ScienceTarget'
            elif isgd==2 : label='ReferenceTarget'
            else:          label='Reference_dither%i' %(isgd-2)
            output_name='PSF_%i_run%i_%s_%i.fits' %(np.round(lambda0*100),n,label,rms)

            if isgd==0: img.image_mask=None
            else:       img.image_mask=mask_coron
            img.pupil_mask=pupil_stop

            if isgd == 0 or isgd == 1: currentPointing=grid[0]
            else: currentPointing = grid[isgd-1]

            offset= np.sqrt(np.sum( currentPointing**2 ) )/1000 
            offset_theta= -np.arctan2(currentPointing[0],currentPointing[1])*180/math.pi
            img.options['source_offset_r']     = offset         # offset in arcseconds
            img.options['source_offset_theta'] = offset_theta   # degrees counterclockwise from instrumental +Y in the science frame

            print isgd,output_name, currentPointing, offset,offset_theta,"\n"

            psf = img.calcPSF(source=src, fft_oversample=4, detector_oversample= 1, fov_arcsec=fovarcsec ) #fov_pixels=fovpixels)
            finalPSF = psf[0].data
            hdr = psf[0].header


            string="LAJOIE: optimal source offset for bar occulter: %f arcseconds" %x0
            hdr.add_history(string)
            string="LAJOIE: actual radial offset is %f arcseconds" %offset
            hdr.add_history(string)
            string="LAJOIE: actual angular offset is %f degrees" %offset_theta
            hdr.add_history(string)
            if img.image_mask==None : hdr.set("CORONMSK",img.image_mask,before='PUPIL') 

            if jitter!=0:
                sumPSF=finalPSF
                for j in xrange(njitter):
                    jitter_pos=np.array([np.random.normal(currentPointing[0],jitter),np.random.normal(currentPointing[1],jitter)])

                    off_jitt       = np.sqrt(np.sum( jitter_pos**2 ) )/1000              # arcseconds
                    off_jitt_theta =-np.arctan2(jitter_pos[0],jitter_pos[1])*180/math.pi # degrees
                    img.options['source_offset_r']     = off_jitt                        # offset in arcseconds
                    img.options['source_offset_theta'] = off_jitt_theta                  # degrees counterclockwise from instrumental +Y in the science frame

                    tmp=img.calcPSF(source=src, fft_oversample=4, detector_oversample= 1, fov_arcsec=fovarcsec) #fov_pixels=fovpixels) 
                    sumPSF+=tmp[0].data

                finalPSF=sumPSF/(njitter+1)
                string="LAJOIE: Jitter (%4.1f mas 1-sigma/axis) included, %d points drawn from Gaussian dist." %(jitter,njitter)
                hdr.add_history(string)

            else: hdr.add_history("LAJOIE: No jitter included")

            ##########################
            # WRITE OUTPUT IMAGE
            ##########################
            outHDU=astropy.io.fits.PrimaryHDU(data=finalPSF, header=hdr )
            outHDU.writeto(outdir+output_name)


    return





#####################################################################
# MAIN Method
#####################################################################
if __name__ == "__main__":
    np.random.seed(0)

    parser = argparse.ArgumentParser(description="Generate PSF grids for LOCI using WebbPSF")

    parser.add_argument("-I"      , type=str, metavar="INSTRUMENT",required=True, help="Instrument to use (MIRI or NIRCam)")
    parser.add_argument("-f"      , type=str, metavar="FILTER"    ,required=True, help="Filter to use")
    parser.add_argument("-mask"   , type=str, metavar="MASK_CORON",required=True, help="Mask coron. to use")
    parser.add_argument("-stop"   , type=str, metavar="PUPIL_STOP",required=True, help="Pupil stop to use")
    parser.add_argument("--rms"   , type=int  , default=136, help="Select OPD rms to use, if available (default: 136 nm)")
    parser.add_argument("--noopd" , action="store_true",     help="Do not use any OPD (default: False)")
    parser.add_argument("--jitter", type=float, default=0  , help="Sigma Jitter (default: 0 mas)")
    parser.add_argument("--fov"   , type=float, default=7.04,help="Field of view (diam.; default=7.04 arcsecond)")
    parser.add_argument("--nruns" , type=int  , default=1  , help="Number of SGD grids to generate (default: 1)")
    parser.add_argument("--gstep" , type=float, default=20 , help="SGD grid steps (default: 20 mas)")
    parser.add_argument("--gnpts" , type=float, default=9  , help="SGD square grid points (default: 9)")

    args = parser.parse_args()

    set_sim_params(args)




        
