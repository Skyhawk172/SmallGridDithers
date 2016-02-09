##################################################################
#
# WHAT: Program that takes normalized coronagraphic PSF and 
# scales them up using pysynphot according to specified
# Teff, logg, z.
#
# The model is based on the JWST Simulator (see Brian
# York) and uses the same unmodified modules.
#
# The scaled PSF are then to be used with mathematica
# notebook JWST_small_grid_dither_reduction_v2.0LOOP.nb
# in order to apply the LOCI algorithm.
# 
# The resulting coronagraphic images can then be averaged
# azimuthally to draw contrast maps using AzimuthalAverage.py
#
# HOW: python try_pysynphotV2.py --help
#
# WHO : C-P Lajoie STScI
# WHEN: March/April 2014
#
##################################################################

import numpy as np
import pysynphot as S
import csv, sys, glob, os, pyfits, astropy
import pylab as P
import JwstImage, ObservationModule, DefaultSettings, MakeCosmicRay
import argparse

#Saturation level for detectors
SATUR_LVL={'MIRI':200000,'NIRCAM':200000} #in photons

#BACKGROUNG levels include QE and Germanium transmission:  photon/sec/pixel. see CWG/ColorPSF/zodiacal.py #
MY_BGLEVEL={
    'MIRI':{'F1065C':2.0743E1, 'F1140C':2.4177E1, 'F1550C':1.244E1,
            'F2300C':3.1405E3,'FND':7.36E-1}
    #,'NIRCam':{'F212N':1.0}
}


class Spectrum(object):
    def __init__(self, instr, filt):
        self.instr= instr
        self.filt = filt

    def getFiltername(self):
        return self.filt

    def getSpectrumInfo(self):
        return ("Teff= %.2f, z=%.2f, log_g=%.2f" %(self.T,self.z,self.g))

    def getSpectrum(self,Teff=5800,z=0.,logg=4.43,vega=False,): #default: Sun at 10pc Teff=5800, z=0, logg=4.4, johnson_v_abs=4.68
        self.johnson_v = S.ObsBandpass('johnson,v')
        self.T=Teff
        self.z=z
        self.g=logg
        if vega: 
            self.star= S.FileSpectrum(S.locations.VegaFile)
            print ' Using default Pysynphot Vega spectrum'
        else: 
            print self.getFiltername(), self.getSpectrumInfo()
            self.star= S.Icat('phoenix',self.T,self.z,self.g) #FLAM surface flux units, i.e. ergs cm-2 s-1 A-1
            #star  = S.BlackBody(6000)
            #star= S.Icat('phoenix',5800,0.,4.4) #FLAM surface flux units, i.e. ergs cm-2 s-1 A-1
            #self.star= self.star.renorm(self.j_v_abs,'vegamag',self.johnson_v) #renorm spectrum to (abs_mag, units system, band)        
        return self.star

    def getBandpass(self): 
        if self.instr=='MIRI':
            filtname='/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/'+self.instr+'/filtersCoro/'+self.filt+'.csv'
            try: 
                self.ifile= open(filtname, 'rU')
                print ' Using filter: ',filtname,'\n'                
            except: 
                print " Wrong or missing filter filename. Choose from:\n"
                filters=glob.glob('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/'+self.instr+'/filtersCoro/F*.csv')
                print [filters[i].split('/')[-1][:-4] for i in xrange(len(filters))],"\n"
                sys.exit()

            self.reader = csv.reader(self.ifile)
            self.fwave=[]
            self.fthru=[]
            self.nrow=0
            for row in self.reader:
                x,y=row
                self.fwave.append(float(x)*1e4)  #Microns in file, so convert to angstroms
                self.fthru.append(float(y)/100.) #Percentage in file, so convert to fraction

            self.filt_wave=np.array(self.fwave)
            self.filt_thru=np.array(self.fthru)
            self.filt_thru[self.filt_thru<0]=0

        elif self.instr=='NIRCam':
            # THESE FILTER FILES INCLUDE: OTE + INSTRUMENT + DBS + QE TRANSMISSION PROFILES, 
            # Provided by J. Stansberry. April 2015
            filtname='/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/'+self.instr+'/filters/'+self.filt+'_dbs_qe.dat'
            try: 
                filt_file = np.loadtxt(filtname)
                print ' Using filter: ',filtname,'\n'
            except: 
                print " Wrong or missing filter filename. Choose from:\n"
                print '/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/'+self.instr+'/filters/'+self.filt+'_dbs_qe.dat'+'\n'
                filters=glob.glob('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/'+self.instr+'/filters/F*.dat')
                print [filters[i].split('/')[-1][:-4] for i in xrange(len(filters))],"\n"
                sys.exit()

            self.filt_wave = filt_file[:,0]*1e4 #Microns in file, so convert to angstroms
            self.filt_thru = filt_file[:,1]
            self.filt_thru[self.filt_thru<0]=0

        # CREATE A "SPECTRUM" FROM THE ARRAYS FILT_WAVE AND FILT_THRU:
        self.bandpass = S.spectrum.ArraySpectralElement(throughput=self.filt_thru, wave=self.filt_wave, waveunits='angstroms')       
        return self.bandpass 



    def getMIRI_QE(self): 
        try: self.ifile= open('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/filtersCoro/MIRI_QE.csv', 'rU')
        except: 
            print " Wrong or missing QE filename."
            sys.exit()

        self.reader = csv.reader(self.ifile)

        self.fwave=[]
        self.fthru=[]
        self.nrow=0
        for row in self.reader:
            x,y=row
            if float(x)*1e4 not in self.fwave: 
                self.fwave.append(float(x)*1e4) #microns in file; convert to angstroms
                self.fthru.append(float(y))

        self.filt_wave=np.array(self.fwave)
        self.filt_thru=np.array(self.fthru)
        self.filt_thru[self.filt_thru<0]=0
        self.ifile.close()
        

        self.bandpass = S.spectrum.ArraySpectralElement(throughput=self.filt_thru,wave=self.filt_wave, waveunits='angstroms')       
        return self.bandpass 


    def getMIRI_OTE_transmission(self):
        #5 to 30 microns; convert to angstroms
        waves = np.arange(5,30+1,0.5)*1e4
        trans = 0.95*np.ones(len(waves)) # see Paul Lightsey figure 4, SPIE 2012 "Optical transmission for the James Webb Space Telescope"

        eOTE = np.asarray( zip(waves,trans) )

        self.filt_wave=eOTE[:,0]
        self.filt_thru=eOTE[:,1]
        self.filt_thru[self.filt_thru<0]=0
        self.bandpass = S.spectrum.ArraySpectralElement(throughput=self.filt_thru,wave=self.filt_wave, waveunits='angstroms') 

        return self.bandpass


    def getMIRI_Germanium(self): 
        try: self.ifile= open('/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/filtersCoro/Germanium.csv', 'rU')
        except: 
            print " Wrong or missing Germanium filename."
            sys.exit()

        self.reader = csv.reader(self.ifile)

        self.fwave=[]
        self.fthru=[]
        self.nrow=0
        for row in self.reader:
            x,y=row
            if float(x)*1e4 not in self.fwave: 
                self.fwave.append(float(x)*1e4) #microns in file; convert to angstroms
                self.fthru.append(float(y))

        self.filt_wave=np.array(self.fwave)
        self.filt_thru=np.array(self.fthru)
        self.filt_thru[self.filt_thru<0]=0
        self.ifile.close()
        

        self.bandpass = S.spectrum.ArraySpectralElement(throughput=self.filt_thru,wave=self.filt_wave, waveunits='angstroms')       
        return self.bandpass 




def AddBKGRDPoissonNoise(data,Bkgrd,absVal=False):
        """
        Generate Poisson noise and add it to the internal image.

        SEE BRIAN YORK'S JwstImage.py routine

        Parameters
        ----------
        absVal: bool, optional
            Take absolute value of `noiseData`.
        """

        # 1. Gaussian distribution - shot noise.
        # 2. Poisson = Gaussian * SQRT(image)
        noiseData = np.random.normal(size=data.shape) * np.sqrt(Bkgrd)

        # Avoid negative noise
        if absVal: noiseData = np.abs(noiseData)
        
        data += noiseData
        del noiseData
        return data



def prepareSpectrum(instr,filt,Teff=5800.,z=0.0,logg=4.44,radius=1.,distance=10.,vega=False,exptime=1.0,webbpsf=False,**kwargs):
    if instr=='NIRCam': filt=filt[:5] # Strip the mask string off, e.g. F210M-MASK210R
    if instr=='MIRI': 
        xdim,ydim = (64,64)
    elif instr=='NIRCam':
        nircamMode = None
        try:
            w_fil = int(filt[1:-1])
        except:
            print filterName, 'has no wavelength info'
            sys.exit()
        if w_fil >= 235:
            xdim,ydim = (109,109)
            nircamMode = 'long'
        else:
            xdim,ydim = (222,222)
            nircamMode = 'short'

    # CREATE SPECTRUM OBJECT:
    spec=Spectrum(instr,filt)

    # CREATE FILTER BANDPASS FOR FILTER NAME:
    bandpass=spec.getBandpass()  

    # CREATE SPECTRUM BASED ON Teff, Z, and LOGG:
    star=spec.getSpectrum(Teff=Teff,z=z,logg=logg,vega=vega) # Returns Angstroms, Flam

    # SCALE FOR DISTANCE & CONVERT UNITS TO: photons s^-1 cm^-2 A^-1
    if vega==False: star*=(radius*2.254e-8/distance)**2 #distance in units of parsec
    star.convert("photlam") 

    # CREATE OBSERVATION (STAR times BANDPASS):
    obs=S.Observation(star,bandpass) 
    #print 'units: ',obs.waveunits, obs.fluxunits


    # ADD QE/GERMANIUM/OTE TRANSMISSION TO OBSERVATION:
    if webbpsf:
        print "Using images created by WebbPSF"
        if instr == 'MIRI': 
            bp_QE   = spec.getMIRI_QE()
            obs     = S.Observation(obs,bp_QE)
            bp_Germ = spec.getMIRI_Germanium()
            obs     = S.Observation(obs,bp_Germ)
            bp_OTE  = spec.getMIRI_OTE_transmission()
            obs     = S.Observation(obs,bp_OTE)
        else: pass # FOR NIRCam, QE/OTE PROFILES ARE INCLUDED IN FILTER TRANSMISSION CURVES
    else: pass     # FOR IMAGES CREATED WITH MATHEMATICA, QE/GERMANIUM/OTE ARE ALREADY INCLUDED


    # INTEGRATE TOTAL NUMBER OF PHOTONS OVER BANDPASS:
    ncounts=np.trapz(obs.flux,obs.wave)*DefaultSettings.JWSTREFS['area']  #multiply by collecting area

    return ncounts, spec, xdim, ydim



def convertnew(spec, ncounts, in_path, out_path, xdim, ydim, Teff=5800.,z=0.0,logg=4.44,radius=1.,distance=10., exptime=1.0, run=1, clobber=False, nonoise=False, **kwargs):

    # LOOP OVER ALL THE ORIGINAL IMAGES: SCALE THEM AND ADD NOISE SOURCES:
    runNumber=run
    string='run'+str(runNumber)+'_'
    listfiles=os.listdir(in_path)
    infiles=[i for i in listfiles if 'PSF' in i]
    inputfiles=[i for i in infiles if 'run'+str(runNumber)+'_' in i]

    ic=0
    for ifile in inputfiles:
        print ifile
        img= JwstImage.JwstImage.FromFitsFile(in_path+ifile,instrument=spec.instr,filterName=spec.filt,out_x=xdim,out_y=ydim,in_path=in_path,out_path='.')
        try: 
            background = DefaultSettings._BGLEV[img.instrument][img.filterName]
        except:
            if ic==0: print '  --> Filter '+img.filterName+" not in DefaultSettings.py for background levels. Using Lajoie's numbers.\n"
            background = MY_BGLEVEL[img.instrument][img.filterName]

        img.exptime=exptime #that actually multiples the image by the exposure time

        if nonoise==False:
            print '   --> Adding noise to PSF'

            img.__imul__( ncounts )            

            img.AddPoissonNoise()

            #print "!!! NO BACKGROUND NOISE !!! "
            #img.__simplyAdd__(background*exptime)
            img._data=AddBKGRDPoissonNoise(img._data,background*exptime)

            if instr=='MIRI':
                img.AddDetectorNoise(instrumentMode='miri_fast')
                img.AddErrorFlat('JWST_Sim/'+DefaultSettings._FILE_ERR_FLAT[instr])
                img.AddErrorDark('JWST_Sim/'+DefaultSettings._FILE_ERR_RDRK[instr])

            elif instr=='NIRCam': 
                img.AddDetectorNoise(instrumentMode=instr.lower())
                img.AddErrorFlat('JWST_Sim/'+DefaultSettings._FILE_ERR_FLAT[instr])
                img.AddErrorDark('JWST_Sim/'+DefaultSettings._FILE_ERR_RDRK[instr+'_'+nircamMode])

            img.AddErrorCosmic(DefaultSettings._PIXSIZE[instr])

        else: 
            img.__imul__( ncounts )

    #   DETECTOR SATURATION???
    #    img._data[img._data>SATUR_LVL[img.instrument]]=SATUR_LVL[img.instrument]

        img._header.add_history("LAJOIE: Spectrum scaled using AddNoisetoPSF.py")
        string = "LAJOIE: Spectrum Teff= %d  Log g=%4.2f z=%5.2f R=%4.2f D=%d exptime=%d" %(Teff, logg, z, radius, distance, exptime)
        img._header.add_history(string)
        try: 
            hdu=astropy.io.fits.PrimaryHDU(data=img._data, header=img._header)
            hdu.writeto(out_path+'Scaled_'+ifile, clobber=clobber)#, img._data, img._header, clobber=clobber)
            #writeto(out_path+'Scaled_'+ifile, img._data, img._header, clobber=clobber)
        except: print '  --> File Scaled_'+ifile+' already exists. Use --clobber to overwrite'
        del img
        ic+=1

#-----------------------------------#
# MAIN Method
#-----------------------------------#
if __name__ == "__main__":

    #-----------------------------------#
    # PARSE COMMAND LINE ARGUMENTS:
    #-----------------------------------#
    parser = argparse.ArgumentParser(description="Add simple noise sources to pre-generated PSFs")

    parser.add_argument("-I"       , metavar="Instrument",required=True, help="JWST instrument")
    parser.add_argument("-f"       , metavar="FILTER",required=True, help="Filter to use")
    parser.add_argument("-step"    , required=True, default=20, help="Dither grid step size (default: 20 mas)")
    parser.add_argument("-jitter"  , required=True, default=0 , help="sigma Jitter (default: 0 mas)")
    parser.add_argument("--webbpsf", action="store_true", default=False,   help="use if images were created with WebbPSF")
    parser.add_argument("--t"      , type=float, default=5800,  help="effective temperature (default: solar T_eff=5800 K)")
    parser.add_argument("--z"      , type=float, default=0.0 ,  help="metallicity (default: solar z=0)")
    parser.add_argument("--g"      , type=float, default=4.44,  help="surface gravity (default: solar log_g=4.44)")
    parser.add_argument("--R"      , type=float, default=1 ,    help="radius of star (default: 1 Rsun )")
    parser.add_argument("--dist"   , type=float, default=10,    help="distance to star parsec (default: 10 pc)")
    parser.add_argument("--exptime", type=float, default=1.0 ,  help="exposure time (default: 1 sec)")
    parser.add_argument("--run"    ,             default="all", help="run number (default: all)")
    parser.add_argument("--vega"   , action="store_true", default=False, help="use Vega spectrum (default: False)")
    parser.add_argument("--clobber", action="store_true", default=False, help="overwrite output files (default: False)")
    parser.add_argument("--nonoise", action="store_true", default=False, help="to turn off noise sources (default: add noise)")
    parser.add_argument("--rms"    , type=int,            default=400,   help="select which OPD rms to use, if available (default: 400 nm)")

    args = parser.parse_args()


    instr=(args.I).upper()
    if instr=='NIRCAM': instr='NIRCam'
    filt=(args.f).upper() 
    kwargs={'webbpsf':args.webbpsf,
            'Teff':args.t,
            'z':args.z,
            'logg':args.g,
            'radius' :args.R,
            'distance':args.dist,
            'vega':args.vega,
            'exptime':args.exptime,
            'clobber':args.clobber,
            'nonoise':args.nonoise
    }


    #-----------------------------------#
    # SET UP PATHS BASED ON ARGUMENTS:
    #-----------------------------------#
    path    = "/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/MIRI/Dither-LOCI/Results/"
    in_path = path+instr+'/'+args.f.upper()+"/DitherStep"+args.step+"/Jitter"+args.jitter+"/Originals"
    out_path= path+instr+'/'+args.f.upper()+"/DitherStep"+args.step+"/Jitter"+args.jitter+"/"

    if   args.vega==False: out_path+='T'+str(args.t)+'_logg'+str(args.g)+'_z'+str(args.z)+'_R'+str(args.R)+'_D'+str(args.dist)+'_t'+str(args.exptime)
    elif args.vega==True : out_path+='Vega'+'_t'+str(args.exptime)

    if args.nonoise==True: out_path+='_NoNoise'

    if args.rms!=400: 
        in_path += str(args.rms)+'nm/' #-Grid25pts/'
        out_path+= '_'+str(args.rms)+'nm'
    else:
        in_path +='/'
        out_path+='/'

    print ' \nInput original PSFs: ',in_path.split("Dither-LOCI")[1]
    if os.path.isdir(in_path)==False: 
        print "   --> Input directory does not exist."
        sys.exit()

    if os.path.exists(out_path): 
        print '\nOutput folder',out_path.split("Dither-LOCI")[1],'already exists\n'
    else: 
        os.mkdir(out_path)
        print '\nCreating folder',out_path.split("Dither-LOCI")[1]+'\n'


    #-----------------------------------#
    # CALCULATE NRUNS IN DIRECTORY
    #-----------------------------------#
    infiles= glob.glob(in_path+"*ScienceTarget*")
    nruns  = len(infiles)


    #-----------------------------------#
    # CONVERSION OF PSFs:
    #-----------------------------------#
    

    scalefactor, spectrum, xdim, ydim = prepareSpectrum(instr, filt, **kwargs)
    print "scale factor phot s^-1 :", scalefactor

    if args.run != 'all':  
        kwargs['run']=int(args.run)
        #convert(in_path,out_path,instr,filt,**kwargs)
        convertnew(spectrum, scalefactor, in_path,out_path, xdim, ydim, **kwargs)
    else:
        for i in xrange(1,nruns+1):
            kwargs['run']=i
            convertnew(spectrum, scalefactor, in_path,out_path, xdim, ydim, **kwargs)
            #convert(in_path,out_path,instr,filt,**kwargs)


