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
import scipy.ndimage.interpolation as scint
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
            self.star= S.Icat('ck04models',self.T,self.z,self.g) #FLAM surface flux units, i.e. ergs cm-2 s-1 A-1
            #star  = S.BlackBody(6000)
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




def findCountRate(instr,filt,Teff=5800.,z=0.0,logg=4.44,radius=1.,distance=10.,vega=False,exptime=1.0,webbpsf=False,**kwargs):

    if instr=='NIRCam': filt=filt[:5] # Strip the mask string off, e.g. F210M-MASK210R

    # CREATE SPECTRUM OBJECT:
    spec=Spectrum(instr,filt)

    # CREATE FILTER BANDPASS FOR FILTER NAME:
    bandpass=spec.getBandpass()  

    # CREATE SPECTRUM BASED ON Teff, Z, and LOGG:
    star=spec.getSpectrum(Teff=Teff,z=z,logg=logg,vega=vega) # Returns Angstroms, Flam

    # SCALE FOR DISTANCE & CONVERT UNITS TO: photons s^-1 cm^-2 A^-1
    if vega==False:  star*=(radius*2.254e-8/distance)**2 #Distance in units of parsec
    star.convert("photlam") 


    # CREATE OBSERVATION (STAR times BANDPASS):
    obs=S.Observation(star,bandpass, binset=star.wave) 
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

    return ncounts, spec



def addNoise(spec, ncounts, in_path, out_path, Teff=5800.,z=0.0,logg=4.44, radius=1.,distance=10., 
             exptime=1.0, run=1, roll=0, clobber=False, nonoise=False, planets=False, vega=False, **kwargs):

    if spec.instr=='MIRI': nircamMode = None
    elif spec.instr=='NIRCam':
        try:
            w_fil = int(spec.filt[1:-1])
            if w_fil >= 235: nircamMode = 'long'
            else: nircamMode = 'short'
        except:
            print spec.filt, 'has no wavelength info'
            sys.exit()


    # LOOP OVER ALL THE ORIGINAL IMAGES. SCALE THEM AND ADD NOISE SOURCES:

    runNumber=run
    string='run'+str(runNumber)+'_'
    listfiles=os.listdir(in_path)
    infiles=[i for i in listfiles if 'PSF' in i]
    inputfiles=[i for i in infiles if 'run'+str(runNumber)+'_' in i]

    unocculted=glob.glob(in_path+"/*run"+str(runNumber)+"_"+"*Unocculted*"+"*.fits")
    hdu0=pyfits.open(unocculted[0])
    xdim, ydim = hdu0[0].data.shape


    # ADD PLANETS IF REQUESTED:
    if planets:
        planet=hdu0[0].data * ncounts * exptime   
        planets_img, planets_mags, planets_seps = addPlanets(spec.instr, spec.filt, planet, -roll)


    for ifile in inputfiles:
        print '\n--> Processing', ifile
        img= JwstImage.JwstImage.FromFitsFile(in_path+ifile,instrument=spec.instr,filterName=spec.filt,out_x=xdim,out_y=ydim,in_path=in_path,out_path='.')
        try: 
            background = DefaultSettings._BGLEV[img.instrument][img.filterName]
        except:
            print '  Filter '+img.filterName+" not in DefaultSettings.py for background levels. Using Lajoie's numbers."
            background = MY_BGLEVEL[img.instrument][img.filterName]


        # MULTIPLY ORIGINAL IMAGE BY NCOUNTS/SECOND:
        img.__imul__( ncounts )

        if "ScienceTarget" in ifile:
            star_throughput = np.sum(img._data)
            print "star_throughput",star_throughput, ncounts

        # MULTIPLY IMAGE BY EXPOSURE TIME:
        img.exptime=exptime

        # ADD PLANETS IF DESIRED:
        if planets and "ScienceTarget" in ifile: 
            print "  Injecting planets to", ifile
            img._data += planets_img

        # ADD "RESIDUAL" NOISE SOURCES:
        if nonoise==False:
            print '  Adding noise to PSF'

            # ADD BACKGROUND TO IMAGE BEFORE CALCULATING TOTAL POISSON NOISE:
            img.__simplyAdd__(background*exptime)

            # ADD POISSON NOISE TO IMAGE:
            img.AddPoissonNoise()

            # ADD "RESIDUAL" DARK, READNOISE, AND FLAT NOISE TO IMAGE:
            if   instr=='MIRI':   
                img.AddErrorDark('JWST_Sim/'+DefaultSettings._FILE_ERR_RDRK[instr])
                img.AddDetectorNoise(instrumentMode='miri_fast')
                img.AddErrorFlat('JWST_Sim/'+DefaultSettings._FILE_ERR_FLAT[instr])
            elif instr=='NIRCam': 
                img.AddErrorDark('JWST_Sim/'+DefaultSettings._FILE_ERR_RDRK[instr+'_'+nircamMode])
                img.AddDetectorNoise(instrumentMode=instr.lower())
                img.AddErrorFlat('JWST_Sim/'+DefaultSettings._FILE_ERR_FLAT[instr])

            # ADD "RESIDUAL" COSMIC RAYS TO IMAGE:
            img.AddErrorCosmic(DefaultSettings._PIXSIZE[instr])
            
            # SUBTRACT NOMINAL BACKGROUND TO MIMICK PROCESSING; POISSON NOISE REMAINS:
            img.__simplyAdd__(-background*exptime)



        # UPDATE IMAGE HEADER & SAVE:
        img._header.add_history("Spectrum scaled using AddNoisetoPSF.py")
        if vega: string = "%s: Using Pysynphot spectrum of Vega, exptime=%d" %(os.getlogin().upper(), exptime)
        else:    string = "%s: Spectrum Teff= %d  Log g=%4.2f z=%5.2f R=%4.2f D=%d exptime=%d" %(os.getlogin().upper(), Teff, logg, z, radius, distance, exptime)
        img._header.add_history(string)

        if planets and "ScienceTarget" in ifile:
            string = "%s: Planets mags B, C, D, E: %5.2f %5.2f %5.2f %5.2f" %(os.getlogin().upper(), planets_mags[spec.filt][0], planets_mags[spec.filt][1], planets_mags[spec.filt][2], planets_mags[spec.filt][3])
            img._header.add_history(string)
            string = "%s: Planet B separation: %s " %(os.getlogin().upper(), str(planets_seps["B"]))
            img._header.add_history(string)
            string = "%s: Planet C separation: %s " %(os.getlogin().upper(), str(planets_seps["C"]))
            img._header.add_history(string)
            string = "%s: Planet D separation: %s " %(os.getlogin().upper(), str(planets_seps["D"]))
            img._header.add_history(string)
            string = "%s: Planet E separation: %s " %(os.getlogin().upper(), str(planets_seps["E"]))
            img._header.add_history(string)
        try: 
            hdu=astropy.io.fits.PrimaryHDU(data=img._data, header=img._header)
            hdu.writeto(out_path+'/Scaled_'+ifile, clobber=clobber)
        except: print '  --> File Scaled_'+ifile+' already exists. Use --clobber to overwrite'
        del img

    return star_throughput



def addPlanets(instr, filt, planet, rotation):
    # VALUES FOR HR 8799, PLANETS B, C, D, and E:
    # - BOCCALETTI PASP PAPER
    # - ApJ 795 : 133 Currie et al. 2015
    mags={"F1065C":[9.73, 9.12, 9.12, 15.0],
          "F1140C":[9.16, 8.82, 8.82, 15.0],
          "F1550C":[9.02, 8.71, 8.71, 15.0],
          "F210M" :[10.3, 9.45, 9.30, 9.30],
          "F430M" :[10.3, 9.45, 9.30, 9.30],
          "F460M" :[10.3, 9.45, 9.30, 9.30]}

    seps={"B":[0.706, -1.563], "C":[0.765, 0.558], "D":[-0.529, 0.323], "E":[-0.09, 0.366] }


    # DETERMINE PIXEL SIZE FOR INSTRUMENT/CHANNEL:
    if instr=="MIRI": pixelSize = 0.11
    if instr=="NIRCam":
        if filt=="F210M": pixelSize = 0.032
        else: pixelSize = 0.065

    # ROTATE THE PLANETS AROUND TO AVOID BAR/4QPM:
    if rotation != 0:
        angle = rotation * np.pi/180.

        delta0, alpha0=[], []
        for key in sorted(seps):
            delta0.append(seps[key][0])
            alpha0.append(seps[key][1])

        delta0=np.asarray(delta0)
        alpha0=np.asarray(alpha0)

        deltas= (alpha0*np.cos(angle) - delta0*np.sin(angle)) / pixelSize
        alphas= (alpha0*np.sin(angle) + delta0*np.cos(angle)) / pixelSize
        
        shiftedB = scint.shift(planet, (alphas[0], deltas[0]) )
        shiftedC = scint.shift(planet, (alphas[1], deltas[1]) )
        shiftedD = scint.shift(planet, (alphas[2], deltas[2]) )
        shiftedE = scint.shift(planet, (alphas[3], deltas[3]) )

    else:
        shiftedB = scint.shift(planet, (seps["B"][0]/pixelSize, seps["B"][1]/pixelSize))
        shiftedC = scint.shift(planet, (seps["C"][0]/pixelSize, seps["C"][1]/pixelSize))
        shiftedD = scint.shift(planet, (seps["D"][0]/pixelSize, seps["D"][1]/pixelSize))
        shiftedE = scint.shift(planet, (seps["E"][0]/pixelSize, seps["E"][1]/pixelSize))


    # SCALE PLANET IMAGE BY MAGNITUDE DIFFERENCE:
    planets = shiftedB*10**(-mags[filt][0]/2.5) + shiftedC*10**(-mags[filt][1]/2.5) + shiftedD*10**(-mags[filt][2]/2.5) + shiftedE*10**(-mags[filt][3]/2.5)

    return planets, mags, seps




#-----------------------------------#
# MAIN Method
#-----------------------------------#
if __name__ == "__main__":
    np.random.seed(0)

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
    parser.add_argument("--planets", action="store_true", default=False, help="Inject HR8799 planets (default: False)")
    parser.add_argument("--roll"   , type=float,          default=0.,    help="Roll the planets clockwise in degrees (default: 0 deg.)")

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
            'nonoise':args.nonoise,
            'planets':args.planets,
            'roll'   :args.roll
    }

    #-----------------------------------#
    # SET UP PATHS BASED ON ARGUMENTS:
    #-----------------------------------#
    path    = "/Users/lajoie/Documents/Work/Projects/JWST/Simulations/Coronagraphs/Dither-LOCI/Results/"
    in_path = path+instr+'/'+args.f.upper()+"/DitherStep"+args.step+"/Jitter"+args.jitter+"/Originals"
    out_path= path+instr+'/'+args.f.upper()+"/DitherStep"+args.step+"/Jitter"+args.jitter+"/"

    if   args.vega==False: out_path+='T'+str(args.t)+'_logg'+str(args.g)+'_z'+str(args.z)+'_R'+str(args.R)+'_D'+str(args.dist)+'_t'+str(args.exptime)
    elif args.vega==True : out_path+='Vega'+'_t'+str(args.exptime)

    if args.nonoise==True: out_path+='_NoNoise'

    if args.rms!=400: 
        in_path += str(args.rms)+'nm/' #_WebbPSF/' #-Grid25pts/'
        out_path+= '_'+str(args.rms)+'nm' #_WebbPSF'
        if args.planets: out_path+="_Planets/"
    else:
        in_path +='/'
        if args.planets: out_path+="_Planets"
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
    scalefactor, spectrum = findCountRate(instr, filt, **kwargs)
    print "Scale factor phot s^-1 :", scalefactor
    
    if args.run != 'all':  
        kwargs['run']=int(args.run)
        star_throughput = addNoise(spectrum, scalefactor, in_path,out_path, **kwargs)
    else:
        for i in xrange(1,nruns+1):
            kwargs['run']=i
            star_throughput = addNoise(spectrum, scalefactor, in_path,out_path, **kwargs)


    with open(out_path+"/InputStar.txt", "w") as text_file:
        text_file.write("Total count/sec: %f\n" % (scalefactor)) #star_throughput))
        text_file.write("Distance: %f\n" % (args.dist))

