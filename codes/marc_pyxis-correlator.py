import Pyxis
import MSUVDISTANCE
import pyrap.tables
from Pyxis.ModSupport import *
import mqt,lsm
# this has some useful templates we can borrow
import imager,ms
import ClassMS
import ClassMS1
import Test
import Source
import os
import sys
import MSResampler
import scipy.special
import pyfits
import pylab
import numpy as np
_MAKEMS_DEFAULTS = dict(
		WriteAutoCorr = True,
		NFrequencies = 1,	
		StepTime = 10,
		StartFreq = 1.4e+9,
		NParts = 1,
		MSDesPath = '.',
		WriteImagerColumns = True,
		AntennaTableName = 'VLACKAT_ANTENNA',
		Declination = '-45.0.0',
		NBands = 1,
		RightAscension = '00:00:00',
		StartTime = '2011/11/16/17:00',
		StepFreq = 125000.0,
		NTimes = 1
 		);


## create MAKEMS_XX globals from dict
for key,value in _MAKEMS_DEFAULTS.iteritems():
  def_global("MAKEMS_%s"%key.upper(),value,"%s option to makems.cfg"%key);

def make_new_ms(Ntimeslots=None, steptim=None, Nfreq=None,StepFreq=None, msname="MeerKAT-snapshot-21cm.MS",**kw):
  	"""runs makems to create a new MS. Keywords can override global MAKEMS_XXX options.""";
  	# load defaults from MAKEMS_xxx variables
  	# make dict of makems.cfg settings using globals and keywords
  	# note that globals and keywords are uppercase
	opts = {}
  	for key in _MAKEMS_DEFAULTS.iterkeys():
    		ku = key.upper();
    		opts[key] = kw.get(ku,globals().get("MAKEMS_%s"%ku,""));
  	# set MS name in options
  	opts['MSName'] = v.MS =msname or "MeerKAT-snapshot-21cm.MS"
	opts['StepTime'] = steptim
	opts['NTimes'] = Ntimeslots 
	opts['NFrequencies']=Nfreq
	opts['StepFreq']=StepFreq
  	# write config filebased on options dict
  	conffile = II('${msname:BASE}.makems.cfg');
  	conf = file(conffile,'w');
  	for kw,value in opts.iteritems():
    		conf.write("%s=%s\n"%(kw,value));
	info(opts)
  	conf.close();
  	x.sh("makems $conffile");
	print "ici avant remove"
  	x.sh("rm -fr $msname");
	print("ici avant move")
 	x.sh("mv ${msname}_p0 $msname");

document_globals(make_new_ms,"MAKEMS_* MS");
def_global("SIM_CONFIG_SECTION","sim","TDL config section");
MS_TARBALL_DIR = os.path.expanduser(".");


def reset_ms(msname="$MS",Ntimslots=None,steptim=None,Nfreq=None,StepFreq=None):
	info(">>> $msname")
        print (">>>>>>>>>>>>>>>>>>>>>>>>>")
  	if os.path.exists(msname):
    		x.sh("rm -fr $msname")
  	make_new_ms(Ntimeslots=Ntimslots, steptim=steptim, Nfreq=Nfreq, StepFreq=StepFreq, msname=msname)
  	# Add Column to store model data with pointing error
  	#xo.sh('ipython addCollToMs.py $MSi MODEL_PERR_DATA')


def_global("SIM_CONFIG_SECTION","sim","TDL config section");

IMGLABEL = "s1"; 
IMGNAME_Template = "plots-${MS:BASE}/${MS:BASE}-$IMGLABEL.dirty.fits"
#MSHI_Template = "${MS:BASE}-hires.MS"
IMGNAME_Template = "plots-${MS:BASE}/${MS:BASE}-$IMGLABEL.dirty.fits"

imager.DIRTY_IMAGE_Template = "${OUTFILE}-${COUNTER}-${METHOD}.dirty.fits"


# ----------All outputs are save in the directory DATA---------------------------------

fov = 1.2
def sinc1 (x,deg=(fov*np.pi**2)/(180.)):
      x1 = x*deg;
      wf = np.sin(x1)/x1;
      wf[x1==0] = 1;
      return wf;
def sinc (x,y,deg=(fov*np.pi**2)/(180.)):
      x1 = np.sqrt(x**2+y**2)*deg;
      wf = np.sin(x1)/x1;
      wf[x1==0] = 1;
      return wf;
def sinc2 (x,y,deg=(fov*np.pi**2)/(180.)):
      x1,y1 = x*deg,y*deg;
      wx,wy = (np.sin(x1)/x1),(np.sin(y1)/y1);
      wx[x1==0] = 1;
      wy[y1==0] = 1;
      return wx*wy;
def airy (x,a=2*np.pi*fov*np.pi/180):
      r = x;
      w = scipy.special.j1(r*a)/(r*a);
      w[r==0] = 0.5;
      return w;
def butter (x1,y1,p=5.,deg=11.):#11                   
          x1 = np.sqrt(x1**2+y1**2);
          wx = (1 + (x1/deg)**(2*p))**(-1);
          wx[x1==0] = 1;
          return wx;
def MoveSrcRadius_normal():
        Ntimeslots=400 # Number of timeslots of the Hight resolution MS
        StepTime=0.1  #StepTime of timeslots of the Hight resolution MS
        NFreq= 150     # Number of Frequency of timeslots of the Hight resolution MS
        StepFreq=125e3  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        dtime=100     # integration or the number of sample where the signal is cut-off
        NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
        NsampleFreq=NFreq # default  equal to Number of frequency hight resolution or less than
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400     # the lenght in arcmin of the x-axis 
        dfreq = 50#50
        list1 = np.arange(180,180.2,180)#(15,121,15.)#(5.5,5.6,5.5)
        steptime = np.arange(0.5,5.6,0.5)
        Lnoise = np.zeros(len(list1))
        StepFreq=np.arange(125000,810000,125000/2.)#StepFreq=np.arange(125000,1250000,125000/2.)#125000.0
        Lmeter = np.zeros(len(steptime))# vector that save the flux of baseline per meter
        
        Rmeter = np.zeros(len(steptime))
        Lav = np.zeros(len(steptime))   # vector that save the flux of simple averaging
        Rav = np.zeros(len(steptime))
        locatsrc = [30,40.,30.]
        ORG_MS = v.MS
        for i  in range(len(list1)):
          for k in range(len(steptime)):
            hires = ORG_MS.split(".MS")[0]+"-hires.MS"

            reset_ms(hires, NsampleTime, steptime[k], NsampleFreq, StepFreq[k])
            reset_ms(ORG_MS,1, steptime[k]*dtime,1,StepFreq[k]*dfreq)  
            options = {}
            options['gridded_sky.grid_m0'] = list1[i]
            options['ms_sel.msname'] = hires;
            mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                    config="tdlconf.profiles",section="Sim_source_radius",
                    options=options);
            mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=100,ntime=100,freq0=50,nfreq=50)
            arrays=mshi.boxcar (dtime,dfreq)
            MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA") 
            center_min = -45*60+list1[i];
            center_deg = math.ceil(center_min/60);
            center_min = abs(center_min - center_deg*60);
            imager.make_image(column = 'CORRECTED_DATA',
                              phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                              restore = False, dirty = True, restore_lsm = False,
                              weight = "natural");
            f,r,sta = Test.giveImage(imager.DIRTY_IMAGE)
            Lav[k] = f
            Rav[k]=StepFreq[k]
            mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=150,ntime=100,freq0=50,nfreq=50)
                        #arrays=mshi.window (sinc2,dtime,dfreq,dump=None,dumpfig=None)
            arrays=mshi.overlap_window (butter,dtime,dfreq,overlap_time=100,overlap_freq=50,dump=None,dumpfig=None)
            MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA")
            imager.make_image(column = 'CORRECTED_DATA',
                              phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                              restore = False, dirty = True, restore_lsm=False,
                              weight = "natural");
            f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
            Lmeter[k] = f
            Rmeter[k] = steptime[k]
          np.save("DATA/Lmeter%f"%(list1[i]/60),Lmeter)
	  np.save("DATA/Rav%f"%(list1[i]/60),Rav)
          np.save("DATA/Lav%f"%(list1[i]/60),Lav)
          np.save("DATA/Rmeter%f"%(list1[i]/60),Rmeter)



def Noise ():
  Ntimeslots=1200 # Number of timeslots of the Hight resolution MS
  StepTime=1.5  #StepTime of timeslots of the Hight resolution MS
  NFreq= 1   # Number of Frequency of timeslots of the Hight resolution MS
  StepFreq=125e3  # Step Frequency of Frequency of timeslots of the Hight resolution MS
  dtime=600     # integration or the number of sample where the signal is cut-off
  NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
  NsampleFreq=NFreq # default  equal to Number of frequency hight resolution or less than
  stepR=10.      # the step radius in arcmin of the x-axis
  lenghtR=400     # the lenght in arcmin of the x-axis 
  dfreq = 1#50

  ORG_MS = v.MS
  hires = ORG_MS.split(".MS")[0]+"-hires.MS"
  reset_ms(hires, NsampleTime, StepTime, NsampleFreq, StepFreq)
  reset_ms(ORG_MS, 1, StepTime*dtime,1, StepFreq*dfreq)
  #reset_ms(ORG_MS, 100/dtime, StepTime*dtime, (50)/dfreq, StepFreq*dfreq)
  list1 = np.arange(0,lenghtR,stepR)
  Lav = np.zeros(len(list1))
  Rav = np.zeros(len(list1))
  options={}
  options["noise_stddev"]=1.
  mqt.run("turbo-sim.py",job="_tdl_job_1_simulate_MS",config="tdlconf.profiles",section="sim_none_source",options=options);
  mshi =MSUVDISTANCE.MSUVDISTANCE(hires+"/",column="DATA",time0=600,ntime=600,freq0=0,nfreq=1)
  arrays , time1,dtimemin=mshi.boxcar(dfreq,dtimemax=600)#boxcar (dfreq,dtimemax=100);
  i =0
  if i ==0:
                reset_ms(ORG_MS, time1, StepTime*dtimemin,dfreq, StepFreq*dfreq)
                os.system("addbitflagcol VLAC-snapshot-21cm.MS")
  mshi.flaggindata(v.MS,arrays,time1,column="CORRECTED_DATA",irgg=i)
  #arrays , time1,dtimemin=mshi.overlap_window (sinc1, dfreq,dtime,50)
  #mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=0,ntime=100,freq0=0,nfreq=0)
  #arrays=mshi.boxcar (dtime,dfreq)
  #MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA") 
  #mshi.flaggindata(v.MS,arrays,time1,column="CORRECTED_DATA",irgg=1)
  #mshi = ClassMS.ClassMS("VLAC-snapshot-21cm-hires.MS/","UVW","DATA",dtime, NsampleTime)
  #dataav,ntim=mshi.NormalAvgALBL();
  #Test.SaveVis("VLAC-snapshot-21cm.MS/",dataav,"DATA");
  imager.make_image(weight="natural",column='CORRECTED_DATA', restore=False, dirty=True, restore_lsm=False);
  nav=Source.getSTD(imager.DIRTY_IMAGE)

  #mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=150,ntime=100,freq0=50,nfreq=50)
  #arrays , time1,dtimemin=mshi.window(sinc1,dfreq,dtimemax=100)#boxcar (dfreq,dtimemax=100);
  arrays , time1,dtimemin=mshi.overlap_window (sinc1, dfreq,600,0)
  if i ==0:
                reset_ms(ORG_MS, time1, StepTime*dtimemin,dfreq, StepFreq*dfreq)
                os.system("addbitflagcol VLAC-snapshot-21cm.MS")
  mshi.flaggindata(v.MS,arrays,time1,column="CORRECTED_DATA",irgg=i)
  #arrays=mshi.overlap_window (butter,dtime,dfreq,overlap_time=150,overlap_freq=50,dump=None,dumpfig=None)
  #MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA")
  #imager.make_image(weight="natural",column='CORRECTED_DATA', restore=False, dirty=True, restore_lsm=False);
 # mshi.flaggindata(v.MS,arrays,time1,column="CORRECTED_DATA",irgg=1)
  imager.make_image(weight="natural",column='CORRECTED_DATA', restore=False, dirty=True, restore_lsm=False);
  nm=Source.getSTD(imager.DIRTY_IMAGE)
  np.save("DATA/avg",nav)
  np.save("DATA/meter",nm)
  info("averaging noise%f"%nav)
  info("meter noise%f"%nm)
  info("nm/na=%f"%(nm/nav))


#--------This function Move a source from the phase center to a maximun position that you precised and averaged over time and frequency-------
def Averaging_OverTimeXFreq():
  # defined your parameters here

        Ntimeslots=300 # Number of timeslots of the Hight resolution MS
        StepTime=1.5  #StepTime of timeslots of the Hight resolution MS
        NFreq=1     # Number of Frequency of timeslots of the Hight resolution MS
        StepFreq=1000  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        dtime=100     # integration or the number of sample where the signal is cut-off
        NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
        NsampleFreq=NFreq # default  equal to Number of frequency hight resolution or less than
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400     # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, NsampleTime, StepTime, NsampleFreq, StepFreq)
        reset_ms(ORG_MS, NsampleTime/dtime, StepTime*dtime, NsampleFreq/dtime, StepFreq*dtime)
        list1 = np.arange(0,lenghtR,stepR)
        Lav = np.zeros(len(list1))
        Rav = np.zeros(len(list1))
        
        for i  in range(len(list1)):
          options = {}
          options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
          center_min = -45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          mshi = ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime, NsampleTime, NsampleFreq)
          dataav,ntim = mshi.NormalAvgALBL();
          Test.SaveVis(v.MS,dataav,"CORRECTED_DATA");
          imager.make_image(column = 'CORRECTED_DATA',
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lav[i] = f
          Rav[i] = list1[i]/60
          np.save("DATA/Lav%i"%(int(dtime*StepTime)),Lav)
          np.save("DATA/Rav%i"%(int(dtime*StepTime)),Rav)

 
#--------This function Move a source from the phase center to a maximun position that you precised and convolve over time and frequency-------
def Convolve_OverTimeXFreq():
  # defined your parameters here

        Ntimeslots=400 # Number of timeslots of the Hight resolution MS
        StepTime=1.5  #StepTime of timeslots of the Hight resolution MS
        NFreq=400     # Number of Frequency of timeslots of the Hight resolution MS
        StepFreq=1000  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        Fov=2.       # Fov of the Baseline dependent filter
        dtime=100     # integration or the number of sample where the signal is cut-off
        NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
        NsampleFreq=NFreq # default  equal to Number of frequency hight resolution or less than
        Overlap=400     # the number of samples extends, Overlap<= Ntimeslots + dtime 
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400     # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, NsampleTime, StepTime, NsampleFreq, StepFreq)
        reset_ms(ORG_MS, NsampleTime/dtime, StepTime*dtime, NsampleFreq/dtime, StepFreq*dtime)
        deg = (Fov*np.pi**2)/(180.)
        sinc1 = lambda x : (np.sin(x*deg))/(x*deg)
        list1 = np.arange(0,lenghtR,stepR)
        Lmeter = np.zeros(len(list1))
        Rmeter = np.zeros(len(list1))
        
        for i  in range(len(list1)):
          options = {}
          options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
          center_min = -45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          mshi = ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime, NsampleTime, NsampleFreq, Overlap)
          datatimeter,time = mshi.ConvolveALBL(sinc1,1)
          Test.SaveVis(v.MS,datatimeter,"CORRECTED_DATA");
          imager.make_image(column = 'CORRECTED_DATA',
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lmeter[i] = f
          Rmeter[i] = list1[i]/60
        np.save("DATA/Lmeter%i"%(int(dtime*StepTime)),Lmeter)
        np.save("DATA/Rmeter%i"%(int(dtime*StepTime)),Rmeter)




#--------This function Move a source from the phase center to a maximun piosition that you precised and average over time-------

def makepsfloresflag (mcor=None,lcor="0", column="DATA",start_time=0,overlap_time=0,Ntime=1200,
						startfreq=0,overlap_freq=0,Nfreq=1,dtime=1200,dfreq=1):
	"""This function make the psf of a source at coordinates (lcor,mcor)"""
        if mcor==None:
          mcor = [0,30,60,90,120,150,180];
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
	i = 0
        for mcor1 in mcor:
          options = {}
          options['gridded_sky.grid_m0'] = mcor1
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS", config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);

	  mshi =MSUVDISTANCE.MSUVDISTANCE(hires+"/",column="DATA",time0=start_time,ntime=Ntime,freq0=startfreq,nfreq=Nfreq)
          arrays , time1,dtimemin=mshi.overlap_window(sinc1,1,1200,0)#window
          if i ==0:
                reset_ms(ORG_MS, time1, 1.*dtimemin,1, 125000.0)
                os.system("addbitflagcol MeerKAT-snapshot-21cm.MS")
          mshi.flaggindata(v.MS,arrays,time1,column="CORRECTED_DATA",irgg=i)
	  i = i + 1;
          # mshi =MSResampler.MSResampler(hires+"/",column="DATA",time0=start_time,ntime=Ntime,freq0=startfreq,nfreq=Nfreq)
          # arrays = mshi.boxcar (dtime,dfreq)
          # #arrays =  mshi.overlap_window (sinc1,dtime,dfreq,overlap_time=overlap_time,overlap_freq=overlap_freq,dump=None,dumpfig=None)
          # MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA") 

          center_min = -45*60+float(mcor1);
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
	  imager.make_image(column = "CORRECTED_DATA", phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = True, dirty = False, restore_lsm = False, weight = "natural");
          data = pyfits.open(imager.RESTORED_IMAGE)[0].data[0][0]
#          imager.make_image(column = "CORRECTED_DATA", phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
#                            restore = False, dirty = True, restore_lsm = False, weight = "natural");
#          data = pyfits.open(imager.DIRTY_IMAGE)[0].data[0][0]
          np.save("DATA/data%fdeg"%(mcor1/60.),data)
        # measure the psf size in arcsec
          absc = np.linspace(-len(data)/2,len(data)/2,len(data));
          absc = absc*float((imager.cellsize).split("arcsec")[0]);
          pylab.plot(absc,data[:,len(data)/2], label="src at:%fdeg"%(mcor1/60.))
        pylab.legend() 
        pylab.grid()
        pylab.savefig("DATA/psflores%i.png"%int(mcor1))


def Averaging_OverTime():
  # defined your parameters here

        Ntimeslots=1000 # Number of timeslots of the Hight resolution MS
        StepTime=1.  #StepTime of timeslots of the Hight resolution MS
        NFreq=100   # Number of Frequency of timeslots of the Hight resolution MS
	dfreq = 25
        StepFreq=1000000.0  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        dtime=4     # integration or the number of sample where the signal is cut-off
        NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
        stepR=4.      # the step radius in arcmin of the x-axis
        lenghtR=1800.1   # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, NsampleTime, StepTime, NFreq, StepFreq)
        #reset_ms(ORG_MS, 1, StepTime*dtime,1, StepFreq*dfreq)i
        list1 = np.arange(0,lenghtR,stepR)
        Lav = np.zeros(len(list1))
        Rav = np.zeros(len(list1))
        c = 0
        for i  in range(len(list1)):
          options = {}
          options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
	  center_min = -45*60+list1[i]
          #center_min = -1*45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          #mshi = ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime, NsampleTime, None)
          #dataav,ntim = mshi.NormalAvgALBL();
          #Test.SaveVis(v.MS,dataav,"CORRECTED_DATA");
	  reset_ms(ORG_MS, 150, StepTime*dtime,2, StepFreq*dfreq)
          os.system("addbitflagcol %s"%v.MS)
	  mshi =MSUVDISTANCE.MSUVDISTANCE(hires+"/",column="DATA",time0=200,ntime=600,freq0=25,nfreq=50)
	  arrays =mshi.boxcar(dfreq,dtime)
	  #arrays = mshi.overlap_window(sinc2,dtime,dfreq,200,25)#window
	  #if i ==0:dtimelpq,dfreq,overlap_time=0,overlap_freq=0,dump=None,dumpfig=None):
	 # reset_ms(ORG_MS, 120, StepTime*dtime,dfreq, StepFreq*dfreq)
	 # os.system("addbitflagcol VLAC-snapshot-21cm.MS")
	  #mshi.flaggindata(v.MS,arrays,time1,column="CORRECTED_DATA",irgg=i)
	  MSUVDISTANCE.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA",i=i)
          imager.make_image(column = "CORRECTED_DATA",
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lav[i] = f
          Rav[i] = list1[i]/60
        np.save("DATA/Lavbox2",Lav)
        np.save("DATA/Ravbox2",Rav)
	import pylab;
	#pylab.plot(Rav,Lav)
	#pylab.show()



#--------This function Move a source from the phase center to a maximun position that you precised and convolve over time-------
def Convolve_OverTime():
  # defined your parameters here

        Ntimeslots=100# Number of timeslots of the Hight resolution MS
        StepTime=0.1  #StepTime of timeslots of the Hight resolution MS
        NFreq=1#50    # Number of Frequency of timeslots of the Hight resolution MS
	dfreq = 1# 50
        StepFreq=125e3  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        Fov=2.       # Fov of the Baseline dependent filter
        dtime=100     # integration or the number of sample where the signal is cut-off
        NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
        Overlap=0     # the number of samples extends, Overlap<= Ntimeslots + dtime
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400    # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, NsampleTime, StepTime, NFreq, StepFreq)
        reset_ms(ORG_MS, (100)/dtime, StepTime*dtime, (1)/dfreq, StepFreq*dfreq)
        deg = (Fov*np.pi**2)/(180.)
        sinc1 = lambda x : (np.sin(x*deg))/(x*deg)
        list1 = np.arange(0,lenghtR,stepR)
        Lmeter = np.zeros(len(list1))
        Rmeter = np.zeros(len(list1))
        
        for i  in range(len(list1)):
          options = {}
          options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
          center_min = -45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          #mshi = ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime, NsampleTime, None, Overlap)
          #datatimeter,time = mshi.ConvolveALBL(sinc1,1)
          #Test.SaveVis(v.MS,datatimeter,"CORRECTED_DATA");
	  mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=0,ntime=100,freq0=0,nfreq=1)#50)
          arrays=mshi.window (sinc2,dtime,dfreq,dump=None,dumpfig=None)
          MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA")
          imager.make_image(column = 'CORRECTED_DATA',
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lmeter[i] = f
          Rmeter[i] = list1[i]/60
          np.save("DATA/Lmeter%i"%(int(dtime*StepTime)),Lmeter)
          np.save("DATA/Rmeter%i"%(int(dtime*StepTime)),Rmeter)

def Overlap_OverTime():
  # defined your parameters here

        Ntimeslots=600#500 # Number of timeslots of the Hight resolution MS
        StepTime=1.5  #StepTime of timeslots of the Hight resolution MS
        NFreq= 1     # Number of Frequency of timeslots of the Hight resolution MS
	dfreq = 1
        StepFreq=125e3  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        Fov=2.       # Fov of the Baseline dependent filter
        dtime=100     # integration or the number of sample where the signal is cut-off
        NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
        Overlap=0     # the number of samples extends, Overlap<= Ntimeslots + dtime
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400    # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, NsampleTime, StepTime, NFreq, StepFreq)
        reset_ms(ORG_MS, (Ntimeslots)/dtime, StepTime*dtime, (NFreq)/dfreq, StepFreq*dfreq)
        deg = (Fov*np.pi**2)/(180.)
        sinc1 = lambda x : (np.sin(x*deg))/(x*deg)
        list1 = np.arange(0,lenghtR,stepR)
        Lmeter = np.zeros(len(list1))
        Rmeter = np.zeros(len(list1))
        
        for i  in range(len(list1)):
          options = {}
         # options['gridded_sky.grid_m0'] = 60
	  options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
          center_min = -45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          #mshi = ClassMS1.ClassMS1(MSHI+"/","UVW","DATA",dtime)
         # datatimeter = mshi.ConvolveALBL(sinc1)
        #  Test.SaveVis(v.MS,datatimeter,"CORRECTED_DATA");
	  import scipy.signal
	  selepian = scipy.signal.slepian(dtime,.1)
	  mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=0,ntime=(Ntimeslots),freq0=0,nfreq=(NFreq))
	  #arrays = mshi.boxcar(dtime,dfreq)
          arrays=mshi.overlap_window (sinc1,dtime,dfreq,overlap_time=0,overlap_freq=0,dump=None,dumpfig=None)
          MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA")	
          imager.make_image(column = 'CORRECTED_DATA',
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lmeter[i] = f
          Rmeter[i] = list1[i]/60
        np.save("DATA/Lmeter%i"%(int(dtime*StepTime)),Lmeter)
        np.save("DATA/Rmeter%i"%(int(dtime*StepTime)),Rmeter)




#--------This function Move a source from the phase center to a maximun position that you precised and average over Frequency-------
def Averaging_OverFrequency():
  # defined your parameters here

        Ntimeslots=100 # Number of timeslots of the Hight resolution MS
        StepTime=0.1  #StepTime of timeslots of the Hight resolution MS
        NFreq=30     # Number of Frequency of timeslots of the Hight resolution MS
        StepFreq=1000  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        dtime=10    # integration or the number of sample where the signal is cut-off
        NsampleFreq=NFreq # default equal to Ntimeslots or less than 
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400     # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, Ntimeslots, StepTime, NsampleFreq, StepFreq)
        reset_ms(ORG_MS, Ntimeslots, StepTime,NsampleFreq/dtime, StepFreq*dtime)
        list1 = np.arange(0,lenghtR,stepR)
        Lav = np.zeros(len(list1))
        Rav = np.zeros(len(list1))
        
        for i  in range(len(list1)):
          options = {}
          options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
          center_min = -45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          mshi = ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime, None, NsampleFreq)
          dataav,ntim = mshi.NormalAvgALBL();
          Test.SaveVis(v.MS,dataav,"CORRECTED_DATA");
          imager.make_image(column = 'CORRECTED_DATA',
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lav[i] = f
          Rav[i] = list1[i]/60
          np.save("DATA/Lav%i"%(int(dtime*StepTime)),Lav)
          np.save("DATA/Rav%i"%(int(dtime*StepTime)),Rav)


#--------This function Move a source from the phase center to a maximun position that you precised and convolve over Frequency-------
def Convolve_OverFrequency():
  # defined your parameters here

        Ntimeslots=100 # Number of timeslots of the Hight resolution MS
        StepTime=0.1  #StepTime of timeslots of the Hight resolution MS
        NFreq=30     # Number of Frequency of timeslots of the Hight resolution MS
        StepFreq=10000  # Step Frequency of Frequency of timeslots of the Hight resolution MS
        Fov=2.       # Fov of the Baseline dependent filter
        dtime=10    # integration or the number of sample where the signal is cut-off
        NsampleFreq=NFreq # default equal to Ntimeslots or less than 
        Overlap=30     # the number of samples extends, Overlap<= Nfreq + dtime
        stepR=10.      # the step radius in arcmin of the x-axis
        lenghtR=400     # the lenght in arcmin of the x-axis 
          
        ORG_MS = v.MS
        hires = ORG_MS.split(".MS")[0]+"-hires.MS"
        reset_ms(hires, Ntimeslots, StepTime, NsampleFreq, StepFreq)
        reset_ms(ORG_MS, Ntimeslots, StepTime, NsampleFreq/dtime, StepFreq*dtime)
        deg = (Fov*np.pi**2)/(180.)
        sinc1 = lambda x : (np.sin(x*deg))/(x*deg)
        list1 = np.arange(0,lenghtR,stepR)
        Lmeter = np.zeros(len(list1))
        Rmeter = np.zeros(len(list1))
        
        for i  in range(len(list1)):
          options = {}
          options['gridded_sky.grid_m0'] = list1[i]
          options['ms_sel.msname'] = hires;
          mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
                  config="tdlconf.profiles",section="Sim_source_radius",
                  options=options);
          center_min = -45*60+list1[i];
          center_deg = math.ceil(center_min/60);
          center_min = abs(center_min - center_deg*60);
          mshi = ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime, None, NsampleFreq, Overlap)
          datatimeter,time = mshi.ConvolveALBL(sinc1,1)
          Test.SaveVis(v.MS,datatimeter,"CORRECTED_DATA");
          imager.make_image(column = 'CORRECTED_DATA',
                            phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
                            restore = False, dirty = True, restore_lsm = False,
                            weight = "natural");
          f,r,stm = Test.giveImage(imager.DIRTY_IMAGE)
          Lmeter[i] = f
          Rmeter[i] = list1[i]/60
          np.save("DATA/Lmeter%i"%(int(dtime*StepTime)),Lmeter)
          np.save("DATA/Rmeter%i"%(int(dtime*StepTime)),Rmeter)


MSHI_Template = "${MS:BASE}-hires.MS"


#----------Averaging over time and frequency of a Full MS---------------
def Averaging_OverTimeXFreq_FullMS():
  Ntimeslots=14400 # Number of timeslots of the Hight resolution MS
  NFreq=60   # Number of Frequency of timeslots of the Hight resolution MS
  dtime=10     # integration or the number of sample where the signal is cut-off
  NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
  NsampleFreq=NFreq # default  equal to Number of frequency hight resolution or less than
  mshi=ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime,NsampleTime,NsampleFreq)
  dataav,ntim=mshi.NormalAvgALBL();
  Test.SaveVis(v.MS,dataav,"DATA");


#----------Convolve over time and frequency of a Full MS---------------
def Convolve_OverTimeXFreq_FullMS():
  Ntimeslots=14400 # Number of timeslots of the Hight resolution MS
  NFreq=60     # Number of Frequency of timeslots of the Hight resolution MS
  dtime=60     # integration or the number of sample where the signal is cut-off
  NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
  NsampleFreq=NFreq # default  equal to Number of frequency hight resolution or less than
  FoV=2.    # Fov
  Overlap=0 # the number of samples extends, Overlap<= 2*Ntimeslots - dtime
  deg=(FoV*np.pi**2)/(180.)
  sinc1=lambda x : (np.sin(x*deg))/(x*deg)
  mshi=ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime,NsampleTime,NsampleFreq, Overlap)
  dataav, time= mshi.ConvolveALBL(sinc1,1)
  Test.SaveVis(v.MS,dataav,"DATA");

#-----------Averaging over Time of a full MS
def Averaging_OverTime_FullMS():
  # defined your parameters here
  Ntimeslots=120# Number of timeslots of the Hight resolution MS
  dtime=120     # integration or the number of sample where the signal is cut-off
  NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
  mshi=ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime,NsampleTime,None)
  dataav, time= mshi.NormalAvgALBL();
  Test.SaveVis(v.MS,dataav,"DATA");


#-----------Convolve over Time of a full MS
def Convolve_OverTime_FullMS():
  # defined your parameters here
  Ntimeslots=43200# Number of timeslots of the Hight resolution MS
  dtime=100   # integration or the number of sample where the signal is cut-off
  NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
  FoV=2.    # Fov
  Overlap=600 # the number of samples extends, Overlap<= 2*Ntimeslots - dtime
  deg=(FoV*np.pi**2)/(180.)
  sinc1=lambda x : (np.sin(x*deg))/(x*deg)
  mshi=ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime,NsampleTime,None,Overlap)
  dataav, time= mshi.ConvolveALBL(sinc1,1)
  Test.SaveVis(v.MS,dataav,"DATA");

#-----------Averaging over Frequency of a full MS
def Averaging_OverFrequency_FullMS():
  # defined your parameters here
  NFreq=300 # Number of timeslots of the Hight resolution MS
  dtime=100     # integration or the number of sample where the signal is cut-off
  NsampleFreq=Nfreq # default equal to Nunmber of frquency or less than 
  mshi=ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime,None,NsampleFreq)
  dataav, time= mshi.NormalAvgALBL();
  Test.SaveVis(v.MS,dataav,"DATA");

#-----------Convolve over Frequency of a full MS
def Convolve_OverFrequency_FullMS():
  # defined your parameters here
  NFreq=300 # Number of timeslots of the Hight resolution MS
  dtime=100     # integration or the number of sample where the signal is cut-off
  NsampleFreq=Nfreq # default equal to Nunmber of frquency or less than 
  FoV=2.    # the number of samples extends, Overlap<= 2*Nfreq - dtime
  Overlap=20 # Numbers of bins extended
  deg=(FoV*np.pi**2)/(180.)
  sinc1=lambda x : (np.sin(x*deg))/(x*deg)
  mshi=ClassMS.ClassMS(MSHI+"/","UVW","DATA",dtime,None,NsampleFreq,Overlap)
  dataav, time= mshi.ConvolveALBL(sinc1,1);
  Test.SaveVis(v.MS,dataav,"DATA");

MSHI_Template = "${MS:BASE}-hires.MS"
IMGNAME_Template = "plots-${MS:BASE}/${MS:BASE}-$IMGLABEL.dirty.fits"
imager.DIRTY_IMAGE_Template = "${OUTFILE}-${COUNTER}-${METHOD}.dirty.fits"
#run whole pipeline
def run_all():
  Ntimeslots=32400#500 # Number of timeslots of the Hight resolution MS 32400
  StepTime=1.  #StepTime of timeslots of the Hight resolution MS
  NFreq= 100     # Number of Frequency of timeslots of the Hight resolution MS
  dfreq = 50
  StepFreq=125e3  # Step Frequency of Frequency of timeslots of the Hight resolution MS
  Fov=2.       # Fov of the Baseline dependent filter
  dtime=90     # integration or the number of sample where the signal is cut-off
  NsampleTime=Ntimeslots # default equal to Ntimeslots or less than 
  Overlap=0     # the number of samples extends, Overlap<= Ntimeslots + dtime
  stepR=10.      # the step radius in arcmin of the x-axis
  lenghtR=400    # the lenght in arcmin of the x-axis 
  
  ORG_MS = v.MS
  hires = ORG_MS.split(".MS")[0]+"-hires.MS"
  #reset_ms(hires, NsampleTime, StepTime, NFreq, StepFreq)
  #reset_ms(ORG_MS, (Ntimeslots-3600)/dtime, StepTime*dtime, (NFreq-100)/dfreq, StepFreq*dfreq)
  #options = {}
  #options['gridded_sky.grid_m0'] = 360
  #options['gridded_sky.grid_l0'] = 360
  #options['ms_sel.msname'] = MSHI;
  #mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
  #        config="tdlconf.profiles",section="Sim_source_radius",
  #        options=options);
  #options = {}
  #options['gridded_sky.grid_m0'] = 10*60
  #options['ms_sel.msname'] = MSHI;
  #mqt.run("turbo-sim.py",job = "_tdl_job_1_simulate_MS",
   #              config="tdlconf.profiles",section="Sim_source_radius",
    #             options=options);
  #center_min = -45*60+360;
  #center_deg = math.ceil(center_min/60);
  #center_min = abs(center_min - center_deg*60);
  #import scipy.signal
  #selepian = scipy.signal.slepian(dtime,.1)
  mshi =MSResampler.MSResampler(MSHI+"/",column="DATA",time0=1800,ntime=(Ntimeslots-3600),freq0=25,nfreq=(NFreq-50))
  arrays = mshi.boxcar(dtime,dfreq)
  #arrays=mshi.overlap_window (sinc2,dtime,dfreq,overlap_time=0,overlap_freq=0,dump=None,dumpfig=None)
  MSResampler.save_visibility_arrays (v.MS,arrays,column="CORRECTED_DATA")
  #imager.make_image(column = 'CORRECTED_DATA',
  #                  phasecenter = "j2000,0h0m,%dd%dm"%(center_deg,center_min),
  #                  restore = False, dirty = True, restore_lsm = False,
  #                  weight = "natural");

 

  
    
