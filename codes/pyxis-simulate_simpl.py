
import numpy as np
from pyrap.tables import table
from Pyxis.ModSupport import *
import mqt,cal,lsm,imager
from lofar import bdsm
import os,sys


v.DESTDIR_Template = "${OUTDIR>/}plots-${MS:BASE}"
v.OUTFILE_Template = "${DESTDIR>/}${MS:BASE}"
v.LSM_Template = '${v.OUTFILE}.lsm.txt'
v.LOG_Template = "${v.OUTFILE}.log"
# number of meq threads
mqt.MULTITHREAD=4

# imaging options
imager.npix = 1024
imager.cellsize = "3arcmin"
imager.mode = "channel"
imager.stokes = "I"
imager.weight = "uniform"
imager.wprojplanes = 0

imager.niter = 100000
imager.gain = .1
imager.threshold ="0.1mJy"

v.MS = 'zen.2455819.67684.uvcRREM.MS'

def simulate(msname="$MS",skymodel=None,section="$SIM_CONFIG_SECTION",args=[],options={},**kw):
     """
        Simulating our skymodel
     """
     # necessary to make '$MS' etc. above work properly
     msname,section = interpolate_locals("msname section");
     options.update(kw);
     # setup TDL options for MS, LSM, etc., borrowing templates from the cal module
     args = [ "${cal.MS_TDL} ${cal.CHAN_TDL} ${cal.LSM_TDL}" ] + list(args);

     # run meqtrees
     options['tiggerlsm.filename']=skymodel
     options['ms_sel.msname']=msname
     mqt.run("turbo-sim.py",job="_tdl_job_1_simulate_MS",config="tdlconf.profiles",section=section,options=options)
     imager.make_image(column='CORRECTED_DATA',dirty=False,restore=True,restore_lsm=False,algorithm='csclean')
     x.sh('mv %s %s'%(imager.RESTORED_IMAGE,'Image.fits')) 

def getData(msname="$MS",outfile=None,imagename=None):
   msname,section = interpolate_locals("msname section")
   ta = table(msname,readonly=False)
   data = ta.getcol('DATA')
   corr_data = ta.getcol('CORRECTED_DATA')
   new_data = data + corr_data
   ta.putcol('CORRECTED_DATA',new_data)
   ta.close()
   # Imaging
   imager.make_image(column='CORRECTED_DATA',dirty=False,restore=True,restore_lsm=False,algorithm='csclean')
   #extractflux(imager.RESTORED_IMAGE,outfile)
   x.sh('mv %s %s'%(imager.RESTORED_IMAGE,imagename))

def generate_pos(fov=None, num_sources=None):
    y=np.random.uniform(low=-1*np.absolute(fov),high=np.absolute(fov),size=num_sources)*(np.pi/180)
    return y

def get_field_center(MS):
    t = table(MS+"/FIELD")
    phase_centre = (t.getcol("PHASE_DIR"))[0,0,:]
    t.close()
    return phase_centre[0], phase_centre[1] #ra0,dec0 in radians

def lm2radec(MS,l,m):#l and m in radians
    rad2deg = lambda val: val * 180./np.pi
    ra0,dec0 = get_field_center(MS) # phase centre in radians
    rho = np.sqrt(l**2+m**2)
    if rho==0:
       ra = ra0
       dec = dec0
    else:
       cc = np.arcsin(rho)
       ra = ra0 - np.arctan2(l*np.sin(cc), rho*np.cos(dec0)*np.cos(cc)-m*np.sin(dec0)*np.sin(cc))
       dec = np.arcsin(np.cos(cc)*np.sin(dec0) + m*np.sin(cc)*np.cos(dec0)/rho)
    return rad2deg(ra), rad2deg(dec)

def radec2lm(MS,ra_d,dec_d):# ra and dec in degrees
    rad2deg = lambda val: val * 180./np.pi
    deg2rad = lambda val: val * np.pi/180
    ra0,dec0 = get_field_center() # phase centre in radians
    ra_r, dec_r = deg2rad(ra_d), deg2rad(dec_d) # coordinates of the sources in radians
    l = np.cos(dec_r)* np.sin(ra_r - ra0)
    m = np.sin(dec_r)*np.cos(dec0) - np.cos(dec_r)*np.sin(dec0)*np.cos(ra_r-ra0)
    return rad2deg(l),rad2deg(m)

def meqskymodel(MS,point_sources,outfile):
     str_out = "#format: name ra_d dec_d i\n"
     for i in range(len(point_sources)):
          amp, l ,m = point_sources[i,0], point_sources[i,1], point_sources[i,2]
          ra_d, dec_d = lm2radec(MS,l,m)
          name = "A"+ str(i)
          str_out += "%s %.12g %.12g %.5g\n"%(name, ra_d, dec_d,amp)

     skymodel = open(outfile,"w")
     skymodel.write(str_out)
     return outfile

def GeneratePointSources(MS, num_sources, fov):
    #contain all the pointsources
    point_sources = np.zeros((num_sources,3))

    #generate positions
    point_sources[:,0] = 1
    point_sources[:,1] = generate_pos(fov=fov, num_sources=num_sources)
    point_sources[:,2] = generate_pos(fov=fov, num_sources=num_sources)

    return point_sources

def GenerateRandomSkies(MS,num_sources,fov,outfile):
    """
      Generating the point sources
    """
    point_sources = GeneratePointSources(MS, num_sources, fov)
    meqskymodel(MS,point_sources,outfile)

def extractflux(image,outfile):
  pybdsmDict={            # define a minimum dictionary for source extraction
   'rms_map': None,
   'thresh': 'hard',
   'thresh_isl': 3.,
   'thresh_pix': 7.0}
 
  img = bdsm.process_image(image,**pybdsmDict)                # run pybdsm
  img.write_catalog(outfile=outfile,format='ascii',catalog_type='gaul',clobber=True) # output the source catalogue
 

def runall():
      skymodel = "skymodel.txt"
      GenerateRandomSkies(v.MS,1,2,skymodel) 
      simulate("$MS",skymodel,section="sim")
      
    
