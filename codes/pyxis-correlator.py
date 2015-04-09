import os
import sys
import imager
import pyfits
Radian2Arcmin = lambda x:x*180*60/np.pi;
Radian2min = lambda x:x*12*60/np.pi;
imager.DIRTY_IMAGE_Template = "plots/dirty.fits"
imager.RESTORED_IMAGE_Template = "plots/restored.fits"
imager.MODEL_IMAGE_Template ="plots/model.fits"

def psf (msname = None):
        """Gets  the psf.
                Return an array of (npix,npix) 
        """;
        msname = msname or v.MS
        imager.make_image(column = "psf", weight = "natural");
 	data = pyfits.open(imager.DIRTY_IMAGE_Template)[0].data[0][0]
	np.save("DATA/SINC/w4800s/datapsf",data)
	return data

def Imagingperbaseline(column = "DATA", baseline = "0:1", center_deg = -49, center_min = 0, mcor1=0):
      """"
        Imaging along the specify baseline
      """

      if baseline == "all":
         antenna_str = ""
      else:
          antennas = baseline.split(':') or baseline.split('-')
          antenna_str ="(ANTENNA1 = "+antennas[0]+" && ANTENNA2 = "+antennas[1]+")"
      imager.npix= 256
      imager.cellsize = "2arcsec"
      dirty = dict(select='%s'%(antenna_str))
      restored = dict(select='%s'%(antenna_str))
      imager.make_image(column = "psf", dirty = dirty, weight = "natural");
      psf = pyfits.open(imager.DIRTY_IMAGE_Template)[0].data[0][0]

      np.save("DATA/psf%s%s"%(antennas[0],antennas[1]),psf)

      #imager.make_image(column='%s'%(column), dirty=dirty, restore=False, restore_lsm=False, algorithm='csclean')

def run_all():
	from pyrap.tables import table
	t = table("KAT7_.MS", readonly=False)
	UVW = t.getcol("UVW")
	A0 = t.getcol("ANTENNA1")
	A1 = t.getcol("ANTENNA2")
	na = np.max(A0)+1
	t.close()
	for p in range(na):
		for q in range(p+1,na):
			Imagingperbaseline(column = "DATA", baseline = "%i:%i"%(p,q))
