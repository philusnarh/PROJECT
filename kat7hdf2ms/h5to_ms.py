import numpy as np
import re
import katdal
from pyrap.tables import table
from optparse import OptionParser
import time


# ms = '/home/narh/Measurement_Sets/KAT7_make300_4hr_1.3GHz.MS'
# fil = '1444074086.h5'

def kath5toms(h5file, msname):
	hdf5 = katdal.open('%s' %h5file)		
	corr_products = hdf5.corr_products
	

	# ++++++++++++ Select number from  string list +++++++++++ 

	s1 = [ map(int, re.findall('\d+', x))  for x in corr_products[:,0] ]
	s2 = [ map(int, re.findall('\d+', x))  for x in corr_products[:,1] ]

	#for iter in 
	#s = map(int, re.findall('\d+', an[:,0][0]))
	#print np.array(s)[0].shape

	#print len(np.squeeze(np.array(s)))
	corr_index1 = np.squeeze(np.array(s1))[::4]
	corr_index2 = np.squeeze(np.array(s2))[::4]


	# Extract visibilities
	#
	vis_data = hdf5.vis.dataset.value
	shape = vis_data.shape

	vis_corr =  np.array([vis_data[:,:, 0,0] + 1j*vis_data[:,:, 0,1], vis_data[:,:, 2,0] + 1j*vis_data[:,:, 2,1],
		              vis_data[:,:, 3,0] + 1j*vis_data[:,:, 3,1], vis_data[:,:, 1,0] + 1j*vis_data[:,:, 1,1]])

	vis_corr = np.swapaxes(vis_corr, 0, 2)
	vis_corr = np.swapaxes(vis_corr, 0, 1)
	
	antnum = np.unique(np.array(s2))
	tab = table("%s" %msname, readonly=False)
	tab.unlock()
	A0 = tab.getcol("ANTENNA1")
	A1 = tab.getcol("ANTENNA2") 
	ant1 = A0.copy() 
	ant2 = A1.copy() 
	
	for iter in xrange(len(antnum)):
		ant1[(A0==iter)] = antnum[iter]
		ant2[(A1==iter)] = antnum[iter]
		
	
	# ++++++++++++ Return index of a sorted list +++++++++++
	#vals = np.array(ant2)
	#sort_index = np.argsort(vals)
	#tab.putcol("ANTENNA1",ant1[sort_index])
	#tab.putcol("ANTENNA2",ant2[sort_index])

	tab.putcol("ANTENNA1",ant1)
	tab.putcol("ANTENNA2",ant2)
	data = tab.getcol("DATA")
	Nb = len(antnum)*(len(antnum) -1)/2 + len(antnum)
	i = 0 
	j = shape[0]/Nb
	 
	for iter in xrange(len(corr_index1)):
	
		data[(ant1 == corr_index1[iter])&(ant2 == corr_index2[iter])] = vis_corr[i:j,:, :] 
		print '\n i', i
		i += shape[0]/Nb
		j += shape[0]/Nb
	
	tab.putcol("DATA", data)
	tab.putcol("MODEL_DATA", data)
	tab.putcol("CORRECTED_DATA", data)	
	tab.close()
	
	return None
	

def run_sim_model():
	#
	#	start-time
	start = time.time()
	startime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Start at %s" % startime
	
	option = OptionParser(usage="usage: %prog [options] filename",
                              version="%prog 1.0")
	option.set_description('Converts kat7 rasterscan from HDF to ms')
	option.add_option('-H','--hdf', dest='h5file', 
	                  default='kat_hdf.h5', action="store", type="string",
	                  help='HDF filename')
	                  
  	option.add_option('-m','--msname', dest='msfile', 
	                  default='kat7msnam.ms', action="store", type="string",
	                  help='Measurement Set name')
	
	options,args = option.parse_args()
	
	kath5toms(h5file = options.h5file, msname = options.msfile)
	
	stoptime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
	print "Stop at %s" % stoptime 
	end = time.time()
	elasped_time = (end - start)/3600.0
	print "Total run time: %7.2f hours" % elasped_time
	

if __name__=="__main__":
	#
	run_sim_model()

#pyxis uv_plot
#pyxis MS=hirax.MS  ms.plot_uvcov
#pyxis MS=MeerKAT-snapshot-21cm.MS ms.plot_uvcov 



#pyxis Imaging
#pyxis hirax.MS im.cellsize=7arcsec im.npix=2048 im.stokes=I im.IMAGER=wsclean im.make_image[psf=True,column=DATA]


#add imaging column
#pyxis ms.add_imaging_columns 'hirax.MS'


