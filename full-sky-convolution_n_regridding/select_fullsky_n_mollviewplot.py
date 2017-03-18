import pyfits
import numpy as np
from matplotlib import pylab as plt
import healpy as hp
from matplotlib import cm
# a function to save plots
def save_fig(prefix):
	"Save current figure in extended postscript and PNG formats."
	plt.savefig('%s.png' % prefix, format='PNG', dpi = 200)
	plt.savefig('%s.eps' % prefix, format='EPS', dpi = 200)
#
#
def figur(fig_num, w_size, h_size, DPI):
    fx = plt.figure(fig_num, figsize=(w_size, h_size), dpi= DPI)
    
    return fx
    	
Tn = 1000   # freq term in MHz
a = 400	    # 1st term
diff = 2.2     # common difference


datamap = []
for iter in [1000, 1100, 1200, 1300, 1400]:	#[0, 49, 99, 149, 199, 249]
     
	hdu = pyfits.open('%s' %('full-sky-map_Jy-p-pixel_500-4-3145728_400-1450MHz.fits'))
	n = int(( iter - a)/diff + 1)  	# channel		     
	datamap.append(hdu[0].data[n]) 
	    
	#
pyfits.writeto('FullSky_map_N512_chanwidth_100MHz_1g_1.1g_1.2g_1.3g_1.4g.fits', np.array(datamap), clobber=True)
sky_map_data = np.array(datamap)#[:,:2,:]
		     
fx = plt.figure(num=1, figsize=(15, 10), dpi= 500) #facecolor = 'red',  edgecolor = 'white')	
#fx.patch.set_alpha(0.7)

# +++++++ sets background to white +++++++++++	  
get_cmap = cm.jet
get_cmap.set_under("w")    
for iter in xrange(5):
	if iter == 0:
		hp.visufunc.mollview(sky_map_data[iter,0,...], coord = ['G', 'C'], norm= 'log', xsize = 3200, unit = 'Jy/pixel', cbar=True,\
        	notext = True, margins=(0.01, 0.01, 0.0, 0.02), sub = (5, 3, 3*iter + 1), cmap = get_cmap)  #	plt.get_cmap('jet')
        	plt.title(r'$ I $', fontsize=20, fontweight="bold") 
		#
		hp.visufunc.mollview(sky_map_data[iter,1, ...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, unit = 'Jy/pixel', cbar=True,\
		notext = True, margins=(0.01, 0.01, 0.0, 0.02), sub = (5, 3, 3*iter + 2), cmap = get_cmap)
		plt.title(r'$ Q $', fontsize=20, fontweight="bold") 
			 
	 	hp.visufunc.mollview(sky_map_data[iter,2, ...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, unit = 'Jy/pixel', cbar=True,\
        	notext = True, margins=(0.01, 0.01, 0.0, 0.02), sub = (5, 3, 3*iter + 3), cmap = get_cmap)
        	plt.title(r'$ U $', fontsize=20, fontweight="bold") 
        	 #
	else:
		hp.visufunc.mollview(sky_map_data[iter,0,...], coord = ['G', 'C'], norm= 'log', xsize = 3200, unit = 'Jy/pixel', cbar=True,\
			 notext = True, title = ' ', margins=(0.01, 0.01, 0.0, 0.02), sub = (5, 3, 3*iter + 1), cmap = get_cmap)
		#
		hp.visufunc.mollview(sky_map_data[iter,1, ...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, unit = 'Jy/pixel', cbar=True,\
			 notext = True, title = ' ',  margins=(0.01, 0.01, 0.0, 0.02), sub = (5, 3, 3*iter + 2), cmap = get_cmap)
			 
	 	hp.visufunc.mollview(sky_map_data[iter,2, ...], coord = ['G', 'C'], norm= 'linear', xsize = 3200, unit = 'Jy/pixel', cbar=True,\
			 notext = True, title = ' ', margins=(0.01, 0.01, 0.0, 0.02), sub = (5, 3, 3*iter + 3), cmap = get_cmap)
#fx.patch.set_facecolor('black')
fx.savefig('stokes-plot')
#plt.show()
plt.close(1)			

#save_fig('stokes-plot')
#plt.close()
