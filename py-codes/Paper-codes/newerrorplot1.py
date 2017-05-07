#
from matplotlib import pyplot as plt
import numpy as np
import healpy as hp
from scipy.stats import sem
from scipy.stats import mstats
from scipy.signal import savgol_filter 
#


def figur(n1, n2, n3, n4):
	    fx = plt.figure(n1, figsize=(n2, n3), dpi= n4)
	    
	    return fx
	    
# Load Error Data
ck = np.load('Tconv_map.npy')
ch = np.load('H5conv_map.npy')
xy_er = np.load('fdhase_error_2.npy')
gp_er = np.load('gphase_error_2.npy')
hh_er = np.load('h56_error_2.npy')
#lmx = hh_er = np.load('h56_error_0.npy')

# 	RMSE
#nm = int(700)
#gb = hp.sphtfunc.anafast(gp_er, lmax = nm) 
#xb = hp.sphtfunc.anafast(xy_er, lmax = nm) 
#hb = hp.sphtfunc.anafast(hh_er, lmax = nm)
#xg = np.sqrt(((gb - hb ) ** 2).mean())
#xf = np.sqrt(((xbr - hb ) ** 2).mean())
#xc = np.sqrt(((hc['IQU'][2]) ** 2).mean())
#xb = np.sqrt(((hb['IQU'][2]) ** 2).mean())
k = ['I', 'Q', 'U']*3
k1 =  sorted(k)



fig = figur(7, 12, 12, 80)
hg = []
hx = []
for i in range(9):
	
	nm = int(400)
	mp = hp.sphtfunc.anafast(ck[i], lmax = nm)
	h5p = hp.sphtfunc.anafast(ch[i], lmax = nm)  
	gb = hp.sphtfunc.anafast(gp_er[i], lmax = nm) 
	xb = hp.sphtfunc.anafast(xy_er[i], lmax = nm) 
	hb = hp.sphtfunc.anafast(hh_er[i], lmax = nm)
	#print np.sqrt(((xb) ** 2).mean())/abs(mp).mean()  #)*100
	#print sem(hb**2/(hb**2).mean())
	#print 
	#print gp_er[i].mean() #conv_map[i].mean()
	#print (gp_er[i].mean())/conv_map[i].mean()*100
	#print sem(gp_er[i]/conv_map[i].mean())
	cc = (gb)**2
	#xg = sem(cc/cc.max())
#	print xg.mean()/mp.mean() * 100
	xg = sem(cc/cc.mean(), axis=None, ddof=0)  #np.sqrt((((gb - hb)/hb ) ** 2).mean())*100
	#xg = sem(cc/mstats.gmean(cc), axis=None, ddof=0)
	cc = (xb)**2
	xf = sem(cc/cc.mean(), axis=None, ddof=0) #np.sqrt((((xb - hb)/hb ) ** 2).mean())*100
	#xf = sem(cc/mstats.gmean(cc), axis=None, ddof=0)
	#print xf*100
	cc = (hb)**2
	xh = sem(cc/cc.mean(), axis=None, ddof=0)
	hg.append(abs(xg-xh))
	hx.append(abs(xf-xh))
	ax = plt.subplot(3,3,i+1)	
	lmax = np.arange(len(gb))/3
	
	if i in [1,2,5,7]:
		mp = savgol_filter(mp, 21, 1, deriv = 0)
		h5p = savgol_filter(h5p, 21, 1, deriv = 0)  
		ax.semilogy(lmax,lmax*(lmax +1.)*mp/2.0*np.pi , label = '$Conv. Pow. Spec^{oskar}$' )
	
		ax.semilogy(lmax,lmax*(lmax +1.)*h5p/2.0*np.pi , label = '$Conv. Pow. Spec^{holo}$' )
	else:
		ax.semilogy(lmax,lmax*(lmax +1.)*mp/2.0*np.pi , label = '$Conv. Pow. Spec^{oskar}$' )
	
		ax.semilogy(lmax,lmax*(lmax +1.)*h5p/2.0*np.pi , label = '$Conv. Pow. Spec^{holo}$' )
	
	ax.semilogy(lmax,lmax*(lmax +1.)*gb/2.0*np.pi , label = '$SE^{GP} = %s $ ' %(round(xg, 4)) + str('%') )
	
	ax.semilogy(lmax,lmax*(lmax +1.)*xb/2.0*np.pi , label = '$SE^{XY} = %s $ ' %(round(xf, 4)) + str('%'))
	
	ax.semilogy(lmax,lmax*(lmax +1.)*hb/2.0*np.pi , label = '$SE^{HB} = %s $ ' %(round(xh, 4))+ str('%'))
	plt.legend(loc='upper right', framealpha = 0.5)
	plt.title(r'$ L_{se}^{%s \longrightarrow %s}$' %(k[i], k1[i]), fontsize=20, fontweight="bold")
fig.text(0.5, 0.04, r" multipole moment $[\, l\, ]$ ", ha='center', fontsize="20")
fig.text(0.04, 0.5, r" Intensity $[Jy]$", va='center', rotation='vertical', fontsize="20")
plt.show() 
np.savetxt('intro-err.txt', np.array([hg,hx]).T, delimiter=' ', fmt='%.4f')  #fmt='%s'
