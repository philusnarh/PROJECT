import pylab, random, os, sys
import healpy as hp
import numpy as np
from astropy import units as u
from astropy import constants as c
from scipy import integrate
from bm_prms import prms


freqs = np.linspace(0.117,0.182,num=131) #aipy likes GHz units. avoiding band edges
nside = 128
npix = hp.nside2npix(nside)
#calculate relevant map parameters
c = 3e8 #m/s
#nside = hp.npix2nside(beam.shape[0])
#npix = beam.shape[0]
ipix = np.arange(npix)
theta,phi = hp.pix2ang(nside,ipix)

#we care about degree scales ~21 degrees
#lmax=9
lmax=3*nside - 1
l,m = hp.Alm.getlm(lmax)

print l
print m
#frequencies
nfreq=freqs.shape[0]
nu = np.outer(np.linspace(100e6,200e6,num=nfreq),np.ones(npix))#*u.Hz


########Verifying Ridhima sims
sky = hp.read_map('stokesI-f100_j2455819.54472.fits')


#hp.orthview(sky)
#pylab.show()

#promote sky to matrix for frequency axis
sky = np.outer(np.ones(nfreq),sky)*pow(nu/150e6,-0.7)

#decompose sky into alm's
n_alm = len(m)
alm = np.zeros((nfreq,n_alm),dtype='complex128')
print 'Calculating sky a_lm values:'
for i in range(nfreq):
	print nu[i,0]/1e6,'MHz'
	alm[i,:] = hp.map2alm(sky[i,:],lmax=lmax,iter=3)

#calculate fringe factor (true for all freqs)	
s  = np.array(hp.pix2vec(nside,ipix))

#
#
bl_length = 100.
#
#

#b = np.resize(np.repeat(np.array([0,bl_length,0]),npix),[3,npix])#*u.meter

#Checking Ridhima's simulations
b = np.resize(np.repeat(np.array([27.322,-239.994,0]),npix),[3,npix])#*u.meter


b_dot_s = np.sum(b*s,axis=0)
factor = np.exp(1.j*np.outer(np.ones(nfreq),b_dot_s)*nu/c) #c.c didn't work

#phasing
#rot_ang = np.linspace(0,5.0592,num=1632)<--LST bin rate
rot_ang = np.linspace(-np.pi,np.pi,num=360*4) #<--24 hours
n = len(rot_ang)

vis = np.zeros([n,nfreq],dtype='complex128')

print 'Calculating visibilities. Time stamp:'
for i in range(n):
	print i
	rotation = np.outer(np.ones(nfreq),np.exp(-1.j*m*rot_ang[i]))
	vis[i,:] = np.sum(alm*rotation,axis=1)

savefile = sys.argv[1]
print 'Saving visibility to %s...'%savefile
np.savez(savefile,vis=vis)	
