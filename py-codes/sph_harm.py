#!/usr/bin/env python
from utilities import *
import scipy.special as sp

class SpheHarm(object):
    def __init__(self, data=[], lmax=10, coeffs=None, deco=True, reco=False, fits_file=''):
        if coeffs!=None: self.coeffs = coeffs
        try: data = fits.getdata(fits_file)
        except: None
        self.deco = deco
        self.reco = reco
        pi = np.pi
        Npix = data.shape[-1]*1j
        phi, theta = np.mgrid[0:2*pi:Npix, 0:pi:Npix]
        self.L, self.M = [], []
        for l in range(lmax+1):
            for m in range(-l,l+1):
                self.L.append(l)
                self.M.append(m)
        self.Nmodes = len(self.L)
        self.bases = [self.sphe_harm_j(j, phi, theta) for j in range(self.Nmodes)]
        self.lmax = lmax

        #print len(data.shape)
        #return
        
        if len(data.shape)==4:
            self.coeffs, self.recons = self.jones_images(data)
        elif len(data.shape)>=5:
            self.all_freqs(fits_file)
            
    def jones_images(self, data):
        coeffs = np.zeros((data.shape[0], data.shape[1], self.Nmodes), dtype=np.complex)
        recons = np.zeros(data.shape, dtype=np.complex)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if self.deco:
                    print "Fitting spherical harmonics to %i %i with maximum l = %i"%(i,j,self.lmax)
                    self.decompose(data[i,j,...])
                    coeffs[i,j,...] = self.coeffs
                if self.reco:
                    print "Reconstructing using the coeffs"
                    self.reconstruct(self.coeffs)
                    recons[i,j,...] = self.recon
        return coeffs, recons

    def all_freqs(self, fits_file):
        start = time.time()
        d = fits.getdata(fits_file)
        h = fits.getheader(fits_file)
        dc = d[0,...] + 1j * d[1,...]
        dc = np.nan_to_num(dc)
        print dc.shape
        recons = np.zeros(dc.shape, dtype=np.complex)
        coeffs = np.zeros((2,2,dc.shape[2],self.Nmodes), dtype=np.complex)
        for f in range(dc.shape[2]):
            print '... Channel %i'%(f+1)
            if self.deco:
                coeffs[:,:,f,:] = self.jones_images(dc[:,:,f,:,:])[0]
            if self.reco:
                coeffs[:,:,f,:,:], recons[:,:,f,:,:] = self.jones_images(dc[:,:,f,:,:])

        np.save(fits_file[:-5]+'_sh_%icoeffs.npy'%self.Nmodes, coeffs)

        if self.reco:
            d[0,...], d[1,...] = recons.real, recons.imag
            fits.writeto(fits_file[:-5]+'_sh_recon_%icoeffs.fits'%self.Nmodes, circular_mask(d), h, overwrite=True)

        end = time.time()
        print "... Time taken: %.2f minutes"%((end-start)/60.)
        
    def sphe_harm_j(self, j, phi, theta):
        l, m = self.L[j], self.M[j]
        return sp.sph_harm(m, l, phi, theta)
    
    def decompose(self, img):
        bases = self.bases
        cov_mat = np.array([[np.sum(i * j) for i in bases] for j in bases])
        cov_mat_in = np.linalg.pinv(cov_mat)
        innerprod = np.array([np.sum(img * basis) for basis in bases])
        self.coeffs = np.dot(cov_mat_in, innerprod)

    def reconstruct(self, coeffs):
        self.recon = np.sum(val * self.bases[i] for (i, val) in enumerate(coeffs))

if __name__=='__main__':
    mod = SpheHarm(fits_file=argv[1], lmax=int(argv[2]))

