"""
Created by Asad,KMB (khmbasad@gmail.com) on 8 May 2017

If you want reconstruct beams of all requencies using Zernike polynbomials:
    python zernike.py allfreq filename:str Nmodes:int threshold:int

"""

#import opticspy
from utilities import *
from matplotlib import ticker
from scipy.misc import factorial as fac
#from holography import load_data

class Zernike(object):
    """
    Decompose and reconstruct a 2d image or 4d Jones matrix using Zernike polynomials
    author: Julein N. Girard (jgirard@ska.ac.za)
    """
    def __init__(self, data=[], Nmodes=50, threshold=None, Npix=None, m=0, n=0, recon=False):
        """
        Inputs
        -----------------
        data : numpy array
            2d image or 4d Jones matrix
        nmodes : int
            Number of modes, i. e. the maximum Noll index to use for decomposition and reconstruction
        threshold : float
            Percentage of highest energy coefficients to keep during truncation

        Outputs
        ----------
        coeffs
        coeffs_trunc
        recon_full
        recon_trunc
        res_full
        res_trunc
        """
        
        self.Nmodes = Nmodes
#        self.threshold = threshold

        # If 'data' is a single 2d image
        try:
            if len(data.shape)==2:
            	self.threshold = threshold
                self.img = data
                self.reconstruct()
            
            # If 'data' is a 2x2 Jones matrix
            elif len(data.shape)==4:
                self.coeffs_J = self.coeffs_trunc_J = np.zeros((data.shape[0], data.shape[1], Nmodes), dtype=np.complex)
                self.recon_full_J = self.recon_trunc_J = np.zeros(data.shape, dtype=np.complex)
                for i in range(data.shape[0]):
                    for j in range(data.shape[1]):
                        print "Fitting Zernike polynomials to Jones %i %i using %i modes"%(i,j,Nmodes)
                        self.img = data[i,j,:,:]
                        if recon:
                            if len(threshold) == 2: 
                            	if i != j: self.threshold = threshold[1]
                            	else: self.threshold = threshold[0]
                            self.reconstruct()
                            self.coeffs_J[i,j,:] = self.coeffs
                            self.coeffs_trunc_J[i,j,:] = self.coeffs_trunc
                            self.recon_full_J[i,j,:,:] = self.recon_full
                            self.recon_trunc_J[i,j,:,:] = self.recon_trunc
                        else:
                            self.reconstruct()
                            self.decompose()
                            self.coeffs_J[i,j,:] = self.coeffs

        except:
            if not Nmodes:
                self.unit_disk(npix=Npix)
                self.basis = self.zernike(m, n, self.grid_rho, self.grid_phi)*self.grid_mask
                self.basis = circular_mask(np.array(self.basis))
            if Nmodes:
                self.unit_disk(npix=Npix)
                self.basis = [self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for i in range(Nmodes)]
                self.basis = circular_mask(np.array(self.basis))

    def zernike_rad(self, m, n, rho):
        """
        Calculate the radial component of Zernike polynomial (m, n) 
        given a grid of radial coordinates rho.
        """
        if (n < 0 or m < 0 or abs(m) > n):
            raise ValueError
        if ((n-m) % 2):
            return rho*0.0
        pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
        return sum(pre_fac(k) * rho**(n-2.0*k) for k in xrange((n-m)/2+1))

    def zernike(self, m, n, rho, phi):
        """
        Calculate Zernike polynomial (m, n) given a grid of radial
        coordinates rho and azimuthal coordinates phi.
        """
        if (m > 0): return self.zernike_rad(m, n, rho) * np.cos(m * phi)
        if (m < 0): return self.zernike_rad(-m, n, rho) * np.sin(-m * phi)
        return self.zernike_rad(0, n, rho)

    def noll_to_zern(self, j):
        """
        Convert linear Noll index to tuple of Zernike indices.
        j is the linear Noll coordinate, n is the radial Zernike index and m is the azimuthal Zernike index.
        @param [in] j Zernike mode Noll index
        @return (n, m) tuple of Zernike indices
        @see <https://oeis.org/A176988>.
        """
        # add 1 to start from 1
        j += 1

        n = 0
        j1 = j-1
        while (j1 > n):
            n += 1
            j1 -= n

        m = (-1)**j * ((n % 2) + 2 * int((j1+((n+1)%2)) / 2.0 ))
        return (n, m)

    def zernikel(self, j, rho, phi):
        """
        Calculate Zernike polynomial with Noll coordinate j given a grid of radial
        coordinates rho and azimuthal coordinates phi.
        """
        nm = self.noll_to_zern(j)
        m, n = nm[1], nm[0]
        return self.zernike(m, n, rho, phi)
    
    def unit_disk(self, npix=None):
        """Create an unit disk and convert to rho, phi"""
        if npix: nx, ny = npix, npix
        else: nx, ny = self.img.shape
        grid = (np.indices((nx, ny), dtype=np.float) - nx/2) / (nx*1./2) # create unit grid [-1,1]
        self.grid_rho = (grid**2.0).sum(0)**0.5 # rho = sqrt(x^2+y^2)
        self.grid_phi = np.arctan2(grid[0], grid[1]) # phi = itan(x/y)
        self.grid_mask = self.grid_rho <= 1 # boolean array specifying where rho<=1
    
    def decompose(self):
        """Decompose using SVD"""
        self.unit_disk()
        # Caculate Zernike bases given the maximum Noll index, N
        N = self.Nmodes
        basis = [self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for i in range(N)]

        # Calculate covariance between all Zernike polynomials
        self.cov_mat = np.array([[np.sum(zerni * zernj) for zerni in basis] for zernj in basis])

        # Invert covariance matrix using SVD (  A x = b  ==>  x =>  x= A^{pseudoinv} b)
        self.cov_mat_in = np.linalg.pinv(self.cov_mat)

        # Inner product between the img and the Zernike bases
        self.innerprod = np.array([np.sum(self.img * zerni) for zerni in basis])
        
        # Dot product between inverse covariance matrix and the innerprod to get the coeffs
        self.coeffs = np.dot(self.cov_mat_in, self.innerprod)

    def truncate(self):
        """Truncate the coefficients upto the given threshold"""
        sortedindex = np.argsort(np.abs(self.coeffs))[::-1]
        Ncoeff = self.coeffs.shape[-1]
        cutoff = np.int(np.round(Ncoeff*self.threshold/100.))
        
        print "Keeping %2.0f %% (N=%s) of the biggest coefficients"%(self.threshold,cutoff)

        self.coeffs_trunc = self.coeffs.copy() # copy of all coeff
        self.coeffs_trunc[sortedindex[cutoff:]] = 0 # put coeff below threshold to 0

    def reconstruct(self):
        """Reconstruct a model image from the coeffcicients"""
        self.decompose()
        self.recon_full = np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs))
        self.res_full = (abs(self.img) - abs(self.recon_full))  * self.grid_mask
        
        # If a threshold is given, truncate the coeffs and create truncated reconstructions
        if self.threshold:
            self.truncate()
            self.recon_trunc = np.sum(val * self.zernikel(i, self.grid_rho, self.grid_phi)*self.grid_mask for (i, val) in enumerate(self.coeffs_trunc))
            self.res_trunc = (abs(self.img) - abs(self.recon_trunc)) * self.grid_mask
            self.diff_full_trunc = (self.recon_full - self.recon_trunc) * self.grid_mask

def zernike_all_freqs(fits_file, Nmodes=50, threshold=20, recon=False):
    """ Use the Zernike class to decompose multi-frequency beam data
    Input:
        fits_file: str, filename, data shape must be (2,Nchan,2,2,Npix,Npix)
        Nmodes: int, maximum Noll index
        threshold: int, percent of zern coeffs to keep during truncation
    Output:
        fits_file+'_zern.fits': Zernike reconstructions as a fits file
        fits_file+'_zern_coeffs.npy': Zernike coeffs
    """
    start = time.time()
    d = fits.getdata(fits_file)
    h = fits.getheader(fits_file)
    dc = d[0,...] + 1j * d[1,...]
    dc = np.nan_to_num(dc)
    full = np.zeros(dc.shape, dtype=np.complex)
    trunc = full
    coeffs = np.zeros((dc.shape[0],2,2,Nmodes), dtype=np.complex)
    for f in range(dc.shape[0]):
        print '... Channel %i'%(f+1)
        mod = Zernike(dc[f,...], Nmodes, threshold, recon=recon)
        coeffs[f,...] = mod.coeffs_J
        if recon:
            full[f,...] = mod.recon_full_J
            trunc[f,...] = mod.recon_trunc_J

    np.save(fits_file[:-5]+'_zern_coeffs_%i.npy'%Nmodes, coeffs)

    if recon:
        d[0,...], d[1,...] = full.real, full.imag
        fits.writeto(fits_file[:-5]+'_zern.fits', circular_mask(d), h, overwrite=True)
        d[0,...], d[1,...] = trunc.real, trunc.imag
        fits.writeto(fits_file[:-5]+'_zern_trunc.fits', circular_mask(d), h, overwrite=True)

    end = time.time()
    print "... Time taken: %.2f minutes"%((end-start)/60.)

def zernike_single_channel(h5, ant, freq, diameter, Npix, bandwidth, stokes='I', Nmodes=50, thresh=20):
    beam, dataset = load_data(h5, ant=ant, freqs=freq, bandwidth=bandwidth, Npix=Npix, diameter=2, save=False)
    beam = beam[0,...]
    
    beam = np.nan_to_num(beam)
    model = Zernike(beam, Nmodes, thresh, recon=True)
    beam = model.recon_trunc_J

    if stokes=='I': d = np.abs((beam[0,0,...]+beam[1,1,...])/2.)

    el = dataset.env_el[0]
    az = np.mean(dataset.scanaz*(180./np.pi))
    ID = dataset.filename.split('/')[-1][:-3]

    hdr = fits.Header()
    hdr['ID'] = ID
    ctypes = ['AZIMUTH', 'ELEVATION']
    crvals = [az, el]
    cdelts = [float(diameter)/d.shape[0], float(diameter)/d.shape[1]]
    cunits = ['deg', 'deg']
    for i in range(len(d.shape)):
        ii = str(i+1)
        hdr['CTYPE'+ii] = ctypes[i]
        hdr['CRVAL'+ii] = crvals[i]
        hdr['CDELT'+ii] = cdelts[i]
        hdr['CUNIT'+ii] = cunits[i]
    hdr['FREQ'] = freq[0]*1e6
    hdr['FWIDTH'] = bandwidth*1e6

    hdu = fits.PrimaryHDU(d, header=hdr)
    filename = h5[:-3]+'_beam_Zernike_model_%iMHz_%ideg.fits'%(freq[0],diameter)
    hdu.writeto(filename, clobber=True)
    print "--> Saved as %s"%filename
    
def rmse(predicted, targets):
    """
    Computes root mean squared error of two numpy ndarrays
    
    Args:
        predicted: an ndarray of predictions
        targets: an ndarray of target values
    
    Returns:
        The root mean squared error as a float
    """
    return (np.sqrt(np.mean((targets-predicted)**2)))



def gof_plot1(data, model, coeffs=None, vrange=[-20,0, -30,-15], extent = [-3,3,-3,3], view=True, title='',  cmap=plt.cm.jet, 
		fontsize = 20, ylabels=['Data', 'Model', 'Residual']):
    #
#    import matplotlib
#    matplotlib.rc('font', family='serif', serif='cm10')
#    matplotlib.rc('text', usetex=True)
#    matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
    #
#    model = model.real
#    data = data.real
    model = circular_mask(model)
    data = circular_mask(data)   
    
    
    fig = plt.figure(1, figsize=(15,10))
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=(3,4), axes_pad=0.0, add_all=True, share_all=False, aspect=True, \
                            label_mode='L', cbar_mode='none')
#    cmap = plt.cm.nipy_spectral
    c = 0
    ims = []
    residual = data - model #abs(data) - abs(model)
    for i in range(2):
        for j in range(2):
            if i==j: vmin, vmax = vrange[0], vrange[1]
            else: vmin, vmax = vrange[2], vrange[3]
#            if i != j:
#            	data[i,j,...]*=150
#            	model[i,j,...]*=150
            im = g[c].imshow(data[i,j,:,:], origin='lower', cmap=cmap, extent=extent)
#            im = g[c].imshow(20*np.log10(abs(data[i,j,:,:])), origin='lower', vmin=vmin, vmax=vmax, cmap=cmap, extent=extent)
            g[c].set_ylabel('%s' %ylabels[0], fontsize=fontsize)
            cax = inset_axes(g[c], loc=1, height='5%', width='85%', bbox_to_anchor=(-0.04,0,1,1.1), bbox_transform=g[c].transAxes)
            cb = plt.colorbar(im, cax=cax, orientation='horizontal')
            cb.ax.xaxis.set_ticks_position('top')
            cb.ax.xaxis.set_label_position('top')
            cb.locator = ticker.MaxNLocator(nbins=4)
            cb.update_ticks()
	    
	    im = g[c+4].imshow(model[i,j,...], origin='lower', cmap=cmap, extent=extent)
#            im = g[c+4].imshow(20*np.log10(abs(model[i,j,...])), origin='lower', vmin=vmin, vmax=vmax, cmap=cmap, extent=extent)
            g[c+4].set_ylabel('%s' %ylabels[1], fontsize=fontsize)
	    
#	    if i != j & (data[i,j,...].max > 0.8):
#            	data[i,j,...]/=150.
#            	model[i,j,...]/=150.
#    	    residual = data[i,j,...] - model[i,j,...]          
	    im = g[c+8].imshow(residual[i,j,...], origin='lower', extent=extent, cmap=cmap)  	    
	    
#    	    im = g[c+8].imshow(residual[i,j,...], origin='lower', vmin=vmin, vmax=vmax, extent=extent, cmap=cmap)
            g[c+8].set_ylabel('%s' %ylabels[2], fontsize=fontsize)
            cax = inset_axes(g[c+8], loc=4, height='5%', width='90%', bbox_to_anchor=(0,-.1,1,1.1), bbox_transform=g[c+8].transAxes)
#            cb = plt.colorbar(im, cax=cax, orientation='horizontal')
            cb = fig.colorbar(im, cax=cax, orientation='horizontal')
            cb.formatter.set_powerlimits((0, 0))
            cb.locator = ticker.MaxNLocator(nbins=5)
            cb.update_ticks()

            c+=1

    for i in range(12):
        g[i].set_xticklabels([])
        g[i].set_yticklabels([])
        
    fig.suptitle(title, fontsize=fontsize)

    if view==True: plt.show()
    plt.close()
    #
    # **************************************************************************************
#    x = np.linspace(extent[0],extent[1], data.shape[-1])
#    fig, ax = plt.subplots(1,4, figsize=(15,3))
#    c = 0
#    for i in range(2):
#        for j in range(2):
#            ax[c].plot(x,data[i,j,:,:].diagonal(), ':', label='Data')
#            ax[c].plot(x,model[i,j,:,:].diagonal(), '-', label='Model')
##            ax[c].plot(x,data[i,j,data.shape[2]/2,:], ':', label='Data')
##            ax[c].plot(x,model[i,j,model.shape[2]/2,:], '-', label='Model')
#            ax[c].legend(loc='upper right', framealpha=0.5)
#            c+=1
#    plt.subplots_adjust(wspace=0.2, hspace=0.05)
#    fig.suptitle('A line through the diagonal of the beams', fontsize=10)
#    if view==True: plt.show()
#    plt.close()
# **********************************************************************************************
    radial_profiles(data, model, extent[1], view, fontsize=fontsize)

    fig, ax = plt.subplots(1,4, figsize=(20,4))
    c=0
    for i in range(2):
        for j in range(2):
            res = residual[i,j,...]
            res = res[~np.isnan(res)]*1e2
            mu, sig = np.mean(res), np.std(res)
            n,b,p = ax[c].hist(res, 60, normed=True, facecolor='green', alpha=0.75, label="%.4f, %.4f"%(mu,sig))
#            print n
            #leg = ax[c].legend(loc='best')
            #leg.get_frame().set_facecolor('none')
            #leg.get_frame().set_linewidth('0.')
            ax[c].set_xlabel('Residual [%]', fontsize=fontsize)
            ax[c].set_title("$\\mu$=%.4f, $\\sigma$=%.4f"%(abs(mu),sig), fontsize=fontsize)
            ax[c].locator_params(nbins=5, axis='x')
            for tick in ax[c].xaxis.get_major_ticks(): tick.label.set_fontsize(18) 
	    for tick in ax[c].yaxis.get_major_ticks(): tick.label.set_fontsize(18)
            c+=1
    if view==True: plt.show()
    plt.close()
    
    if coeffs is not None:
	    fig, ax = plt.subplots(1,4, figsize=(20,4))
	    c=0
	    for i in range(2):
		for j in range(2):
		    ax[c].plot(range(coeffs.shape[-1]), coeffs[i,j,:], 'o:', markersize=5)
		    ax[c].set_xlim([-1,coeffs.shape[-1]+1])
		    ax[c].set_xlabel('Mode number',fontsize=fontsize)
		    ax[c].set_title('Coefficients',fontsize=fontsize)
		    ax[c].ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#		    ax[c].yaxis.get_major_formatter().set_powerlimits((0, 1))
		    for tick in ax[c].xaxis.get_major_ticks(): tick.label.set_fontsize(18) 
		    for tick in ax[c].yaxis.get_major_ticks(): tick.label.set_fontsize(18)
		    
		    c+=1
	    if view==True: plt.show()
	    plt.close()

def radial_profiles(data, model, extent=3, view=True,fontsize=15):
    residual = abs(data)-abs(model)
#    N = data.shape[-1]/2-2
    
    fig, ax = plt.subplots(1,4, figsize=(25,5))
    c = 0
    for i in range(2):
        for j in range(2):
            res = residual[i,j,...]
            res = np.nan_to_num(res)
            mu, sig = np.mean(res)*1e2, np.std(res)*1e1
#            mu, sig = np.mean(res), np.std(res)
            N = len(abs(azimuthalAverage(data[i,j,...])))
            x = np.linspace(0,extent, N)
#            if i !=j:  print abs(azimuthalAverage(model[i,j,...])); np.save('wow',abs(azimuthalAverage(model[i,j,...])))
            ax[c].loglog(x, abs(azimuthalAverage(data[i,j,...])), 'x', markersize=3, label='Measurement')
            ax[c].loglog(x, abs(azimuthalAverage(model[i,j,...])), '-', label='Model')
            ax[c].loglog(x, abs(azimuthalAverage(res)), '-', label='Mean error')
            
#            ax[c].loglog(x, abs(radial_profile_x(data[i,j,...])[:N]), 'x', markersize=3, label='Measurement')
#            ax[c].loglog(x, abs(radial_profile_x(model[i,j,...])[:N]), '-', label='Model')
#            ax[c].loglog(x, abs(radial_profile_x(res)[:N]), '-', label='Mean error')
            leg = ax[c].legend(loc='best',fontsize=fontsize)
            leg.get_frame().set_facecolor('none')
            leg.get_frame().set_linewidth('0.')
            ax[c].set_xlabel('Radius [deg]', fontsize=fontsize)
            ax[c].set_title("$\\mu$=%.2f, $\\sigma$=%.2f"%(abs(mu),sig), fontsize=fontsize)
            ax[c].set_xlim([1e-1,3.1])
            if i==j: ax[c].set_ylim([1e-5,1.2])
            if i!=j: ax[c].set_ylim([1e-4,5e-2])
            for tick in ax[c].xaxis.get_major_ticks(): tick.label.set_fontsize(18) 
	    for tick in ax[c].yaxis.get_major_ticks(): tick.label.set_fontsize(18)
            c+=1
    plt.subplots_adjust(wspace=0.2, hspace=0.05)
    if view==True: plt.show()
    plt.close()	
    
#def radial_profiles(data, model, extent=3, view=True):
#    residual = abs(data)-abs(model)
#    N = data.shape[-1]/2-2
#    x = np.linspace(0,extent, N)
#    fig, ax = plt.subplots(1,4, figsize=(25,5))
#    c = 0
#    for i in range(2):
#        for j in range(2):
#            res = residual[i,j,...]
#            res = np.nan_to_num(res)
##            mu, sig = np.mean(res)*1e2, np.std(res)*1e2
#            mu, sig = np.mean(res), np.std(res)
#            ax[c].loglog(x, abs(azimuthalAverage(data[i,j,...])), 'x', markersize=3, label='Measurement')
#            ax[c].loglog(x, abs(azimuthalAverage(model[i,j,...])), '-', label='Model')
#            ax[c].loglog(x, abs(azimuthalAverage(res)), '-', label='Mean error')
#            
##            ax[c].loglog(x, abs(radial_profile_x(data[i,j,...])[:N]), 'x', markersize=3, label='Measurement')
##            ax[c].loglog(x, abs(radial_profile_x(model[i,j,...])[:N]), '-', label='Model')
##            ax[c].loglog(x, abs(radial_profile_x(res)[:N]), '-', label='Mean error')
#            leg = ax[c].legend(loc='best')
#            leg.get_frame().set_facecolor('none')
#            leg.get_frame().set_linewidth('0.')
#            ax[c].set_xlabel('Radius [deg]')
#            ax[c].set_title("$\\mu$=%.2f, $\\sigma$=%.2f"%(mu,sig))
#            ax[c].set_xlim([1e-1,3.1])
#            if i==j: ax[c].set_ylim([1e-5,1.2])
#            if i!=j: ax[c].set_ylim([1e-4,5e-2])
#            c+=1
#    plt.subplots_adjust(wspace=0.2, hspace=0.05)
#    if view==True: plt.show()
#    plt.close()																
    				

def inspect_coeffs(coeffs, nu, idxs, show=True, Noll=False, ylim=[None,None], xlim=None, marker='.', original=[]):
    plt.rcParams['figure.figsize'] = (30,15)
    plt.rcParams['font.size'] = 12
    f, axes = plt.subplots(2,2)
    colrs = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'steelblue', 'gray', 'maroon']
    corrs = [['XX', 'XY'], ['YX', 'YY']]
    
    for i in range(2):
        for j in range(2):
            coeffs_avg = np.mean(coeffs[:,i,j,:], axis=0)
            sortidx = np.argsort(coeffs_avg)[::-1]
            ax = axes[i,j]
            count = 0
            for idx in idxs:
                colr = colrs[count]
                if Noll:
                    ind = idx-1; label = idx
                else: ind = sortidx[idx]; label = sortidx[idx]+1
                ax.plot(nu, coeffs[:,i,j,ind], marker, color=colr, label="%i"%label)
                try: ax.plot(nu, original[:,i,j,ind], 'o', markersize=2, color=colr, label="%i data"%label)
                except: None
                if i==j: ax.set_ylim(ylim[0])
                elif i!=j: ax.set_ylim(ylim[1])
                count += 1
            leg = ax.legend(ncol=5, loc='best')
            leg.get_frame().set_facecolor('none')
            leg.get_frame().set_linewidth('0.')
            ax.set_xlabel('Frequency [MHz]')
            ax.set_ylabel('Energy of a coefficient by Noll index')
            ax.set_title(corrs[i][j])
            if xlim==None: ax.set_xlim([min(nu), max(nu)+1])
            else: ax.set_xlim(xlim)
    #plt.suptitle('Real (solid) and Imag (dashed) highest-energy coeffs: %i to %i'%(idxs[0], idxs[-1]), fontsize=14)
    
    if show: plt.show()

class Zernike_opticspy(object):
    def __init__(self, data, Nm_co=12, Nm_cross=12):
        self.data = data
        self.model = np.zeros(data.shape)
        self.residual = np.zeros(data.shape)
        self.coeffs = np.zeros((2,2,37))
        l = data.shape[-1]
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if i==j: Nmodes = Nm_co
                elif i!=j: Nmodes = Nm_cross
                fl, C = opticspy.zernike.fitting(self.data[i,j,...], Nmodes, removepiston=False)
                self.model[i,j,:,:] = circular_mask(C.zernikematrix(l=l))
                self.coeffs[i,j,:] = fl
                self.residual[i,j,...] = abs(self.data[i,j,...]) - abs(self.model[i,j,...])
                maxres = np.nanmax(abs(self.residual[i,j,...]))
                print "Fitted %ith order Zernike to data %i%i, maximum residual %f"%(Nmodes,i,j,maxres)


from scipy.interpolate import splrep, splev

def coeffs_interpolate(coeffs, freqs):
    coeffs_interp = np.zeros(coeffs.shape, dtype=np.complex)
    for i in range(2):
        for j in range(2):
            for c in range(coeffs.shape[-1]):
                spl = splrep(freqs, coeffs[:,i,j,c].real)
                coeffs_interp[:,i,j,c].real = splev(freqs, spl)
                spl = splrep(freqs, coeffs[:,i,j,c].imag)
                coeffs_interp[:,i,j,c].imag = splev(freqs, spl)
    return coeffs_interp


if __name__=='__main__':
    if sys.argv[1]=='allfreq': zernike_all_freqs(sys.argv[2], Nmodes=int(sys.argv[3]), threshold=int(sys.argv[4]))

