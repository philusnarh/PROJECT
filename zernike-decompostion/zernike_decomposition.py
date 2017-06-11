#
"""This script was written to deal with the decomposition of a 2-D image by set of Zernike functions"""

# ++++++++ Import python modules ++++++++++++++
from scipy.misc import factorial as fac
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys
#

# ++++++++++++++  Zernike decomposition functions  +++++++++++++++++++++++++

def zernike_num_coeff(n):
	"""
	computes the number of zernike coefficients
	given the radial polynomial order/level n.
	-----------------------------------------

	parameter::

	n: scalar (i.e. must be >=0). 
	   Number of radial order 

	return:: 
		output int value, representing the number of zernike coeffcients
	"""
    
	if not (n>=0):
		print('Input parameter must be >= 0')
		raise AssertionError()    
    
	return sum(xrange(n+1)) + n+1
	
	

def zernike_rad(m, n, rho): 
	"""
	 Function to compute Zernike Radial Polynomials:
	 -----------------------------------------------
	 
	 parameters::	 
	 
	   rho: float, scalar or array-like (i.e. 0 <= rho <= 1).
	   	Radius describing a unit disk
	   	
	   m: scalar (i.e. must be <= n).
	   	the repetition of Zernike moment
	   	
	   n: scalar (i.e. must be >= 0).
	   	the order of Zernike polynomial	   
	"""
   
	rad = np.zeros(rho.shape, rho.dtype)
	P = (n - abs(m)) / 2
	Q = (n + abs(m)) / 2

	for s in xrange(P + 1):
		c = (-1) ** s * fac(n - s)
		c /= fac(s) * fac(Q - s) * fac(P - s)
		rad += c * rho ** (n - 2 * s)
        
	return rad
	
	
	
def zernike_poly(m, n, rho, phi):
	"""
	Function to compute the Zernike polynomial (n, m) given a grid of radial
	coordinates rho and azimuthal coordinates phi.
	----------------------------------------------------------

	>>> zernike(3,5, 0.12345, 1.0)
	0.0073082282475042991
		
	parameters::
	
	m: scalar (i.e. must be <= n).
	   	the repetition of Zernike moment
	   	
   	n: scalar (i.e. must be >= 0).
	   	the order of Zernike polynomial
	   	
   	rho: float, scalar or array-like (i.e. 0 <= rho <= 1).
	   	Radius describing a unit disk
	   	
   	phi: float or array-like (i.e. 0 <= phi <= 2*pi)
   		the azimuthal angle.
	
	"""
	
	if (m > 0): return zernike_rad(m, n, rho) * np.cos(m * phi)
	if (m < 0): return zernike_rad(-m, n, rho) * np.sin(-m * phi)

	return zernike_rad(0, n, rho)
	


def  zernike_Double_Index(nlevels):
    
    """    
	Computes 2 indices (n,m) to specify Zernike functions, starting
	from the top, shifts to the left and then right.
	--------------------------------------------------------------
	parameter::

	nlevels: scalar (i.e. must be >=0). Number of radial order  

	return:

	n: array_like. the radial index
	 
	m: array like. the meridional index

	NB: m <= n

    """
    
	    
    if not (nlevels>=0):
        print('Input parameter nlevels must be >= 0')
        raise AssertionError()
    
    if (nlevels == 0):
        
        m = 0
        n = 0
        
        return n, m
    
    else:
        
        # ++++ Defining layout for row number n and colunmn number m ++++++++

        row_n = nlevels+1
        col_m = 2*nlevels +1
        x = np.arange(row_n)
        y = np.arange(-(col_m-1)//2, (col_m+1)//2,1)
        Q = [(i,j) for i in x for j in y]
        #


        nm_index = []
    
        top = (col_m + 1)/2
        leftside = row_n*col_m - col_m + 1
        rightside = row_n*col_m    

        k1 = 0; k2 = 0

        for i in xrange(top,row_n*col_m+1, 2*col_m):

            nm_index.append(Q[i-1])
            s1 = i + col_m + 1
            s2 = i + col_m - 1 
            jj1 = k1
            jj2 = k2


            while (s2 <= leftside): 

              nm_index.append(Q[s2-1])
              s2 +=col_m - 1
              jj1 += 1
              jj2 -= 1

            leftside +=2

            jj1 = k1
            jj2 = k2

            while (s1 <= rightside):    

        #       
              nm_index.append(Q[s1-1])
              s1 +=col_m + 1
              jj1 += 1
              jj2 += 1

            rightside -=2
            k1 = 0; k2 += 2

        n = np.array(nm_index)[:,0]
        m = np.array(nm_index)[:,1]

        return n, m

        
def cart2pol(x, y):
	"""
	Transform Cartesian coordinates to polar 
	
	parameters::
	
	x: array_like.
		corresponding entries in matrices x 
	
	y: array_like.
		corresponding entries in matrices y
	
	return::
	
	rho: array_like.  (i.e. 0 <= rho <= 1).
	   	Radius describing a unit disk
	   	
   	phi: array-like (i.e. 0 <= phi <= 2*pi)
   		the azimuthal angle.
	""" 

	rho = np.sqrt(x**2 + y**2)
	phi = np.arctan2(y, x)

	return rho, phi
	

def unit_disk(imgSize):
        """Create an unit disk and convert to rho, phi and 
            and masking the grid.
           
           parameter::
           
           imgSize: scalar, image Nside"""
           
#        src = np.nan_to_num(imgsrc)
	if not (imgSize > 0):
		print('Nside must be > 0')
		raise AssertionError()
		
        nx = imgSize  #, ny = src.shape
        grid = (np.indices((nx, nx), dtype=np.float) - nx/2) / (nx*1./2) # create unit grid [-1,1]        
        grid_rho, grid_phi = cart2pol(x=grid[0], y=grid[-1])
#        grid_rho = (grid**2.0).sum(0)**0.5 # rho = sqrt(x^2+y^2)
#        grid_phi = np.arctan2(grid[0], grid[-1]) # phi = itan(x/y)
        grid_mask = grid_rho <= 1 
        
        return  grid_rho, grid_phi, grid_mask
        
               
def zernike_basis_n_coeffs(imgsrc, nlevels, return_basis=True):
	"""Description: Represent a 2D image as a sum of Zernike polynomials using
	      a matrix inversion scheme.
	      
	% 
	% This function attempts to solve the C_{n}^{m}'s in equation,
	% 
	%         	                          
	%                    
	%		      M__
	%                     \
	% 	 W(rho,phi) = /__  C_{n}^{m} * Z_{n}^{m}(rho,phi)
	%                    m,n
	% 
	% where the Z_{n}^{m}(rho,phi)'s are the Zernike polynomials from the 'zernike_poly'
	% function, W(rho,phi) is the 2D image (i.e. wavefront) to be represented as a sum of 
	% Zernike polynomials, the C_{n}^{m}'s are the Zernike coefficients, and M (i.e. M=nlevels)
	% is the number of Zernike polynomials to use.
	%
	% Input:    
	% 
	%	imgsrc:- n X n (square) array. 2D image to be represented as a sum of Zernike polynomials (i.e. W(rho,phi))
	%                 
	%       nlevels: Number of Zernike polynomials to use (i.e. M)
	% Output:   
	%  	C - Zernike coefficients (C_{n}^{m}'s) as a vector
	%      NB: also returns an (optional) basis of zernike polynomials
	%	
	"""	

	im = np.nan_to_num(imgsrc)	

	
	# Computes 2 indices (n,m) to specify Zernike functions
	n,m = zernike_Double_Index(nlevels)
	
	# Generate (rho, phi) grids and masking grid	
	grid_rho, grid_phi, grid_mask = unit_disk(im.shape[0])
	
	# Generate the Zernike basis to the nlevels 
	basis=[zernike_poly(i, j, grid_rho, grid_phi)*grid_mask for i,j in zip(m,n)]	

	# Compute covariance between all Zernike polynomials
	cov_mat = np.array([[np.sum(zerni * zernj) for zerni in basis] for zernj in basis])

	# Invert covariance matrix using SVD (  A x = b  ==>  x =>  x= A^{pseudoinv} b)
	cov_mat_in = np.linalg.pinv(cov_mat)	
	
	# +++ Compute inner product between the img and the Zernike bases
	innerprod = np.array([np.sum(im * zerni) for zerni in basis])
	
	# Computing beam coefficients in the Zernike basis to the Nth order 
	# Given the inner product vector of the input image with Zernike basis,
	# calculate the Zernike polynomial coefficients
	# (A^T A)^(-1) A^T b = x
	# innerprod === A^T b
	# coefficients == (A^T A)^(-1)
	coefficients = np.dot(cov_mat_in, innerprod)
	
	if return_basis is True: return coefficients, np.array(basis)
	    
	else: return coefficients
	

def truncate(coeffs, threshold=99):
        """Truncate the coefficients upto the given threshold
           parameter::
             threshold: scalar"""
        sortedindex = np.argsort(np.abs(coeffs))[::-1]
        Ncoeff = coeffs.shape[-1]
        cutoff = np.int(np.round(Ncoeff*threshold/100.))
        
#        print "Keeping %2.0f %% (N=%s) of the biggest coefficients"%(threshold,cutoff)

        coeffs_trunc = coeffs.copy() 			# copy of all coeff
        coeffs_trunc[sortedindex[cutoff:]] = 0 		# put coeff
        
        return coeffs_trunc
        	
	    
def reconstruct_image(imgsrc, nlevels, trunc_threshold=None):
	"""Reconstruct a model image from the coeffcicients
	"""
	
	coeffs = zernike_basis_n_coeffs(imgsrc, nlevels, return_basis=False)
	
	if trunc_threshold is not None:
		coeffs = truncate(coeffs, threshold=trunc_threshold)
		
	# Generate (rho, phi) grids and masking grid	
	grid_rho, grid_phi, grid_mask = unit_disk(imgsrc.shape[0])
	
	# Computes 2 indices (n,m) to specify Zernike functions
	n,m = zernike_Double_Index(nlevels)
	
	reconstr_im = np.sum(val * zernike_poly(m[i], n[i], grid_rho, grid_phi)*grid_mask for (i, val) in enumerate(coeffs))
	
	return reconstr_im


def  zernike_visuo__pyramid(zbasis, n, m, nlevels, figsize=(12, 12), cmap='jet', fontsize=20, colorbar_labelsize=10):
    
    """
    figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    2Dpolynomials
    Computes 2 indices (n,m) to specify Zernike functions, starting
    from the top, shifts to the left and then right.
    --------------------------------------------------------------
    parameter::
    
    nlevels: int value (i.e. must be >=0). Number of radial order  
    
    return:
    
    n: array_like. the radial index 
    m: array like. the meridional index
    
    NB: m <= n
    
    """
    
    cmap = plt.get_cmap('%s' %cmap)
    
    index = 0
    if not (nlevels>=0):
        print('Input parameter must be >= 0')
        raise AssertionError()        
    
    axlist = []
    if (nlevels == 0):
        
        fig = plt.figure(num = 1, figsize=figsize)
        ax = fig.add_subplot(1,1,1)
        axlist.append(ax)
        im = ax.imshow(zbasis, cmap=cmap, interpolation='lanczos')
        ax.set_title(r'$Z_{%d}^{%d}$' %(n,m), fontsize=fontsize)
	ax.axis('off')

    
    else:
        
        # ++++ Defining layout for row number n and colunmn number m ++++++++
        
        fig = plt.figure(1, figsize=figsize)
        row_n = nlevels + 1
        col_m = 2*nlevels + 1

        top = (col_m + 1)/2
        leftside = row_n*col_m - col_m + 1
        rightside = row_n*col_m    

        k1 = 0; k2 = 0
        

        for i in xrange(top,row_n*col_m+1, 2*col_m):

            ax = fig.add_subplot(row_n,col_m,i)
            axlist.append(ax)
            im=ax.imshow(zbasis[index], cmap=cmap, interpolation='lanczos', alpha=None)
            ax.set_title(r'$Z_{%d}^{%d}$' %(n[index],m[index]), fontsize=fontsize)
            ax.axis('off')
            index += 1
            s1 = i + col_m + 1
            s2 = i + col_m - 1 
            jj1 = k1
            jj2 = k2


            while (s2 <= leftside): 

              ax = fig.add_subplot(row_n,col_m,s2)
              axlist.append(ax)
              im=ax.imshow(zbasis[index], cmap=cmap, interpolation='lanczos')
              ax.set_title(r'$Z_{%d}^{%d}$' %(n[index],m[index]), fontsize=fontsize)
              ax.axis('off')
              index += 1
              s2 +=col_m - 1
              jj1 += 1
              jj2 -= 1

            leftside +=2

            jj1 = k1
            jj2 = k2

            while (s1 <= rightside):
                
              ax = fig.add_subplot(row_n,col_m,s1)
              axlist.append(ax)
              im=ax.imshow(zbasis[index], cmap=cmap, interpolation='lanczos')
              ax.set_title(r'$Z_{%d}^{%d}$' %(n[index],m[index]), fontsize=fontsize)
              ax.axis('off')
              index += 1
              s1 +=col_m + 1
              jj1 += 1
              jj2 += 1

            rightside -=2
            k1 = 0; k2 += 2


    cbar = fig.colorbar(im, ax=axlist,fraction=0.04, orientation='horizontal')        
    cbar.ax.tick_params(labelsize=colorbar_labelsize)
    fig.subplots_adjust(wspace=0,hspace=0, right=0.78, bottom=0.2)
    fig.savefig('zernike_orders.png', dpi=300)

    return None


