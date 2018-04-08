import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate




def plot2d_multi(d, sh=True, vrange=None, diameter=6, cbticks=None, 
		cblabel='', title=['I', 'Q', 'U', 'V'], figsize=(8,8), svfigname=None):
    
#    cm = plt.get_cmap('%s' %cmp) #'jet')#

    fig = plt.figure(1, figsize=figsize)
    nrc = (d.shape[0], d.shape[1])
    ngrids = d.shape[0] * d.shape[1]
    extent = [-diameter/2, diameter/2, -diameter/2, diameter/2]
    g = axes_grid.ImageGrid(fig, 111, nrows_ncols=nrc, axes_pad=0.0, add_all=True,\
                            share_all=False, aspect=True, label_mode='1', cbar_mode='none')

    c = 0
    ims = []
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            # Fix vmin, vmax and color bar ticks
            if vrange!=None:
                if i!=j: vmin, vmax = vrange[2], vrange[3]
                elif i==j: vmin, vmax = vrange[0], vrange[1]
            elif vrange==None: vmin, vmax, cticks = None, None, None
# plt.cm.cubehelix
            im = g[c].imshow(d[i,j,:,:], extent=extent, origin='lower',  vmin=vmin, vmax=vmax)
            if title is not '':
		    if i == 0:
		    	g[j].set_title(r'%s' %title[j], fontsize=15)
		    
#             print i
#            cmap=plt.cm.nipy_spectral  cmap=plt.cm.jet,
            ims.append(im)

            # Color bar
            if ngrids==4: rights=[1,3]
            elif ngrids==16: rights=[3,7,11,15]
            elif ngrids==12: rights=[3,7,11]
            if c in rights:
                cax = inset_axes(g[c], loc=5, width='5%', height='96%', bbox_to_anchor=(0.1,0,1,1), \
                    bbox_transform=g[c].transAxes)
                cb = plt.colorbar(im, cax=cax, format = '%.0e', orientation='vertical') # shrink=0.9, pad = 0.05)
                from matplotlib import ticker
                tick_locator = ticker.MaxNLocator(nbins=3)
                cb.locator = tick_locator
                cb.update_ticks()
                if vrange != None: cb.set_ticks(range(vmin, vmax+1,5))
                cb.ax.xaxis.set_ticks_position('top')
#                 if cblabel is not None:
                cb.ax.set_ylabel('%s' %cblabel)   #%Power [decibel] Power [dB]

            # Ticks and circles
            g[c].set_ylabel('Angular distance [deg]', fontsize = 15)
            g[c].set_xlabel('Angular distance [deg]', fontsize = 15)
#            if i == 0:
#            	g[i].set_title(r'%s' %title[j])
            cn, o = (extent[0]+extent[1])/2, abs(extent[0])/3.
            g[c].add_artist(plt.Circle((cn,cn), o*1, color='black', linestyle='dashed', fill=False))
            g[c].add_artist(plt.Circle((cn,cn), o*2, color='black', linestyle='dashed', fill=False))
            g[c].add_artist(plt.Circle((cn,cn), o*3, color='black', linestyle='dashed', fill=False))
            c += 1
#    g[0].text(0,d.shape[3]+15, title, fontsize=13)

#    fig.subplots_adjust(wspace=0,hspace=0,left=0.05,right=.87,bottom=-.13,top=1.1)
#    matplotlib.tight_layout.auto_adjust_subplotpars(fig, renderer, nrows_ncols, num1num2_list, subplot_list, ax_bbox_list=None, pad=1.08, h_pad=None, w_pad=None, rect=None)
#     plt.tight_layout()
#    fig.tight_layout(pad=0, h_pad=None, w_pad=None, rect=None) #  h_pad=None, w_pad=None
    plt.subplots_adjust(hspace=0.0, wspace=0.0)
    if svfigname is not None: fig.savefig('%s' %svfigname, dpi=100)
    
    if sh==True: plt.show(); plt.close()



def jones2mueller(xx_real, xx_imag, xy_real, xy_imag, yx_real, yx_imag, yy_real, yy_imag):
        #
        '''
        Extracts the Jones components into XX, XY, YX, YY
        to produce Mueller components'''
        
        #
        #   Generating Jones Terms
        #
        xx = xx_real + 1j*xx_imag
        xy = xy_real + 1j*xy_imag
        yx = yx_real + 1j*yx_imag
        yy = yy_real + 1j*yy_imag
        
        M = []
        #    
        #  Generating Mueller Terms
        #   
        m_ii = 0.5*(xx*np.conjugate(xx) + xy*np.conjugate(xy) + yx*np.conjugate(yx) + yy*np.conjugate(yy))
        m_iq = 0.5*(xx*np.conjugate(yx) - xy*np.conjugate(xx) + yx*np.conjugate(yy) + yy*np.conjugate(xy)) 
        m_iu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) + yx*np.conjugate(yy)+ yy*np.conjugate(yx)) 
        m_iv = 0.5*1j*(xx*np.conjugate(xy) + yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx))
        
        M.append([m_ii, m_iq, m_iu, m_iv])
        
        m_qi = 0.5*(xx*np.conjugate(yx) - xy*np.conjugate(xx) + yx*np.conjugate(yy) + yy*np.conjugate(xy))   
        m_qq = 0.5*(xx*np.conjugate(xx) - xy*np.conjugate(xy) - yx*np.conjugate(yx)+ yy*np.conjugate(yy)) 
        m_qu = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(yy) - yy*np.conjugate(yx)) 
        m_qv = 0.5*1j*(xx*np.conjugate(xy) - yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx))
        
        M.append([m_qi, m_qq, m_qu, m_qv])
        
        m_ui = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) + yx*np.conjugate(yy)+ yy*np.conjugate(yx))
        m_uq = 0.5*(xx*np.conjugate(xy) + xy*np.conjugate(xx) - yx*np.conjugate(yy) - yy*np.conjugate(yx))  
        m_uu = 0.5*(xx*np.conjugate(yy) + yy*np.conjugate(xx) + xy*np.conjugate(yx) + yx*np.conjugate(xy))         
        m_uv = 0.5*1j*(xx*np.conjugate(yy) - yy*np.conjugate(xx) - xy*np.conjugate(yx) - xy*np.conjugate(xy))
        
        M.append([m_ui, m_uq, m_uu, m_uv])
        
        m_vi = 0.5*1j*(xx*np.conjugate(xy) + yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx))  
        m_vq = 0.5*1j*(xx*np.conjugate(xy) - yx*np.conjugate(yy) - xy*np.conjugate(xx) - yy*np.conjugate(yx))       
        m_vu = 0.5*1j*(-xx*np.conjugate(yy) + yy*np.conjugate(xx) - xy*np.conjugate(yx) + xy*np.conjugate(xy))    
        m_vv = 0.5*(xx*np.conjugate(yy) - yx*np.conjugate(xy) + yy*np.conjugate(xx) - xy*np.conjugate(yx))
        M.append([m_vi, m_vq, m_vu, m_vv])
        #
        
        return   np.array(M).real 


#def jones_to_mueller(f1,f2=None, sv=True):
#    """
#    Convert 2x2xNxN Jones matrix into 4x4xNxN Mueller matrix.
#    If f1=f2, compute the autocorrelation version for single dishes
#    """
#    if f2==None: f2 = f1
#    # a, b = np.array(f1) #np.load(f1), np.load(f2)
#    a = np.array(f1)
#    b = np.array(f1)
#    M = np.zeros((a.shape[0],4,4,a.shape[3],a.shape[4]), dtype='c16')
#    S = 0.5*np.matrix('1 1 0 0; 0 0 1 1j; 0 0 1 -1j; 1 -1 0 0')
#    for f in range(a.shape[0]):
##        print f
#        for i in range(a.shape[3]):
##            if i in [50,100,200,300,400,500]: print i
#            for j in range(a.shape[4]):
#                ab = np.kron(a[f,:,:,i,j],b[f,:,:,i,j].conj())
#                M[f,:,:,i,j] = np.dot( np.dot(np.linalg.inv(S), ab), S )
#    if sv==True: np.save(f1[:-4]+'_M.npy', M)
#    return M

def jones_to_mueller(f1,f2=None, sv=False):
    """
    Convert 2x2xNxN Jones matrix into 4x4xNxN Mueller matrix.
    If f1=f2, compute the autocorrelation version for single dishes
    """
    if f2==None: f2 = f1
    # a, b = np.array(f1) #np.load(f1), np.load(f2)
    a = np.array(f1)
    b = np.array(f1)
    M = np.zeros((a.shape[0],4,4,a.shape[3],a.shape[4]), dtype=complex) #dtype='c16'
    
    S = 0.5*np.matrix('1 1 0 0; 0 0 1 1j; 0 0 1 -1j; 1 -1 0 0')
    for f in range(a.shape[0]):
#        print f
        for i in range(a.shape[3]):
#            if i in [50,100,200,300,400,500]: print i
            for j in range(a.shape[4]):
                ab = np.kron(a[f,:,:,i,j],b[f,:,:,i,j].conj())
                M[f,:,:,i,j] = np.dot( np.dot(np.linalg.inv(S), ab), S )
    if sv==True: np.save(f1[:-4]+'_M.npy', M)
    return M

def zeromask(channel, eps = 3.8800000000000004e-9):
    
    return np.ma.masked_where(abs(channel) < eps,channel , copy=True)

# ++++++++++ interpolating to have same size ++++++++

#
def rescale_data_size(data, newsizex, newsizey):
    dshape = data.shape
    # define new size
    outKSize_x = newsizex
    outKSize_y = newsizey
    
    # Rescale Data Size
    x_old = np.linspace(-dshape[0]/2., dshape[0]/2., dshape[0])      
    y_old = np.linspace(-dshape[-1]/2., dshape[-1]/2., dshape[-1])
    xnew = np.linspace(x_old.min(), x_old.max(), outKSize_x)
    ynew =  np.linspace(y_old.min(), y_old.max(), outKSize_y)
    
    # Perform Interpolation
    interp_Fxn = interpolate.RectBivariateSpline(np.sort(x_old),np.sort(y_old),data, kx=3,ky=3) 

    return interp_Fxn(xnew,ynew)
