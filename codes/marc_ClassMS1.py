import numpy as np
from pyrap.tables import table
import os
import pylab
from pyrap.images import image
import sys
class ClassMS1():

    """
    

    Measurement set class
    """

    def __init__(self,MSname,UVW_colname="UVW",work_colname="DATA",dltime=100,overlap=0):
      
        table_all=table(MSname,ack=False,readonly=False)
        A0=table_all.getcol('ANTENNA1')
        A1=table_all.getcol('ANTENNA2')
        vis_all=table_all.getcol(work_colname)
        na=np.max(A1)+1
        nbl=(na*(na-1))/2+na
        NChan=vis_all.shape[1]
        NPol=vis_all.shape[2]
        Nbin=vis_all.shape[0]
        UVW=table_all.getcol(UVW_colname)
        Ntimes=(vis_all[(A0==0)&(A1==1)]).shape[0]
        table_all.close()
        t2=table(MSname+"SPECTRAL_WINDOW",readonly=False)
	wave=3e8/t2.getcol("REF_FREQUENCY")[0]
        t2.close()
        
        self.MSName=MSname
        self.ColNameIn=work_colname
        self.na=na
        self.NChan=NChan
        self.NPol=NPol
        self.Nbin=Nbin
        self.nbl=nbl
        self.data=vis_all
        self.A0=A0
        self.A1=A1
        self.UVW=UVW
        self.Ntimes=Ntimes
        self.dltime=dltime
    	self.wave=wave
	self.overlap=overlap
    def NormalAvgALBL(self):
         #try:
             _inddata=(self.A0==0)&(self.A1==1)
             _data=(self.data[_inddata].copy())
             _data = _data[_data.shape[0]/2 -self.dltime/2 :_data.shape[0]/2+self.dltime/2,:,:]
        
        
             data_result=np.zeros_like(self.data[0:self.nbl,:,:]).copy()
             _A0=self.A0[0:self.nbl].copy()
             _A1=self.A1[0:self.nbl].copy()
             for i in range(self.na):
                 for j in range(i,self.na):
                     
                     _ind_data=(_A0==i)&(_A1==j)
                     
                     _inddata=(self.A0==i)&(self.A1==j)
                     _data=(self.data[_inddata].copy())
                     _data = (_data[_data.shape[0]/2 -self.dltime/2 :_data.shape[0]/2+self.dltime/2,:,:]).copy()
                
                     data_result[_ind_data,:,:]=_data.mean(axis=0)

             return data_result, data_result[(_A0==0)&(_A1==1)].shape[0]
         #except:
          #   print "absolute value intergration time equal to zeros or gratter thant ",self.data.shape[0]
          #   sys.exit()

    def ConvolveALBL(self,function,var=None):
        _inddata=(self.A0==0)&(self.A1==1)
        _data=(self.data[_inddata].copy())
        _data = _data[_data.shape[0]/2 -(self.dltime/2) :_data.shape[0]/2+(self.dltime/2),:,:]
        
        
        data_result=np.zeros_like(self.data[0:self.nbl,:,:]).copy()
        _A0=self.A0[0:self.nbl].copy()
        _A1=self.A1[0:self.nbl].copy()

        if var!=None:
            for i in range(self.na):
                for j in range(i,self.na):
                    ind=(_A0==i)&(_A1==j)
                    ind1=(self.A0==i)&(self.A1==j)
                    dat=self.data[ind1].copy()
                    uvw=self.UVW[ind1].copy()
                    dat= (dat[dat.shape[0]/2 -(self.dltime/2 +self.overlap/2) : dat.shape[0]/2+(self.dltime/2 +self.overlap/2),:,:]).copy()
                    uvw = (uvw[uvw.shape[0]/2 -(self.dltime/2 +self.overlap/2) : uvw.shape[0]/2+(self.dltime/2 +self.overlap/2),:]).copy()
              
                    data_result[ind,:,:]=self.__convolBlUvwNtimMult__(dat,function,uvw, 2*(self.dltime/2 +self.overlap/2)   )
        else:
            weight=(self.WeightForALBL(function, 2*(self.dltime/2)  )).copy()
            for i in range(self.na):
                for j in range(i+1,self.na):
                    ind=(_A0==i)&(_A1==j)
                    ind1=(self.A0==i)&(self.A1==j)
                    dat=self.data[ind1].copy()
                    dat = (dat[dat.shape[0]/2 -(self.dltime/2) :dat.shape[0]/2+(self.dltime/2),:,:]).copy()
                    data_result[ind,:,:]=self.__convolBlSameNtimMult__(dat,function,weight,  2*(self.dltime/2))
        return data_result

    def __fragmentData__(self,__data):
         __data_c=(__data[0:(__data.shape[0]/self.dltime)*self.dltime,:,:]).copy()
         __data_r=(__data[__data_c.shape[0]:__data.shape[0],:,:]).copy()
         return  __data_c,__data_r

    def __defragmentData__(self,__data_c,__data_r):
        __data_c_r=np.zeros((__data_c.shape[0]+1,self.NChan,self.NPol),dtype=complex)
        __data_c_r[0:__data_c.shape[0],:,:]=__data_c.copy()
        __data_c_r[__data_c.shape[0]:self.Ntimes,:,:]=__data_r.copy()
        return __data_c_r

    def __fragmentUvw__(self,__uvw):
         __uvw_c=(__uvw[0:(__uvw.shape[0]/self.dltime)*self.dltime,:]).copy()
         __uvw_r=(__uvw[__uvw_c.shape[0]:__uvw.shape[0],:]).copy()
         return  __uvw_c,__uvw_r
   
         
    def __convolBlUvwNtimMult__(self,data,function,uvw,dltime):
        weightList=self.WeightIvProBL(function,uvw,dltime)
        dat=data.reshape(dltime,4)
        dat=dat.T
        dat=dat.reshape(4,1,dltime)
        __data=((dat*weightList).sum(axis=2)/weightList.sum(axis=2))
        return __data.T

    def __convolBlUvwNtimNmult__(self,data,function,uvw):
        __uvw_c,__uvw_r=self.__fragmentUvw__(uvw)
        __uvw_c=__uvw_c.reshape(__uvw_c.shape[0]/self.dltime,self.dltime,__uvw_c.shape[1]).copy()
        weightList=np.array([self.WeightIvProBL(function,__uvw_c[i,:,:],self.dltime) for i in range(__uvw_c.shape[0])]).copy()
        __data_c,__data_r=self.__fragmentData__(data)
        __data_c=(__data_c.reshape(__data_c.shape[0]/self.dltime,self.dltime,self.NChan,self.NPol)).copy()
        __data=np.array([__data_c[i,:,:,:].T for i in range(__data_c.shape[0])]).copy()
        __data_c=((__data*weightList).sum(axis=3)/weightList.sum(axis=3)).reshape(__data.shape[0],__data.shape[2],__data.shape[1]).copy()
        __uvw_r=__uvw_r.reshape(__uvw_r.shape[0]/__uvw_r.shape[0],__uvw_r.shape[0],__uvw_r.shape[1])
        __data_r=(__data_r.reshape(__data_r.shape[0]/__data_r.shape[0],__data_r.shape[0],self.NChan,self.NPol)).copy()      
        weightList=np.array([self.WeightIvProBL(function,__uvw_r[i,:,:],__uvw_r.shape[1]) for i in range(__data_r.shape[0])]).copy() 
        __data=np.array([__data_r[i,:,:,:].T for i in range(__data_r.shape[0])]).copy()
        __data_r=((__data*weightList).sum(axis=3)/weightList.sum(axis=3)).reshape(__data.shape[0],__data.shape[2],__data.shape[1]).copy()
        __data=self.__defragmentData__(__data_c,__data_r)
        return __data

    def __convolBlSameNtimMult__(self,data,function,weight,dltime):
        dat=data.reshape(dltime,4)
        dat=dat.T
        dat=dat.reshape(4,1,dltime)
        __data=((dat*weight).sum(axis=2))/weight.sum(axis=2)
        return __data.T

    def __convolBlSameNtimNmult__(self,data,function,weight_c,weight_r):
        __data_c,__data_r=self.__fragmentData__(data)
        #weight=self.WeightForALBL(function,self.dltime)
        __data_c=(__data_c.reshape(__data_c.shape[0]/self.dltime,self.dltime,self.NChan,self.NPol)).copy()
        weightList=np.array([weight_c for i in range(__data_c.shape[0])]).copy()
        __data=np.array([__data_c[i,:,:,:].T for i in range(__data_c.shape[0])]).copy()
        __data_c=((__data*weightList).sum(axis=3)/weightList.sum(axis=3)).reshape(__data.shape[0],__data.shape[2],__data.shape[1])
        #weight=self.WeightForALBL(function,__data_r.shape[0])
        __data_r=(__data_r.reshape(__data_r.shape[0]/__data_r.shape[0],__data_r.shape[0],self.NChan,self.NPol)).copy()
        weightList=np.array([weight_r for i in range(__data_r.shape[0])]).copy()
        __data=np.array([__data_r[i,:,:,:].T for i in range(__data_r.shape[0])]).copy()
        __data_r=((__data*weightList).sum(axis=3)/weightList.sum(axis=3)).reshape(__data.shape[0],__data.shape[2],__data.shape[1])
        __data=self.__defragmentData__(__data_c,__data_r)
        return __data

    def WeightIvProBL(self,function,uvw,dl):
        if dl%2==0:
            u0=(uvw[:,0][dl/2-1]+uvw[:,0][dl/2])/2.
            v0=(uvw[:,1][dl/2-1]+uvw[:,1][dl/2])/2.
            w0=(uvw[:,2][dl/2-1]+uvw[:,2][dl/2])/2.
        else:
            u0=(uvw[:,0][dl/2])
            v0=(uvw[:,1][dl/2])
            w0=(uvw[:,2][dl/2])
        vect=np.zeros((self.NPol,self.NChan,dl),dtype=complex)
        d=np.sqrt((uvw[:,0]-u0)**2+(uvw[:,1]-v0)**2+(uvw[:,2]-w0)**2)
        vect[:,:,:]=function(d/self.wave)
        return vect
    
    def WeightForALBL(self,function,dl):
        vect=np.zeros((self.NPol,self.NChan,dl),dtype=complex)
        ref_pixel=(dl-1)/2.
        for i in range(dl):
            vect[:,:,i]=function((i-ref_pixel)/self.wave)
        if dl==self.dltime:
            print "gau vect shape dl=====",vect.shape
        else:
            print "gau vect shape et dl",vect.shape,dl
	import  scipy.signal
	v=scipy.signal.slepian(dl,0.01)
	for i in range(dl):
		vect[:,:,i]=v[i]
        return vect
        
    
if __name__ == "__main__":
        ms3=sys.argv[1]
        MS=ClassMS(ms3,"UVW","DATA",10)
        print MS.na, MS.data.shape
        data,t=MS.NormalAvgALBL()
        print data.shape,t
        # #print MS.WeightForALBL(lambda x: np.exp(x)).shape
       #  data1=MS.ConvolveALBL(lambda x: np.exp((-x**2)/2.))
       #  print data.shape[0]
       #  print data1.shape[0]
       #  print "shape normL",np.std(data)
       #  print "shape gaussi",np.std(data1)
       # # print "shape weigth",np.std(data1)


       #  # # print "shape",data.shape
       #  # # print "sdt",np.std(data)
       #  # t_3s=table(ms,readonly=False,ack=False)
       #  # A1=t_3s.getcol("ANTENNA1") 
       #  # A2=t_3s.getcol("ANTENNA2")
       #  # vis1=data[(A1==0)&(A2==1)]
