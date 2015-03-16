
from Pyxis.ModSupport import *
from optparse import OptionParser
from pyrap.tables import table
import ClusteredSources
import pickle
import Tigger
import mqt,cal,lsm,imager
import os,sys
import pyfits
import ms_extract
import GenerateRandomFields
import AddingTag

# ouptut files
v.DESTDIR_Template = "${OUTDIR>/}testplots-${MS:BASE}"
v.OUTFILE_Template = "${DESTDIR>/}${MS:BASE}"
v.LSM_Template = '${v.OUTFILE}.lsm.txt'
v.LOG_Template = "${v.OUTFILE}.log"
# number of meq threads
mqt.MULTITHREAD=4

# imaging options
imager.npix = 2048
imager.cellsize = "4arcsec"
imager.mode = "channel"
imager.stokes = "I"
imager.weight = "natural"
imager.filter = "1arcsec,1arcsec,0deg"
imager.wprojplanes = 0
# cleaning options:
imager.niter = 100000
imager.gain = .1
imager.threshold ="0.1mJy"


def SortingModel(LSM):

     """
       Sorting the skymodel in descending order of brightness
     """
     f_skymodel="skymodel.txt"

     try:
        ClusteredSources.Cluster(LSM,f_skymodel)
        skymodel_arr = np.loadtxt(f_skymodel,dtype='str')
     except AttributeError:
        print " AttributeError: Source object has no attribute cluster"
        if LSM.endswith(".lsm.html"): LSM = txtconvert(LSM)
        skymodel_arr = np.loadtxt(LSM,dtype='str')

     skymodel_sortedarr = sorted(skymodel_arr, key=lambda skymodel_arr_entry: skymodel_arr_entry[3])
     sorted_skyarr = np.copy(skymodel_sortedarr)
     numpyconvert(sorted_skyarr,f_skymodel)
     return f_skymodel, sorted_skyarr



def simulate(msname="$MS",LSM="$LSM",section="$SIM_CONFIG_SECTION",args=[],options={},**kw):
     """
        Simulating our skymodel
     """
     # necessary to make '$MS' etc. above work properly
     msname,section = interpolate_locals("msname section");
     options.update(kw);
     # setup TDL options for MS, LSM, etc., borrowing templates from the cal module
     args = [ "${cal.MS_TDL} ${cal.CHAN_TDL} ${cal.LSM_TDL}" ] + list(args);

     # run meqtrees
     print "LSM=",LSM
     f_skymodel,sorted_skyarr=SortingModel(LSM)
     options['tiggerlsm.filename']=f_skymodel
     options['ms_sel.msname']=MS
     options['noise_stddev']=0.001
     mqt.run("turbo-sim.py",job="_tdl_job_1_simulate_MS",config="tdlconf.profiles",section=section,options=options)
     dirty = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1)
     if section=="empty_field":
         imager.make_image(column='CORRECTED_DATA', dirty=dirty, restore=False, restore_lsm=False, algorithm='csclean')
         empty_field=imager.DIRTY_IMAGE
         return empty_field
     



def calibrate(msname="$MS",LSM="$LSM",section="$SIM_CONFIG_SECTION",percentage=None,smoothtimeinterval=None,outfile=None,gaulfile=None,imgname=None,**kw):
     """
       Calibrating with the calibration skymodel
     """
     
     # creating an empty image
     empty_field=simulate("$MS",LSM=LSM,section="empty_field")

     rad2deg = lambda val: val * 180./np.pi    # conversion from radians to degrees
     f_skymodel,sorted_skyarr=SortingModel(LSM)
     nsources = len(sorted_skyarr)
     unmodelsrc = percentage
     unmodel_skyarr = sorted_skyarr[0:unmodelsrc]
     cal_skyarr = sorted_skyarr[unmodelsrc::,:]
     if section=="dE_cal":
       if len(sorted_skyarr)==2:
          phase_center = phasecenter(MS)
          newrow=np.array(['A2',str(rad2deg(phase_center[0])-0.5),str(rad2deg(phase_center[1])-0.4),'1e-20','+cal'])
          cal_skyarr = np.vstack([cal_skyarr,newrow])

     # simulating the visibilities
     simulate("$MS",LSM=LSM,section="turbo-sim:sim")

     # writing numpy array to textfile
     unmodel_skymodel = "unmodel_skymodel.txt"
     cal_skymodel = "cal_skymodel.txt"
     numpyconvert(cal_skyarr, cal_skymodel),numpyconvert(unmodel_skyarr, unmodel_skymodel)

     #converting txt file to lsm.html
     lsmskymodel =lsmconvert(f_skymodel)
     cal_lsmskymodel=lsmconvert(cal_skymodel)
     unmodel_lsmskymodel=lsmconvert(unmodel_skymodel)
     


     # calibrating with skymodel
     options={}
     options['tiggerlsm.filename']=cal_lsmskymodel
     options['ms_sel.ms_ifr_subset_str']="all"
     options['ms_sel.msname']=MS
     options['stefcal_diffgain.timesmooth']=smoothtimeinterval
     mqt.run("calico-stefcal.py",job="stefcal",config="tdlconf.profiles",section=section,args=args,options=options)
   
     # imaging the calibrated visibilities
     dirty = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1)
     restored = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1)
     #imager.make_image(column='CORRECTED_DATA', dirty=dirty, restore=restored, restore_lsm=False, algorithm='csclean', lsm='$LSM')

     # cleaning using a mask#    
     maskIm = MS + "_masked.fits"
     x.sh("tigger-restore %s %s %s -b 20.0'' -f" %(empty_field,unmodel_lsmskymodel,maskIm))
     peak_flux = pyfits.open(maskIm)[0].data[...].max()
     rms = pyfits.open(maskIm)[0].data[...].std()
     imager.make_threshold_mask(input=maskIm,threshold=0.5*rms)
     imager.threshold = str((float(unmodel_skyarr[len(unmodel_skyarr)-1,3]))/1000)+"Jy"
     imager.make_image(column='CORRECTED_DATA',dirty=dirty,restore=restored,restore_lsm=False,algorithm='csclean',mask=True)
     

     # Importing the extracte flux as a string
     calibrationflux = (TotalFlux(unmodel_skyarr,3)/TotalFlux(sorted_skyarr,3))*100
     str_out ="#   ****** %s Sources in skymodel ******\n#   ***** %s Sources in calibration model *****\n#   ****** %.2g%s flux removed from calibration model ******\n" %(str(nsources),str(unmodelsrc),calibrationflux,'%')
     
     # Importing the residual stats
     str_out2 ,rms= fitsinfo(imager.RESIDUAL_IMAGE)

     str_out += SourceDetector(imager.RESTORED_IMAGE,unmodel_lsmskymodel,5,3,rms,gaulfile)
     
     str_out +=str_out2
     # Witing the extrated stats to txtfile
     Writetotxt(str_out,outfile)
     # Saving the fits image
     SaveMap(imager.RESTORED_IMAGE, imgname)



      
def SaveMap(fitsname, newfitsname):
    """ Saving the fits image """
   
    os.system("mv "+imager.RESTORED_IMAGE+" "+newfitsname)


def MaskImaging(maskIm, emptyfield,skymodel):
     x.sh("tigger-restore %s %s %s -b 90.0'' -f" %(empty_field,skymodel,maskIm))
     skyarr = txtconvert(skymodel)
     peak_flux = pyfits.open(maskIm)[0].data[...].max()
     rms = pyfits.open(maskIm)[0].data[...].std()
     imager.make_threshold_mask(input=maskIm,threshold=0.5*rms)
     #imager.threshold = 0.1*str((float(skyarr[len(skyarr)-1,3]))/1000)+"Jy"
     imager.make_image(column='CORRECTED_DATA',dirty=dirty,restore=restored,restore_lsm=False,algorithm='csclean',mask=True)


def Extractflux(restored_image,thresh_pix,thresh_isl,threshold,gaulfile):
        """                                   
            Extracting flux using pybdsm   
        """                                  
        
        lsm.pybdsm_search(image=imager.RESTORED_IMAGE,thresh_pix=thresh_pix,thresh_isl=thresh_isl,blank_limit=threshold,gaul=gaulfile,clobber=True)

        pybdsmsky = lsm.PYBDSM_OUTPUT
        pybdsmmodtxtfile = txtconvert(pybdsmsky)
        pybdsmarray = np.copy(np.loadtxt(pybdsmmodtxtfile,dtype='str'))
        return pybdsmsky,pybdsmarray


def SourceDetector(restored_image,unmodelsky,thresh_pix,thresh_isl,threshold,gaulfile):
     """
         Picking the gaussian fit closest to the actual source
     """
     

     deg2rad = lambda val: val/60.* np.pi/180.
     rad2deg = lambda val: val * 180./np.pi
     unmodeltxt = txtconvert(unmodelsky)
     unmodelarr = np.loadtxt(unmodeltxt,dtype='str') 
     unmodelmod = Tigger.load(unmodelsky)

     
     pybdsmsky,pybdsmarr = Extractflux(restored_image,thresh_pix,thresh_isl,threshold,gaulfile)   
     pybdsmmod = Tigger.load(pybdsmsky)
     
          
     
     # extracting sources above a given sigma level
     str_out = "#   name actualflux retrievedflux suppression \n"

     # getting the sources within a certain tolerance value
     str_out = "#   name pybdsmname ra_d dec_d actualflux retrievedflux E_retrievedflux suppression \n"
     gaullist = np.loadtxt(gaulfile,dtype='str')
     for src in unmodelmod.sources:
         fluxarray = np.array([])
         name,ra, dec ,amp = src.name,src.pos.ra,src.pos.dec, src.flux.I
         sources=pybdsmmod.getSourcesNear(ra,dec,deg2rad(1))
         if sources!=[]:
            for nsrc in sources:
                clustername,clusterflux = nsrc.cluster, nsrc.cluster_flux
                fluxarray = np.append(fluxarray,clusterflux)
                retrievedflux = np.max(fluxarray)
                nsrc.setTag('real',1)
                src_error=np.array([])
                for k in range(len(gaullist)):
                  if gaullist.ndim !=1:
                     if float(gaullist[k,8])==float(nsrc.flux.I):
                         src_error = np.append(src_error,float(gaullist[k,9]))
                  else:
                     if float(gaullist[8])==float(nsrc.flux.I):
                         src_error = np.append(src_error,float(gaullist[9]))

                error_val = np.sum(src_error)/len(src_error)                             
            

            suppression = (amp-retrievedflux)/amp
            E_suppression = error_val*suppression/retrievedflux
            str_out +="   %s %s %.6g %.6g %.4g %.4g %.4g %.4g \n" %(name,clustername,rad2deg(ra),rad2deg(dec),amp,retrievedflux,error_val,suppression)
     
     
     
     newmodel="pybdsm_"+str(len(unmodelmod.sources))+".lsm.html"
     #Taglsm(pybdsmsky,newmodel,'real')
     return str_out


def Writetotxt(str_out,filename):
    """
      Writing to text file
    """
    std_out = open(filename,'w')
    std_out.write(str_out)
    std_out.close()




def TotalFlux(array,index):
     """
       Computing sum of a specific row
     """
     # converting string to float
     floatarray = np.array(array[:,index],dtype='float')
     totalsum = np.sum(floatarray)
     
     return totalsum

def fitsinfo(fitsname):
        """
        Get fits info
        """
        stats = "\n#   INFO: getting stats from residual image \n#   sum min max mean rms \n"
        fitsdata = pyfits.open(fitsname)[0].data[...]
        total = fitsdata.sum()
        min_flux = fitsdata.min()
        max_flux = fitsdata.max()
        mean_flux = fitsdata.mean()
        rms = fitsdata.std()
        stats += "#   %.4g %.4g %.4g %.4g %.4g \n"%(total,min_flux,max_flux,mean_flux,rms)
        return stats,rms



def lsmconvert(skymodel):
          """
             Converting textfile to lsm.html format
          """
          os.system("tigger-convert %s -f"%(skymodel))
          return "%s.lsm.html"%(skymodel[:-4])

def txtconvert(skymodel):
          """
             Converting lsm.html to txt
          """

          txtfile = skymodel[:-9] + ".txt"
          os.system("tigger-convert %s %s -f"%(skymodel,txtfile))
          return txtfile

def Taglsm(skymodel=None,newmodel=None,tag=None):
          """
             Converting textfile to lsm.html format
          """
          os.system("tigger-convert %s %s --select %s==1 -f"%(skymodel,newmodel,tag))
          return newmodel


def numpyconvert(array,txtfile):
      """
         Converting numpy array to txtfile
      """

      np.savetxt(open(txtfile,"w"), array, fmt='%s')
      with open(txtfile, "r+") as f: s = f.read(); f.seek(0); f.write("#format: name ra_d dec_d i tags...\n" + s)
      return txtfile

def readfirstline(textfile):
      """
        Reading forst line of a textfile
      """
      with open("skymodel.txt","r") as f:
          first_line=f.readline()
      return first_line


def extract(MS):
    tl=table(MS,readonly=False)
    corr_data=tl.getcol("CORRECTED_DATA")
    model_data=tl.getcol("MODEL_DATA")
    data = tl.getcol("DATA")
    tl.close()
    return data, corr_data, model_data


def TagSources(skymodel='$LSM',outfile=None,tag=None):
    """
      Tracking the sources with a specific tag
    """
    print LSM
    if skymodel.endswith(".txt"): skymodel = lsmconvert(II(skymodel)) # converting to lsm.html format
    f = open(outfile,'w')
    mod = Tigger.load(skymodel)
    for src in mod.sources:
        if src.getTag('%s'%str(tag))==True: f.write("%s\n"%src.name)


def RemovingTag(skymodel='$LSM',dElsm=None,index=None):
    """ 
      Removing Individual Sources with 'dE' Tag
    """
    mod = Tigger.load(II(skymodel))
    dEmod = Tigger.load(dElsm)
    dESources = [src2 for src2 in dEmod.sources]
    if index == 0:
        return skymodel
    else:
        dESources = dESources[0:index]
        for src2 in dESources:
           for src1 in mod.sources:
             if src1.name == src2.name:
                 src1.setTag('dE',None)
        outfile = "new"+II(skymodel)
        mod.save(outfile) 
        return outfile


def phasecenter(MS):
    tab = table(MS+"/FIELD")
    phase_centre = (tab.getcol("PHASE_DIR"))[0,0,:]
    tab.close()
    return phase_centre

def write(MS,data,column):
    tl=table(MS, readonly=False)
    tl.putcol(column ,data)
    tl.close()


def lm2raddec(msname='$MS',l=None,m=None):
    rad2deg = lambda val: val * 180./np.pi
    ra0,dec0 = phasecenter(II(msname)) # phase centre in radians
    rho = np.sqrt(l**2+m**2)
    if rho==0:
       ra = ra0
       dec = dec0
    else:
       cc = np.arcsin(rho)
       ra = ra0 - np.arctan2(l*np.sin(cc), rho*np.cos(dec0)*np.cos(cc)-m*np.sin(dec0)*np.sin(cc))
       dec = np.arcsin(np.cos(cc)*np.sin(dec0) + m*np.sin(cc)*np.cos(dec0)/rho)
    return rad2deg(ra), rad2deg(dec)


def meqskymodel(msname='$MS',point_sources=None,outfile=None):
     str_out = "#format: name ra_d dec_d i tags... \n"
     for i in range(len(point_sources)):
          amp, l ,m = point_sources[i,0], point_sources[i,1], point_sources[i,2]
          ra_d, dec_d = lm2raddec(II(msname),l,m)
          name = "A"+ str(i)
          tags = "+dE" if amp>=1 else "+cal"
          str_out += "%s %.12g %.12g %.5g %s\n"%(name, ra_d, dec_d,amp,tags)

     skymodel = open(outfile,"w")
     skymodel.write(str_out)
     return outfile


def ExtractAntenna(MS):
    tl=table(MS,readonly=False)
    A1=tl.getcol("ANTENNA1")
    A2=tl.getcol("ANTENNA2")
    ta=table(MS +"/ANTENNA")
    na=len(ta.getcol("NAME")) #number of antennas
    nb=na*(na-1)/2+na # number of baselines
    ns=(len(A1))/nb # timeslots
    tl.close()
    ta.close()
    return na ,nb, ns


def ExtractdE(msname='$MS',calmodel=None):
    """
      Extracting Differential Gains
    """

    import cPickle
    from Timba.Meq import meq
    dgains = cPickle.load(file('diffgain.cp'))
    na, nb, ns = ExtractAntenna(II(msname))
    temp = np.ones((na,na),)
    dE = np.zeros((na,ns),)
    dEsoln = np.zeros((na,na,ns),)
    dEcalmodel=lsmdEconvert(calmodel)
    mod = Tigger.load(dEcalmodel)
    dEsources = mod.sources
    antenna_str="0123456789ABCD"
    for src in mod.sources:
        for s in range(len(antenna_str)):
           dEname = 'dE:'+str(src.name)
           dEsolnx = dgains['gains'][str(dEname)]['solutions'][(str(antenna_str[s]),0)]
           dEsolny = dgains['gains'][str(dEname)]['solutions'][(str(antenna_str[s]),1)]
           dE[s,:] = (dEsolnx[:,0]+dEsolny[:,0])/2

    for t in range(ns):
          dEsoln[:,:,t] = np.dot(np.diag(dE[:,t]),temp)
          dEsoln[:,:,t] = np.dot(dEsoln[:,:,t] ,np.diag(dE[:,t].conj()))
    return dEsoln, dEsources


def ImagingdE(msname='$MS',R=None,base=None,column=None,model=None,calmodel=None,unmodel=None):
       """
          Convolving differential gains with the visibilities
       """
       rad2deg = lambda val: val * 180./np.pi    # conversion from radians to degrees
       data ,corr_data, model_data = extract(II(msname))
       na ,nb, ns = ExtractAntenna(II(msname))
       dE, dEsources = ExtractdE(msname,calmodel)
       print dEsources
       print "dE=",dE
       mod = Tigger.load(model)
       print mod.sources
       for src in mod.sources:
             str_out = "#format: name ra_d dec_d i \n"
             name,ra,dec,amp = src.name,src.pos.ra,src.pos.dec,src.flux.I
             str_out += "%s %.8g %.8g %.4g \n"%(name,rad2deg(ra),rad2deg(dec),amp)
             with open('onesrcsky.txt', 'w') as f:f.write(str_out)
             if src.name=="A0":
                    dE_vis,base=sim_modellledsrc("$MS",section="simulation",model='onesrcsky.txt')

             else:               
                    vis,base=sim_modellledsrc("$MS",section="simulation",model='onesrcsky.txt')



             if column == "CORRECTED_DATA":
                 t_data = corr_data
             else:
                 t_data = data
             t_data.fill(0.)

       for t in range(ns):
           for rows in range(dE.shape[0]):
                for col in range(dE.shape[0]):
                       t_data[base[rows,col,t],0,3]=dE[rows,col,t]
                       t_data[base[rows,col,t],0,0]=dE[rows,col,t]
       write(II(msname),t_data,column)
       S, base = ms_extract.read_R(II(msname),"CORRECTED_DATA")
       print "S=",S
       dirty = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1)
       restored = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1)
       imager.make_image(column='%s'%(column), dirty=dirty, restore=False, restore_lsm=False, algorithm='csclean')


def convertstr(antenna):
      """
        Converting the antennas bearing an alphabetical string label to integer string
      """
      if antenna=="A":
         ant = '10'
      elif antenna=="B":
         ant = '11'
      elif antenna=="C":
         ant = '12'
      elif antenna=="D":
         ant = '13'
      else:
         ant = str(antenna)
      return ant


def Imagingperbaseline(column,baseline):
      """"
        Imaging along specific baseline
      """

      if baseline == "all":
         antenna_str = ""
      else:
          antennas = baseline.split(':') or baseline.split('-')
          antenna1,antenna2=convertstr(baseline[0]),convertstr(baseline[1])
          antenna_str ="(ANTENNA1 = "+antenna1+" && ANTENNA2 = "+antenna2+")"

      dirty = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1,select='%s'%(antenna_str))
      restored = dict(nchan=1,chanstart=0, chanstep=1, img_nchan=1, img_chanstart=0, img_chanstep=1,select='%s'%(antenna_str))
      imager.make_image(column='%s'%(column), dirty=dirty, restore=False, restore_lsm=False, algorithm='csclean')



def runall():
     nsources = 100
     percentage = np.arange(5,96,5)
     nsamples=100
     solnint = np.arange(1,62,10)
     
     
#     dELsm = Taglsm("clustered3C147.lsm.html",newmodel="dESources.lsm.html",tag="dE")
#     dEmod = Tigger.load(dELsm)
#     dESources = [srcs for srcs in dEmod]
#     dESourceNum = len(dESources)
#     print dESourceNum
#     for ind, src in enumerate(dESources):
#         LSM = RemovingTag("clustered3C147.lsm.html",dElsm=dELsm,index=ind)
#         for p in percentage:
#            for s in solnint:
#                outfile = "3C147dE_%s-%sint-%sdE.txt"%(str(p),str(s),str(dESourceNum-ind))
#                gaulfile = "3C147dE_%s-%sint-%sdE.gaul"%(str(p),str(s),str(dESourceNum-ind))
#                imgname = "3C147dE_%s-%sint-%sdE.fits"%(str(p),str(s),str(dESourceNum-ind))
#                unmodsrc = int((p/100.)*173)
#                calibrate('$MS',LSM=LSM,section="dE_cal",outfile=outfile, percentage=unmodsrc,smoothtimeinterval=s,gaulfile=gaulfile,imgname=imgname)
          


     for k in range(2,nsamples):
        for p in percentage:
            LSM = "200source_%s.lsm.html"%str(k+1)
            outfile = "200GsourceExp%s-natural_%s.txt"%(str(k+1),str(p))
            gaulfile = "200GsourceExp%s-natural_%s.gaul"%(str(k+1),str(p))
            imgname = "200GsourceExp%s-natural_%s.fits"%(str(k+1),str(p))
            unmodsrc = int((p/100.)*nsources)
            calibrate('$MS',LSM=LSM,section="calico-stefcal:GCal",outfile=outfile, percentage=unmodsrc,smoothtimeinterval=1,gaulfile=gaulfile,imgname=imgname)
 
          
   






