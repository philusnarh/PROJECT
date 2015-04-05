from Pyxis.ModSupport import *

import imager,std

import pyfits
import Tigger
import im.argo

register_pyxis_module(superglobals="OUTFILE");

tigger_restore  = x("tigger-restore")
tigger_convert  = x("tigger-convert")
tigger_tag      = x("tigger-tag")

v.define("LSM","lsm.lsm.html",
  """current local sky model""");

define("LSM_TDL_Template","tiggerlsm.filename=$LSM",
  """TDL option for selecting current lsm""");
define("LSMREF","",
  """reference LSM (for transferring tags, etc.)""");
  
define('PYBDSM_OUTPUT_Template',"${OUTFILE}_pybdsm.lsm.html","""output LSM file for pybdsm search""");    
define('PYBDSM_POLARIZED',None,"""set to True to run pybdsm in polarized mode""");
define('PYBDSM_OPTIONS',{},"Extra options given to pybdsm");
_pybdsm = x.pybdsm;

define('CLUSTER_DIST',0,
  """source clustering distance, arcsec. If 0, then CLUSTER_DIST_BEAMS is used instead.""");
define('CLUSTER_DIST_BEAMS',3,
  """source clustering distance, in units of PSF size (measured as (BMAJ+BMIN)/2). If BMAJ/BMIN is not defined,
  this falls back to 60 arcsec.""");
define('MIN_EXTENT',0,
  """minimum Gaussian source extent, arcsec; sources smaller than this will be converted to point sources""");

def pybdsm_search (image="${imager.RESTORED_IMAGE}",output="$PYBDSM_OUTPUT",pol='$PYBDSM_POLARIZED',
  select=None,
  threshold=None,pbexp=None,**kw):
  """Runs pybdsm on the specified 'image', converts the results into a Tigger model and writes it to 'output'.
  Use 'threshold' to specify a non-default threshold (thresh_isl and thresh_pix).
  Use 'pol' to force non-default polarized mode.
  Use 'pbexp' to supply a primary beam expression (passed to tigger-convert), in which case the output model will contain
  intrinsic fluxes.
  Use 'select' to apply a selection string on the new model (e.g. "I.gt.0.001")
  """
  image,output,pol = interpolate_locals("image output pol");
  makedir(v.DESTDIR);
  # setup parameters
  gaul = II("${output:BASEPATH}.gaul");
  # info("PyBDSM filenames are $output $gaul");
  # start with default PYBDSM options
  opts = PYBDSM_OPTIONS.copy();
  opts.update(kw);
  # override with explicit arguments
  if threshold:
    opts['thresh_pix'] = threshold;
  if pol is not None:
    opts['polarisation_do'] = is_true(pol);
  pol = opts.get('polarisation_do',False);
  opts['quiet'] = True;
  # run pybdsm
  info("running PyBDSM process_image($image,%s)"%",".join(sorted([ "%s=%s"%x for x in opts.iteritems() ])));
  from lofar import bdsm
  img = bdsm.process_image(image,**kw);
  info("writing PyBDSM gaul catalog");
  img.write_catalog(outfile=gaul,format='ascii',catalog_type='gaul',clobber=True);
  # add log to output
  logfile = II("${output:BASEPATH}.pybdsm.log");
  if exists(logfile):
    info("PyBDSM log output follows:");
    for line in file(logfile):
      print "     ",line;
  else:
    warn("PyBDSM log $logfile not found");
  # set clustering parameter from beam size
  cluster = CLUSTER_DIST;
  if not cluster:
    hdr = pyfits.open(image)[0].header;
    # BMAJ/BMIN is in degrees -- convert to seconds, or fall back to 60" if not set
    cluster = 1800*(hdr.get('BMAJ',0)+hdr.get('BMIN',0))*CLUSTER_DIST_BEAMS or 60;
  # convert catalog
  if pbexp:
    args = [ "--primary-beam",pbexp,"--app-to-int" ]
  else:
    args = []
  if select:
    args += [ "--select",select ];
  verifyGaulModel(gaul)
  
  #Dictionary for establishing correspondence between parameter names in gaul files produced by pybdsm, and pyxis parameter names
  dict_gaul2lsm = {'Gaus_id':'name', 'Isl_id':'Isl_id', 'Source_id':'Source_id', 'Wave_id':'Wave_id', 'RA':'ra_d', 'E_RA':'E_RA', 'DEC':'dec_d', 'E_DEC':'E_DEC', 'Total_flux':'i', 'E_Total_flux':'E_Total_flux', 'Peak_flux':'Peak_flux', 'E_Peak_flux':'E_Peak_flux', 'Xposn':'Xposn', 'E_Xposn':'E_Xposn', 'Yposn':'Yposn', 'E_Yposn':'E_Yposn', 'Maj':'Maj', 'E_Maj':'E_Maj', 'Min':'Min', 'E_Min':'E_Min', 'PA':'PA', 'E_PA':'E_PA', 'Maj_img_plane':'Maj_img_plane', 'E_Maj_img_plane':'E_Maj_img_plane', 'Min_img_plane':'Min_img_plane', 'E_Min_img_plane':'E_Min_img_plane', 'PA_img_plane':'PA_img_plane', 'E_PA_img_plane':'E_PA_img_plane', 'DC_Maj':'emaj_d', 'E_DC_Maj':'E_DC_Maj', 'DC_Min':'emin_d', 'E_DC_Min':'E_DC_Min', 'DC_PA':'pa_d', 'E_DC_PA':'E_DC_PA', 'DC_Maj_img_plane':'DC_Maj_img_plane', 'E_DC_Maj_img_plane':'E_DC_Maj_img_plane', 'DC_Min_img_plane':'DC_Min_img_plane', 'E_DC_Min_img_plane':'E_DC_Min_img_plane', 'DC_PA_img_plane':'DC_PA_img_plane', 'E_DC_PA_img_plane':'E_DC_PA_img_plane', 'Isl_Total_flux':'Isl_Total_flux', 'E_Isl_Total_flux':'E_Isl_Total_flux', 'Isl_rms':'Isl_rms', 'Isl_mean':'Isl_mean', 'Resid_Isl_rms':'Resid_Isl_rms', 'Resid_Isl_mean':'Resid_Isl_mean', 'S_Code':'S_Code', 'Total_Q':'q', 'E_Total_Q':'E_Total_Q', 'Total_U':'u', 'E_Total_U':'E_Total_U', 'Total_V':'v', 'E_Total_V':'E_Total_V', 'Linear_Pol_frac':'Linear_Pol_frac', 'Elow_Linear_Pol_frac':'Elow_Linear_Pol_frac', 'Ehigh_Linear_Pol_frac':'Ehigh_Linear_Pol_frac', 'Circ_Pol_Frac':'Circ_Pol_Frac', 'Elow_Circ_Pol_Frac':'Elow_Circ_Pol_Frac', 'Ehigh_Circ_Pol_Frac':'Ehigh_Circ_Pol_Frac', 'Total_Pol_Frac':'Total_Pol_Frac', 'Elow_Total_Pol_Frac':'Elow_Total_Pol_Frac', 'Ehigh_Total_Pol_Frac':'Ehigh_Total_Pol_Frac', 'Linear_Pol_Ang':'Linear_Pol_Ang', 'E_Linear_Pol_Ang':'E_Linear_Pol_Ang'}

  #Dictionary for classifying a parameter as a general parameter or a polarization-specific parameter
  dict_pol_flag = {'Gaus_id':0, 'Isl_id':0, 'Source_id':0, 'Wave_id':0, 'RA':0, 'E_RA':0, 'DEC':0, 'E_DEC':0, 'Total_flux':0, 'E_Total_flux':0, 'Peak_flux':0, 'E_Peak_flux':0, 'Xposn':0, 'E_Xposn':0, 'Yposn':0, 'E_Yposn':0, 'Maj':0, 'E_Maj':0, 'Min':0, 'E_Min':0, 'PA':0, 'E_PA':0, 'Maj_img_plane':0, 'E_Maj_img_plane':0, 'Min_img_plane':0, 'E_Min_img_plane':0, 'PA_img_plane':0, 'E_PA_img_plane':0, 'DC_Maj':0, 'E_DC_Maj':0, 'DC_Min':0, 'E_DC_Min':0, 'DC_PA':0, 'E_DC_PA':0, 'DC_Maj_img_plane':0, 'E_DC_Maj_img_plane':0, 'DC_Min_img_plane':0, 'E_DC_Min_img_plane':0, 'DC_PA_img_plane':0, 'E_DC_PA_img_plane':0, 'Isl_Total_flux':0, 'E_Isl_Total_flux':0, 'Isl_rms':0, 'Isl_mean':0, 'Resid_Isl_rms':0, 'Resid_Isl_mean':0, 'S_Code':0, 'Total_Q':1, 'E_Total_Q':1, 'Total_U':1, 'E_Total_U':1, 'Total_V':1, 'E_Total_V':1, 'Linear_Pol_frac':1, 'Elow_Linear_Pol_frac':1, 'Ehigh_Linear_Pol_frac':1, 'Circ_Pol_Frac':1, 'Elow_Circ_Pol_Frac':1, 'Ehigh_Circ_Pol_Frac':1, 'Total_Pol_Frac':1, 'Elow_Total_Pol_Frac':1, 'Ehigh_Total_Pol_Frac':1, 'Linear_Pol_Ang':1, 'E_Linear_Pol_Ang':1}

  lines = [line.strip() for line in open(gaul)]
  
  for line in range(len(lines)):
    if lines[line]:
      if lines[line].split()[0] is not '#': 
        gaul_params = lines[line-1].split()[1:] #Parameter list is last line in gaul file that begins with a '#'
        break
  
  # Initialize lists for general and polarization parameters 
  lsm_params_general = []
  lsm_params_polarization = []

  for param in gaul_params:
    if dict_pol_flag[param] is 0:
     lsm_params_general.append(dict_gaul2lsm[param])
    if dict_pol_flag[param] is 1:
     lsm_params_polarization.append(dict_gaul2lsm[param])
  
  general_params_string = ' '.join(lsm_params_general)
  pol_params_string = ' '.join(lsm_params_polarization)

  tigger_convert(gaul,output,"-t","ASCII","--format", general_params_string + (pol_params_string if pol else ""),
    "-f","--rename",
    "--cluster-dist",cluster,
    "--min-extent",MIN_EXTENT,
    split_args=False,
    *args);
    
document_globals(pybdsm_search,"PYBDSM_* imager.RESTORED_IMAGE CLUSTER_* MIN_EXTENT");

def verifyGaulModel(gaullsm):
  """Check all sources in a gaul file are in valid locations before running tigger
  convert. Useful when images are 'all-sky' and have undefined regions.
  """
  falseSources=0
  olsm=''
  fh=open(gaullsm,'r')
  for ll in fh.readlines():
    cll=' '.join(ll.split())
    if cll=='' or cll.startswith('#'):
      olsm+=ll
      continue
    lineArray=cll.split(' ')
    if math.isnan(float(lineArray[4])): falseSources+=1
    else: olsm+=ll
  fh.close()

  fh=open(gaullsm,'w')
  fh.write(olsm)
  fh.close()


def transfer_tags (fromlsm="$LSMREF",lsm="$LSM",output="$LSM",tags="dE",tolerance=60*ARCSEC):
  """Transfers tags from a reference LSM to the given LSM. That is, for every tag
  in the given list, finds all sources with those tags in 'fromlsm', then applies 
  these tags to all nearby sources in 'lsm' (within a radius of 'tolerance'). 
  Saves the result to an LSM file given by 'output'.
  """
  fromlsm,lsm,output,tags = interpolate_locals("fromlsm lsm output tags");
  # now, set dE tags on sources
  tagset = frozenset(tags.split());
  info("Transferring tags %s from %s to %s (%.2f\" tolerance)"%(",".join(tagset),fromlsm,lsm,tolerance/ARCSEC));
  import Tigger
  refmodel = Tigger.load(fromlsm);
  model = Tigger.load(lsm);
  # for each dE-tagged source in the reference model, find all nearby sources
  # in our LSM, and tag them
  for src0 in refmodel.getSourceSubset(",".join(["="+x for x in tagset])):
    for src in model.getSourcesNear(src0.pos.ra,src0.pos.dec,tolerance=tolerance):
      for tag in tagset:
        tagval = src0.getTag(tag,None);
        if tagval is not None:
          if src.getTag(tag,None) != tagval:
            src.setTag(tag,tagval);
            info("setting tag %s=%s on source %s (from reference source %s)"%(tag,tagval,src.name,src0.name))
  model.save(output);

  

CC_RESCALE = 1.  
CC_IMAGE_Template = "${OUTFILE}_ccmodel.fits"

def add_ccs (lsm="$LSM",filename="${imager.MODEL_IMAGE}",
             cc_image="$CC_IMAGE",srcname="ccmodel",output="$LSM",zeroneg=True,scale=None,pad=1):
  """Adds clean components from the specified FITS image 'filename' to the sky model given by 'lsm'.
  Saves the result to an LSM file given by 'output'.
  The CC image is copied to 'cc_image', optionally rescaled by 'scale', and optionally has negative pixels reset to zero (if zeroneg=True).
  'srcname' gives the name of the resulting LSM component.
  'pad' gives the padding attribute of the LSM component, use e.g. 2 if CC image has significant signal towards the edges.
  """;
  lsm,filename,cc_image,srcname,output = interpolate_locals("lsm filename cc_image srcname output");
  info("adding clean components from $filename ($cc_image), resulting in model $output");
  # rescale image
  ff = pyfits.open(filename);
  ff[0].data *= (scale if scale is not None else CC_RESCALE);
  if zeroneg:
    ff[0].data[ff[0].data<0] = 0;
  ff.writeto(cc_image,clobber=True);
  tigger_convert(lsm,output,"-f","--add-brick","$srcname:$cc_image:%f"%pad);

document_globals(add_ccs,"MODEL_CC_*");  


def pointify (lsm="$LSM",output="$LSM",name=""):
  """Replaces names sources with point sources""";
  lsm,output,name = interpolate_locals("lsm output name");
  model = Tigger.load(lsm);
  src = model.findSource(name);
  info("Setting source $name in model $lsm to point source, saving to $output");
  src.shape = None;
  model.save(output);

document_globals(add_ccs,"MODEL_CC_*");  


define('SOFIA_CFG','${OUTFILE}_sofia_conf.txt','Name of auto generated sofia config file')
define('SOFIA_PATH','sofia_pipeline.py','SoFia pipeline path (sofia_pipeline.py)')
define('SOFIA_OUTDIR','${v.DESTDIR}','SoFiA output directory')
define('SOFIA_OUT_PREFIX','${OUTFILE:BASE}sofia','Prefix for sofia outputs')
define('_SOFIA_DEFAULTS',{},'Dictionary of default SoFiA options')

_SOFIA_DEFAULTS = {'writeCat': {'writeASCII': 'true', 'parameters': "['*']", 'compress': 'false', 'basename': '', 'writeXML': 'false', 'writeSQL': 'false', 'outputDir': ''}, 'parameters': {'dilateThreshold': '0.02', 'optimiseMask': 'false', 'fitBusyFunction': 'false', 'dilatePixMax': '10', 'dilateChan': '1', 'dilateMask': 'false'}, 'smooth': {'kernelZ': '3.0', 'kernel': 'gaussian', 'kernelX': '3.0', 'kernelY': '3.0', 'edgeMode': 'constant'}, 'merge': {'radiusY': '3', 'radiusX': '3', 'radiusZ': '3', 'minSizeY': '3', 'minSizeX': '3', 'minSizeZ': '2'}, 'flag': {'regions': '[]'}, 'optical': {'storeMultiCat': 'false', 'sourceCatalogue': '', 'specSize': '1e+5', 'spatSize': '0.01'}, 'steps': {'doReliability': 'false', 'doWriteFilteredCube': 'false', 'doMom1': 'false', 'doParameterise': 'true', 'doFlag': 'false', 'doSmooth': 'false', 'doDebug': 'false', 'doCNHI': 'false', 'doWavelet': 'false', 'doCubelets': 'false', 'doWriteMask': 'false', 'doOptical': 'false', 'doSubcube': 'false', 'doMerge': 'true', 'doMom0': 'false', 'doScaleNoise': 'false', 'doSCfind': 'true', 'doWriteCat': 'true', 'doThreshold': 'false'}, 'wavelet': {'threshold': '5.0', 'scaleXY': '-1', 'positivity': 'false', 'iterations': '3', 'scaleZ': '-1'}, 'threshold': {'threshold': '4.0', 'rmsMode': 'std', 'clipMethod': 'relative', 'verbose': 'false'}, 'import': {'weightsFunction': '', 'subcubeMode': 'pixel', 'subcube': '[]', 'maskFile': '', 'inFile': '', 'weightsFile': ''}, 'CNHI': {'verbose': '1', 'qReq': '3.8', 'maxScale': '-1', 'medianTest': 'true', 'minScale': '5', 'pReq': '1e-5'}, 'reliability': {'threshold': '0.9', 'kernel': '[0.15,0.05,0.1]', 'parSpace': "['SNRsum','SNRmax','NRvox']", 'fMin': '0.0', 'makePlot': 'false'}, 'scaleNoise': {'scaleX': 'false', 'scaleY': 'false', 'scaleZ': 'true', 'statistic': 'mad', 'edgeX': '0', 'edgeY': '0', 'edgeZ': '0'}, 'SCfind': {'kernels': "[[ 0, 0, 0,'b'],[ 0, 0, 3,'b'],[ 0, 0, 7,'b'],[ 0, 0, 15,'b'],[ 3, 3, 0,'b'],[ 3, 3, 3,'b'],[ 3, 3, 7,'b'],[ 3, 3, 15,'b'],[ 6, 6, 0,'b'],[ 6, 6, 3,'b'],[ 6, 6, 7,'b'],[ 6, 6, 15,'b']]", 'maskScaleXY': '2.0', 'verbose': 'true', 'rmsMode': 'negative', 'edgeMode': 'constant', 'kernelUnit': 'pixel', 'threshold': '4.0', 'maskScaleZ': '2.0', 'sizeFilter': '0.0'}}

def sofia_search(fitsname='${im.RESTORED_IMAGE}',sofia_conf=None,
                 threshold=4,do_reliability=True,reliability=0.9,
                 merge=3,basename='$SOFIA_OUT_PREFIX',makeplot=True,
                 outdir='$SOFIA_OUTDIR',options={}):
    """ Runs SoFiA source finding pipeline. Only a few options are provided here. 
        For more eleborate settings, add options (as you would in a SoFiA configuarion file) 
        via the [options] dictionary or provide a SoFiA configuration file via [sofia_conf]
       -------
       fitsname : Name of fits map on which to run the source finder
       sofia_conf : SoFiA configuration file
       threshold : Peak threshold for source finder [in units of sigma above noise rms]
       do_reliability : Do reliability caltulations
       reliability : Reliability threshold. (0,1]
       merge : merge 
       options : extra options which will directly to sofia configuration file.
    """

    fitsname,basename,outdir = interpolate_locals('fitsname basename outdir')
    makedir(v.DESTDIR);

    # swap freq and stokes in fits hdr
    im.argo.swap_stokes_freq(fitsname,freq2stokes=True)

    # use default SoFiA config if not specified
    if not sofia_conf:

        sofia_conf = II(SOFIA_CFG)
        _SOFIA_DEFAULTS['import']['inFile'] = fitsname
        _SOFIA_DEFAULTS['writeCat']['outputDir'] = outdir
        _SOFIA_DEFAULTS['writeCat']['basename'] = basename

        if os.path.exists(sofia_conf):
            x.sh(II('mv $sofia_conf ${sofia_conf}.old'))
        sofia_std = open(sofia_conf,'w')

        if threshold!=4 and isinstance(threshold,(int,float)):
            _SOFIA_DEFAULTS['threshold']['threshold'] = threshold
        if do_reliability and isinstance(do_reliability,bool):
            _SOFIA_DEFAULTS['steps']['doReliability'] = 'true'
            if isinstance(reliability,(int,float)):
                _SOFIA_DEFAULTS['reliability']['threshold'] = reliability
                _SOFIA_DEFAULTS['reliability']['makePlot'] = makeplot
        if merge!=3 and isinstance(merge,int):
            _SOFIA_DEFAULTS['merge']['mergeX'] = merge
            _SOFIA_DEFAULTS['merge']['mergeY'] = merge
            _SOFIA_DEFAULTS['merge']['mergeZ'] = merge

      # Update default SoFia configuration dictionary with user options
        for key,val in options.iteritems():
            a,b = key.split('.')
            if a not in _SOFIA_DEFAULTS.keys():
                abort('Option ${a}.${b} is not recognisable.')
            _SOFIA_DEFAULTS[a][b] = val

      # Generate sofia configuration file
        sofia_std.write('#Sofia autogen configuration file [pyxis]')
        for a,b in  _SOFIA_DEFAULTS.iteritems():
            for key,val in b.iteritems():
                sofia_std.write('\n%s.%s = %s'%(a,key,val))
        sofia_std.close()

    x.sh('$SOFIA_PATH $sofia_conf')

document_globals(sofia_search," im.RESTORED_IMAGE _SOFIA_DEFAULTS SOFIA_CFG SOFIA_OUTDIR SOFIA_PATH SOFIA_OUT_PREFIX")


