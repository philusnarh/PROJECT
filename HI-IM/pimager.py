##


# filenames for images 
define("BASENAME_IMAGE_Template","${OUTFILE}${-<IMAGER}","default base name for all image filenames below");
define("DIRTY_IMAGE_Template", "${BASENAME_IMAGE}.dirty.fits","output filename for dirty image");
define("PSF_IMAGE_Template", "${BASENAME_IMAGE}.psf.fits","output filename for psf image");
define("RESTORED_IMAGE_Template", "${BASENAME_IMAGE}.restored.fits","output filename for restored image");
define("RESIDUAL_IMAGE_Template", "${BASENAME_IMAGE}.residual.fits","output filename for deconvolution residuals");
define("MODEL_IMAGE_Template", "${BASENAME_IMAGE}.model.fits","output filename for deconvolution model");
define("FULLREST_IMAGE_Template", "${BASENAME_IMAGE}.fullrest.fits","output filename for LSM-restored image");
define("MASK_IMAGE_Template", "${BASENAME_IMAGE}.mask.fits","output filename for CLEAN mask");

# How to channelize the output image. 0 for average all, 1 to include all, 2 to average with a step of 2, etc.
# None means defer to 'imager' module options
define("IMAGE_CHANNELIZE",0,"image channels selection: 0 for all, 1 for per-channel cube")
# passed to tigger-restore when restoring models into images. Use e.g. "-b 45" for a 45" restoring beam.
define("RESTORING_OPTIONS","","extra options to tigger-restore for LSM-restoring")
# default clean algorithm
define("CLEAN_ALGORITHM","clark","CLEAN algorithm (clark, hogbom, csclean, etc.)")

def fits2casa (input,output):
    """Converts FITS image to CASA image."""
    im.argo.fits2casa(input,output)

def make_image (msname="$MS",column="$COLUMN",imager='$IMAGER',
                dirty=True,restore=False,restore_lsm=True,psf=False,
                dirty_image="$DIRTY_IMAGE",
                restored_image="$RESTORED_IMAGE",
                residual_image="$RESIDUAL_IMAGE",
                psf_image="$PSF_IMAGE",
                model_image="$MODEL_IMAGE",
                algorithm="$CLEAN_ALGORITHM",
                fullrest_image='${FULLREST_IMAGE}',
                restoring_options='${RESTORING_OPTIONS}',
                channelize=None,lsm="$LSM",**kw0):
    """Makes image(s) from MS. Set dirty and restore to True or False to make the appropriate images. You can also
    set either to a dict of options to be passed to the imager. If restore=True and restore_lsm is True and 'lsm' is set, 
    it will also make a full-restored image (i.e. will restore the LSM into the image) with tigger-restore. Use this when 
    deconvolving residual images. Note that RESTORING_OPTIONS are passed to tigger-restore.
  
    'channelize', if set, overrides the IMAGE_CHANNELIZE setting. If both are None, the options in the 'imager' module take effect.
  
    'algorithm' is the deconvolution algorithm to use (hogbom, clark, csclean, multiscale, entropy) 
  
    'dirty_image', etc. sets the image names, with defaults determined by the globals DIRTY_IMAGE, etc.
    """

