import numpy as np
from pyrap.tables import table
#tab = table("SKA1REF2-2h60s_1400MHz_50MHchan.MS",readonly=False)
# pyxis k7ms.ms im.make_image[column=CORRECTED_DATA,restore=True,psf=True]


def AddImagingColumn(msname,add_nonzero_modeldata = False):  #None
	"""
	msname::  str
	modeldata:: array
	"""
	
	ms = table('%s' %msname,readonly=False)
	# This function returns the names of all the columns in the MS.
	desc = ms.getcoldesc("DATA")
	
	if (add_nonzero_modeldata == False):		#if (add_nonzero_modeldata == 'NO') or 'no': #if modeldata == None:
		print "\n >> Add Imaging data column with zeros"
		msData = ms.getcol("DATA")
		modeldata = np.zeros_like(msData, dtype=None, order='K', subok=True)
	else:
		print "\n >> Add Imaging data column with DATA values"
		modeldata = ms.getcol("DATA")
		print modeldata
		
	if 'MODEL_DATA' in ms.colnames():
		pass
	else:		
		
		desc["name"]="MODEL_DATA"
		ms.addcols(desc)
		ms.putcol("MODEL_DATA",modeldata)
		
	if 'CORRECTED_DATA' in ms.colnames():
		pass
	else:		
		
		desc["name"]="CORRECTED_DATA"
		ms.addcols(desc)
		ms.putcol("CORRECTED_DATA",modeldata)
	ms.close()
	
	return None
#
msname = 'k7ms_skymodel.ms' #'k7ms.ms'
AddImagingColumn(msname,add_nonzero_modeldata = True)


print " DONE !!!!!!!!!!"
