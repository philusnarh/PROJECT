#!/usr/bin/env python
"""Functions for astronomy coversions between x,y pixel locations and 
world coordinates on the sky with various projections."""

# worldpos.c -- WCS Algorithms from Classic AIPS.
# Copyright (C) 1994
# Associated Universities, Inc. Washington DC, USA.
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS should be addressed as follows:
#        Internet email: aipsmail@nrao.edu
#        Postal address: AIPS Group
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
#              -=-=-=-=-=-=-
#
# These two ANSI C functions, worldpos() and xypix(), perform
# forward and reverse WCS computations for 8 types of projective
# geometries ("-SIN", "-TAN", "-ARC", "-NCP", "-GLS", "-MER", "-AIT"
# and "-STG"):
#
#     worldpos() converts from pixel location to RA,Dec
#     xypix()    converts from RA,Dec to pixel location
#
# where "(RA,Dec)" are more generically (long,lat). These functions
# are based on the WCS implementation of Classic AIPS, an
# implementation which has been in production use for more than ten
# years. See the two memos by Eric Greisen
#
#     ftp://fits.cv.nrao.edu/fits/documents/wcs/aips27.ps.Z
#     ftp://fits.cv.nrao.edu/fits/documents/wcs/aips46.ps.Z
#
# for descriptions of the 8 projective geometries and the
# algorithms.  Footnotes in these two documents describe the
# differences between these algorithms and the 1993-94 WCS draft
# proposal (see URL below). In particular, these algorithms support
# ordinary field rotation, but not skew geometries (CD or PC matrix
# cases). Also, the MER and AIT algorithms work correctly only for
# CRVALi=(0,0). Users should note that GLS projections with yref!=0
# will behave differently in this code than in the draft WCS
# proposal.  The NCP projection is now obsolete (it is a special
# case of SIN).  WCS syntax and semantics for various advanced
# features is discussed in the draft WCS proposal by Greisen and
# Calabretta at:
#      ftp://fits.cv.nrao.edu/fits/documents/wcs/wcs.all.ps.Z
#              -=-=-=-
#
# The original version of this code was Emailed to D.Wells on
# Friday, 23 September by Bill Cotton <bcotton@gorilla.cv.nrao.edu>,
# who described it as a "..more or less.. exact translation from the
# AIPSish..". Changes were made by Don Wells <dwells@nrao.edu>
# during the period October 11-13, 1994:
# 1) added GNU license and header comments
# 2) added testpos.c program to perform extensive circularity tests
# 3) changed float-->double to get more than 7 significant figures
# 4) testpos.c circularity test failed on MER and AIT. B.Cotton
#    found that "..there were a couple of lines of code [in] the wrong
#    place as a result of merging several Fortran routines."
# 5) testpos.c found 0h wraparound in xypix() and worldpos().
# 6) E.Greisen recommended removal of various redundant if-statements,
#    and addition of a 360d difference test to MER case of worldpos().
#
# Work done by Nicholas Chapman <nchapman@u.northwestern.edu> on
# 6-7 April 2011:
# I started with worldpos.c found in NEMO (Peter Teuben's nemo), and
# converted it to python.  Along the way I pythonized the code to raise
# exceptions rather than return error codes.  I also changed some stuff like
# 2pi/4 became pi/2.  I also added three new functions that are high-level
# for user calling.  The _worldpos() and _xypix() are now hidden by default,
# though you can still get them through python.  Lastly, I borrowed the 
# cd matrix code from Doug Mink's wcstools and his version of worldpos.c.
# Note, this seems to work in the limited cases I tested it, but it hasn't
# been resolved what should happen for -ait and -mer projections (doug mink's
# code doesn't address this shortcoming).  Both -ait and -mer projections use
# the cosr and sinr variables, but they are only defined for regular
# rotations.  Also, xinc and yinc are normally set to cdelt1 and cdelt2, but
# don't make sense with a cd matrix.  For now, I have error checking for 
# the use of a cd_matrix with the -ait and -mer projections and throw an
# error in such cases.
#
# I tested the -sin and -tan projections for regular rotations and with
# a cd matrix that has some rotation (but no rotation angle defined!), and
# it seemed to work fine.
#
# 18 July 2011 - I added the SFL projection, which is the same as GLS,
# according to Calabretta & Greisen, 2002, A&A, 395, 1077
#
# 11 Dec 2011 - I fixed a bug with sky2xy() and _xypix().  Converting to 
# pixels using the cd matrix was wrong because we need pixels/degree, but
# cd_matrix is degrees/pix.  So, I took the inverse and pass in those values.
#
# 24 April 2012 - Fixed a stupid bug when trying to compute a single ra/dec
# or x/y coordinate.  Also some better error-trapping when reading the wcs
# with getwcs()
#
# 25 April 2012 - Changed public function to a single class with the available
# methods sky2xy() and xy2sky().  This seems to be more logical.
#
# 31 October 2012 - Added code so that sky2xy() will know about sexagesimal
# (given as a string) and convert the values to degrees.
#
import math   # used by _worldpos() and _xypix()
import sys    # used by _checkproj() and wcs class
import pyfits # used by wcs class

cond2r = math.pi/180.
deps   = 1.0e-5
# allowed projections
ctypes = ("-sin","-tan","-arc","-ncp", "-gls", "-mer", "-ait", "-stg","-sfl")
      
def _checkproj(proj):
   """Checks to see if input projection is a valid one"""
   if proj not in ctypes:
      sys.stderr.write("### Warning! Invalid projection %s.  Assuming linear.\n" %proj)

def _worldpos(xpix, ypix, xref, yref, xrefpix, yrefpix, xinc, yinc, proj,
   rot=None, cd=None):
   """Convert x,y pixel value to world coordinates in degrees

      xpix    : x pixel number
      ypix    : y pixel number
      xref    : x reference coordinate value (deg) (CRVAL1)
      yref    : y reference coordinate value (deg) (CRVAL2)
      xrefpix : x reference pixel (CRPIX1)
      yrefpix : y reference pixel (CRPIX2)
      xinc    : x coordinate increment (deg) (CDELT1)
      yinc    : y coordinate increment (deg) (CDELT2)
      proj    : projection type code e.g. "-SIN"
      rot     : rotation (deg)  (from N through E)
      cd      : list/tuple of four values (the cd_matrix)
      
      Note, you MUST define either rot or cd.  If both are given, cd is 
      used by default.
      
      returns the two coordinates, raises ValueError if the angle is too
      large for projection, and ValueError if xinc or yinc is zero."""

   proj = proj.lower()
   _checkproj(proj)
   # check axis increments - bail out if either 0
   if (xinc == 0.0) or (yinc == 0.0):
      raise ValueError("Input xinc or yinc is zero!")

   # Offset from ref pixel
   dx = xpix - xrefpix
   dy = ypix - yrefpix

   if cd is not None:
      if len(cd) != 4:
         raise IndexError("You must give four values for the cd matrix!")
      if proj in ('-ait','-mer'):
         raise ValueError('cd matrix cannot be used with -AIT or -MER projections!')
      temp = dx*cd[0] + dy*cd[1]
      dy   = dx*cd[2] + dy*cd[3]
      dx   = temp
   elif rot is not None:
      # scale by xinc
      dx = dx * xinc
      dy = dy * yinc
      
      # Take out rotation
      cosr = math.cos(rot*cond2r)
      sinr = math.sin(rot*cond2r)
      if (rot != 0.0):
         temp = dx * cosr - dy * sinr
         dy   = dy * cosr + dx * sinr
         dx   = temp
   else: # both are None
      raise ValueError("You must define either rot or cd keywords!")

   # convert to radians
   ra0    = xref * cond2r
   dec0   = yref * cond2r
   l      = dx * cond2r
   m      = dy * cond2r
   sins   = l*l + m*m
   decout = 0.0
   raout  = 0.0
   cos0   = math.cos(dec0)
   sin0   = math.sin(dec0)

   if proj == '-sin':
      if (sins > 1.0):
         raise ValueError("Angle too large for projection!")
      coss = math.sqrt(1.0 - sins)
      dt = sin0 * coss + cos0 * m
      if abs(dt) > 1:
         raise ValueError("Angle too large for projection!")
      dect = math.asin(dt)
      rat = cos0 * coss - sin0 * m
      if ((rat==0.0) and (l==0.0)):
         raise ValueError("Angle too large for projection!")
      rat = math.atan2(l, rat) + ra0;
   elif proj == '-tan':
      if (sins > 1.0):
         raise ValueError("Angle too large for projection!")
      dect = cos0 - m * sin0
      if (dect == 0.0):
         raise ValueError("Angle too large for projection!")
      rat = ra0 + math.atan2(l, dect)
      dect = math.atan(math.cos(rat-ra0) * (m * cos0 + sin0) / dect)
   elif proj == '-arc':
      if (sins >= math.pi**2):
         raise ValueError("Angle too large for projection!")
      sins = math.sqrt(sins)
      coss = math.cos(sins)
      if (sins != 0.0):
         sins = math.sin(sins) / sins
      else:
         sins = 1.0
      dt = m * cos0 * sins + sin0 * coss
      if abs(dt) > 1:
         raise ValueError("Angle too large for projection!")
      dect = math.asin(dt)
      da = coss - dt * sin0
      dt = l * sins * cos0
      if (da == 0.0) and (dt == 0.0):
         raise ValueError("Angle too large for projection!")
      rat = ra0 + math.atan2(dt, da)
   elif proj == '-ncp': # north celestial pole
      dect = cos0 - m * sin0
      if dect == 0.0:
         raise ValueError("Angle too large for projection!")
      rat = ra0 + math.atan2(l, dect)
      dt = math.cos(rat-ra0)
      if dt == 0.0:
         raise ValueError("Angle too large for projection!")
      dect = dect / dt
      if abs(dect) > 1.0:
         raise ValueError("Angle too large for projection!")
      dect = math.acos(dect)
      if dec0 < 0.0:
         dect = -dect
   elif proj == '-gls' or proj == '-sfl': # global sinusoid, Samson-Flamsteed
      dect = dec0 + m
      if abs(dect) > math.pi/2.:
         raise ValueError("Angle too large for projection!")
      coss = math.cos(dect)
      if abs(l) > math.pi*coss:
         raise ValueError("Angle too large for projection!")
      rat = ra0
      if coss > deps:
         rat = rat + l / coss
   elif proj == '-mer': # mercator
      dt = yinc * cosr + xinc * sinr
      if dt == 0.0:
         dt = 1.0
      dy = (yref/2.0 + 45.0) * cond2r
      dx = dy + dt / 2.0 * cond2r
      dy = math.log(math.tan(dy))
      dx = math.log(math.tan(dx))
      geo2 = dt * cond2r / (dx - dy)
      geo3 = geo2 * dy
      geo1 = math.cos(yref*cond2r)
      if geo1 <= 0.0:
         geo1 = 1.0
      rat = l / geo1 + ra0
      if abs(rat - ra0) > 2*math.pi:
         raise ValueError("Angle too large for projection!") # added 10/13/94 DCW/EWG
      dt = 0.0
      if geo2 != 0.0:
         dt = (m + geo3) / geo2
      dt = math.exp(dt)
      dect = 2.0 * math.atan(dt) - math.pi/2.
   elif proj == '-ait': # Aitoff
      dt = yinc*cosr + xinc*sinr
      if dt == 0.0:
         dt = 1.0
      dt = dt * cond2r
      dy = yref * cond2r
      dx = math.sin(dy+dt)/math.sqrt((1.0 + math.cos(dy+dt))/2.0) - \
           math.sin(dy)/math.sqrt((1.0 + math.cos(dy))/2.0)
      if dx == 0.0:
         dx = 1.0
      geo2 = dt / dx
      dt = xinc*cosr - yinc* sinr
      if dt == 0.0:
         dt = 1.0
      dt = dt * cond2r
      dx = 2.0 * math.cos(dy) * math.sin(dt/2.0)
      if dx == 0.0:
         dx = 1.0
      geo1 = dt * math.sqrt((1 + math.cos(dy)*math.cos(dt/2.0))/2.0) / dx
      geo3 = geo2 * math.sin(dy) / math.sqrt((1 + math.cos(dy))/2.0)
      rat  = ra0
      dect = dec0
      if (l == 0.0) and (m==0.0):
         pass
      else:
         dz = 4.0 - l*l/(4.0*geo1*geo1) - ((m+geo3)/geo2)**2
         if (dz > 4.0) or (dz < 2.0):
            raise ValueError("Angle too large for projection!")
         dz = 0.5 * math.sqrt(dz)
         dd = (m + geo3) * dz / geo2
         if abs(dd) > 1.0:
            raise ValueError("Angle too large for projection!")
         dd = math.asin(dd)
         if abs(math.cos(dd)) < deps:
            raise ValueError("Angle too large for projection!")
         da = l * dz / (2.0 * geo1 * math.cos(dd))
         if abs(da) > 1.0:
            raise ValueError("Angle too large for projection!")
         da = math.asin(da)
         rat = ra0 + 2.0 * da
         dect = dd
   elif proj == '-stg': # Sterographic
      dz = (4.0 - sins) / (4.0 + sins)
      if abs(dz) > 1.0:
         raise ValueError("Angle too large for projection!")
      dect = dz * sin0 + m * cos0 * (1.0+dz) / 2.0
      if abs(dect) > 1.0:
         raise ValueError("Angle too large for projection!")
      dect = math.asin(dect)
      rat  = math.cos(dect)
      if abs(rat) < deps:
         raise ValueError("Angle too large for projection!")
      rat = l * (1.0+dz) / (2.0 * rat)
      if abs(rat) > 1.0:
         raise ValueError("Angle too large for projection!")
      rat = math.asin(rat)
      mg = 1.0 + math.sin(dect) * sin0 + math.cos(dect) * cos0 * math.cos(rat)
      if abs(mg) < deps:
         raise ValueError("Angle too large for projection!")
      mg = 2*(math.sin(dect)*cos0 - math.cos(dect)*sin0*math.cos(rat)) / mg
      if abs(mg-m) > deps:
         rat = math.pi - rat
      rat = ra0 + rat
   else: # default is linear
      rat  =  ra0 + l
      dect = dec0 + m

   # return ra in range
   raout = rat;
   decout = dect;
   if raout-ra0 > math.pi:
      raout = raout - 2*math.pi
   elif raout-ra0 < -math.pi:
      raout = raout + 2*math.pi
   if raout < 0.0:
      raout += 2*math.pi # added by DCW 10/12/94

   # correct units back to degrees
   xpos  = raout  / cond2r
   ypos  = decout  / cond2r
   return xpos,ypos

def _xypix(xpos, ypos, xref, yref, xrefpix, yrefpix, xinc, yinc, proj, 
   rot=None, dc=None):
   """Convert input ra,dec to x,y pixels

      xpos    : x (RA) coordinate (deg)
      ypos    : y (dec) coordinate (deg)
      xref    : x reference coordinate value (deg) (CRVAL1)
      yref    : y reference coordinate value (deg) (CRVAL2)
      xrefpix : x reference pixel (CRPIX1)
      yrefpix : y reference pixel (CRPIX2)
      xinc    : x coordinate increment (deg) (CDELT1)
      yinc    : y coordinate increment (deg) (CDELT2)
      proj    : projection type code e.g. -SIN (CTYPE1/2)
      rot     : rotation (deg)  (from N through E)
      dc      : list/tuple of four values that are the inverse of the
                cd_matrix in the FITS header
      
      Note, you MUST define either rot or dc.  If both are given, dc is 
      used by default.
      
      returns the x,y pixel positions, or raises ValueError for angles too
      large for projection, ValueError if xinc or yinc is zero, and
      Arithmetic error in one instance. """

   proj = proj.lower()
   _checkproj(proj)
   # check axis increments - bail out if either 0
   if (xinc == 0.0) or (yinc == 0.0):
      raise ValueError("Input xinc or yinc is zero!")

   # 0h wrap-around tests added by D.Wells 10/12/94:
   dt = (xpos - xref)
   if (dt > +180):
      xpos -= 360
   if (dt < -180):
      xpos += 360

   if dc is not None and proj in ('-ait','-mer'):
      raise ValueError('cd matrix cannot be used with -AIT or -MER projections!')
   elif rot is not None:
      cosr = math.cos(rot * cond2r)
      sinr = math.sin(rot * cond2r)

   # Non linear position
   ra0  = xref * cond2r
   dec0 = yref * cond2r
   ra   = xpos * cond2r
   dec  = ypos * cond2r

   # compute direction cosine
   coss = math.cos(dec)
   sins = math.sin(dec)
   l    = math.sin(ra-ra0) * coss
   sint = sins * math.sin(dec0) + coss * math.cos(dec0) * math.cos(ra-ra0)

   if proj == '-sin':
      if sint < 0.0:
         raise ValueError("Angle too large for projection!")
      m = sins * math.cos(dec0) - coss * math.sin(dec0) * math.cos(ra-ra0)
   elif proj == '-tan':
      if sint <= 0.0:
         raise ValueError("Angle too large for projection!")
      m = sins * math.sin(dec0) + coss * math.cos(dec0) * math.cos(ra-ra0)
      l = l / m
      m = (sins*math.cos(dec0) - coss*math.sin(dec0) * math.cos(ra-ra0)) / m
   elif proj == '-arc':
      m = sins * math.sin(dec0) + coss * math.cos(dec0) * math.cos(ra-ra0)
      if m < -1.0:
         m = -1.0
      elif m > 1.0:
         m = 1.0
      m = math.acos(m)
      if m != 0:
         m = m / math.sin(m)
      else:
         m = 1.0
      l = l * m
      m = m*(sins*math.cos(dec0) - coss*math.sin(dec0) * math.cos(ra-ra0))
   elif proj == '-ncp': # North celestial pole
      if dec0 == 0.0:
         raise ValueError("Angle too large for projection!") # can't stand the equator
      else:
         m = (math.cos(dec0) - coss * math.cos(ra-ra0)) / math.sin(dec0)
   elif proj == '-gls' or proj == '-sfl': # global sinusoid, samson-flamsteed
      dt = ra - ra0
      if abs(dec) > math.pi/2.:
         raise ValueError("Angle too large for projection!")
      if abs(dec0) > math.pi/2.:
         raise ValueError("Angle too large for projection!")
      m = dec - dec0
      l = dt * coss
   elif proj == '-mer': # mercator
      dt = yinc * cosr + xinc * sinr
      if dt == 0.0:
         dt = 1.0
      dy = (yref/2.0 + 45.0) * cond2r
      dx = dy + dt / 2.0 * cond2r
      dy = math.log(math.tan(dy))
      dx = math.log(math.tan(dx))
      geo2 = dt * cond2r / (dx - dy)
      geo3 = geo2 * dy
      geo1 = math.cos(yref*cond2r)
      if geo1 <= 0.0:
         geo1 = 1.0
      dt = ra - ra0
      l  = geo1 * dt
      dt = dec / 2.0 + math.pi/4.
      dt = math.tan(dt)
      if dt < deps:
         raise ArithmeticError("dt < %f" %deps)
      m = geo2 * math.log(dt) - geo3
   elif proj == '-ait': # Aitoff
      l = 0.0
      m = 0.0
      da = (ra - ra0) / 2.0
      if abs(da) > math.pi/2.:
         raise ValueError("Angle too large for projection!")
      dt = yinc*cosr + xinc*sinr
      if dt == 0.0:
         dt = 1.0
      dt = dt * cond2r
      dy = yref * cond2r
      dx = math.sin(dy+dt)/math.sqrt((1+math.cos(dy+dt))/2.) - \
           math.sin(dy)/math.sqrt((1+math.cos(dy))/2.)
      if dx == 0.0:
         dx = 1.0
      geo2 = dt / dx
      dt = xinc*cosr - yinc* sinr
      if dt == 0.0:
         dt = 1.0
      dt = dt * cond2r
      dx = 2.0 * math.cos(dy) * math.sin(dt/2.0)
      if dx == 0.0:
         dx = 1.0
      geo1 = dt * math.sqrt((1.0+math.cos(dy)*math.cos(dt/2.0))/2.) / dx
      geo3 = geo2 * math.sin(dy) / math.sqrt((1+math.cos(dy))/2.)
      dt = math.sqrt ((1 + math.cos(dec) * math.cos(da))/2.)
      if abs(dt) < deps:
         raise ZeroDivisionError("dt < %f" %deps)
      l = 2.0 * geo1 * math.cos(dec) * math.sin(da) / dt
      m = geo2 * math.sin(dec) / dt - geo3
   elif proj == '-stg': # Sterographic
      da = ra - ra0
      if abs(dec) > math.pi/2.:
         raise ValueError("Angle too large for projection!")
      dd = 1.0 + sins * math.sin(dec0) + coss * math.cos(dec0) * math.cos(da)
      if abs(dd) < deps:
         raise ValueError("Angle too large for projection!")
      dd = 2.0 / dd
      l = l * dd
      m = dd * (sins * math.cos(dec0) - coss * math.sin(dec0) * math.cos(da))
   else: # linear
      l = cond2r*(xpos - xref)
      m = cond2r*(ypos - yref)

   # back to degrees 
   dx = l / cond2r
   dy = m / cond2r

   if dc is not None:
      if len(dc) != 4:
         raise IndexError("You must give four values for the cd matrix!")
      dz = dx*dc[0] + dy*dc[1]
      dy = dx*dc[2] + dy*dc[3]
      dx = dz
   elif rot is not None:
      # correct for rotation
      dz = dx*cosr + dy*sinr
      dy = dy*cosr - dx*sinr
      dx = dz

      # correct for xinc,yinc
      dx = dx / xinc
      dy = dy / yinc
   else: # both are None
      raise ValueError("You must define either rot or cd keywords!")

   # convert to pixels
   xpix = dx + xrefpix
   ypix = dy + yrefpix
   return xpix,ypix

class wcs:
   def __init__(self,filename,ext=0,rot=0,cd=False):
      """Get the WCS geometry
      
         filename - name of the fits file
         ext      - extension number to read for header info (start from zero)
         rot      - name for rotation angle header keyword.  Can also be an angle
                    in degrees
         cd       - set to true to use CD matrix instead of rot."""
      self.header = {'rot' : None, 'cd' : None, 'dc' : None}
      if isinstance(filename,str):
         self.filename = filename
         img = pyfits.open(filename)
      else:
         img = filename
         self.filename = None
         # assume a pyfits.open hdu was passed as sys.stdin
         # TODO: Doesn't work currently
         #self.filename = None
         #img = pyfits.open(filename,mode='readonly')
      self.ext = ext
      hdrkeys = img[ext].header.keys() # list of existing keywords ascardlist().
      if 'CRVAL1' not in hdrkeys:
         raise KeyError("CRVAL1 missing from header!")
      if 'CRVAL2' not in hdrkeys:
         raise KeyError("CRVAL2 missing from header!")
      if 'CRPIX1' not in hdrkeys:
         raise KeyError("CRPIX1 missing from header!")
      if 'CRPIX2' not in hdrkeys:
         raise KeyError("CRPIX2 missing from header!")
      if 'NAXIS1' not in hdrkeys:
         raise KeyError("NAXIS1 missing from header!")
      if 'NAXIS2' not in hdrkeys:
         raise KeyError("NAXIS2 missing from header!")
      if 'CTYPE1' not in hdrkeys:
         raise KeyError("CTYPE1 missing from header!")
      self.header['crval1'] = float(img[ext].header['crval1'])
      self.header['crval2'] = float(img[ext].header['crval2'])
      self.header['crpix1'] = float(img[ext].header['crpix1'])
      self.header['crpix2'] = float(img[ext].header['crpix2'])
      self.header['naxis1'] = float(img[ext].header['naxis1'])
      self.header['naxis2'] = float(img[ext].header['naxis2'])
      self.header['proj']   = img[ext].header['ctype1'][-4:] # just want projection

      if cd is True: # use CD matrix
         #wcs['cdelt1'] = float(img[ext].header['cd1_1'])
         #wcs['cdelt2'] = float(img[ext].header['cd2_2'])
         cd = [0.,0.,0.,0.]
         dc = [0.,0.,0.,0.] # inverse of cd
         missing = 0
         if 'CD1_1' in hdrkeys:
            cd[0] = float(img[ext].header['cd1_1'])
         else:
            sys.stderr.write("### Warning! Cannot CD1_1 keyword\n")
            missing += 1
         if 'CD1_2' in hdrkeys:
            cd[1] = float(img[ext].header['cd1_2'])
         else:
            sys.stderr.write("### Warning! Cannot CD1_2 keyword\n")
            missing += 1
         if 'CD2_1' in hdrkeys:
            cd[2] = float(img[ext].header['cd2_1'])
         else:
            sys.stderr.write("### Warning! Cannot CD2_1 keyword\n")
            missing += 1
         if 'CD2_2' in hdrkeys:
            cd[3] = float(img[ext].header['cd2_2'])
         else:
            sys.stderr.write("### Warning! Cannot CD2_2 keyword\n")
            missing += 1
         if missing == 4:
            raise KeyError("Cannot find CD Matrix!")

         self.header['cd'] = cd
         for i in range(len(cd)):
            try:
               dc[i] = 1./cd[i]
            except ZeroDivisionError:
               dc[i] = 0.0
         self.header['dc'] = dc
      else: # use rotation
         try:
            self.header['cdelt1'] = float(img[ext].header['cdelt1'])
         except KeyError:
            raise KeyError("CDELT1 missing from header!")
         try:
            self.header['cdelt2'] = float(img[ext].header['cdelt2'])
         except KeyError:
            raise KeyError("CDELT2 missing from header!")

         if isinstance(rot,str): # is a keyword
            if rot in hdrkeys:
               self.header['rot'] = float(img[ext].header[rot])
            else:
               raise KeyError("Cannot find rotation angle keyword %s!" %rot)
         else:
            self.header['rot'] = rot
      if self.filename is None:
         pass
      else:
         img.close()

   def sky2xy(self,ra,dec):
      """Convert ra,dec into x,y pixels

         ra,dec - can be either single values or iterable list/tuples/etc"""

      wcs = self.header
      # if ra or dec are strings (presumably sexagesimal, convert to degrees)
      if isinstance(ra,str):
         ra = self._sex2deg(ra,'ra')
      if isinstance(dec,str):
         dec = self._sex2deg(dec,'dec')
      try: # if ra,dec are iterable, then iterate over all values
         n1 = len(ra)
         n2 = len(dec)
         if n1 != n2:
            raise IndexError("number of values in ra and dec are not equal!")
         x = []
         y = []
         # ensure all values in iterable are in degrees
         ra  = map(lambda a: self._sex2deg(a,'ra'),ra)
         dec = map(lambda a: self._sex2deg(a,'dec'),dec)
         if wcs['cd'] is not None:
            for a,b in zip(ra,dec):
               t1,t2 = _xypix(a, b, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
                  wcs['crpix2'], wcs['cd'][0], wcs['cd'][3], wcs['proj'], 
                  dc=wcs['dc'])
               x.append(t1)
               y.append(t2)
         elif wcs['rot'] is not None:
            for a,b in zip(ra,dec):
               t1,t2 = _xypix(a, b, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
                  wcs['crpix2'], wcs['cdelt1'], wcs['cdelt2'], wcs['proj'],
                  rot=wcs['rot'])
               x.append(t1)
               y.append(t2)
         else:
            raise KeyError("Either rot or cd must be specified with wcs!")
      except TypeError: # ra,dec are single values
         if wcs['cd'] is not None:
            x,y = _xypix(ra, dec, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
               wcs['crpix2'], wcs['cd'][0], wcs['cd'][3], wcs['proj'], 
               dc=wcs['dc'])
         elif wcs['rot'] is not None:
            x,y = _xypix(ra, dec, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
               wcs['crpix2'], wcs['cdelt1'], wcs['cdelt2'], wcs['proj'],
               rot=wcs['rot'])
         else:
            raise KeyError("Either rot or cd must be specified with wcs!")      
      return x,y

   def xy2sky(self,x,y):
      """Convert x,y into ra,dec

         x,y - can be either single values or iterable list/tuples/etc"""

      wcs = self.header
      try: # if x,y are iterable, then iterate over all values
         n1 = len(x)
         n2 = len(y)
         if n1 != n2:
            raise IndexError("number of values in x and y are not equal!")
         ra = []
         dec = []
         for a,b in zip(x,y):
            if wcs['cd'] is not None:
               t1,t2 = _worldpos(a, b, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
                  wcs['crpix2'], wcs['cd'][0], wcs['cd'][3], wcs['proj'], 
                  cd=wcs['cd'])
            elif wcs['rot'] is not None:
               t1,t2 = _worldpos(a, b, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
                  wcs['crpix2'], wcs['cdelt1'], wcs['cdelt2'], wcs['proj'],
                  rot=wcs['rot'])
            else:
               raise KeyError("Either rot or cd must be specified with wcs!")
            ra.append(t1)
            dec.append(t2)
      except TypeError: # x,y are single values
         if wcs['cd'] is not None:
            ra,dec = _worldpos(x, y, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
               wcs['crpix2'], wcs['cd'][0], wcs['cd'][3], wcs['proj'], 
               cd=wcs['cd'])
         elif wcs['rot'] is not None:
            ra,dec = _worldpos(x, y, wcs['crval1'], wcs['crval2'], wcs['crpix1'], 
               wcs['crpix2'], wcs['cdelt1'], wcs['cdelt2'], wcs['proj'],
               rot=wcs['rot'])
         else:
            raise KeyError("Either rot or cd must be specified with wcs!")      
      return ra,dec

   def _sex2deg(self,value,coord):
      '''Convert sexagesimal to degrees, unless it is already in degrees'''

      if isinstance(value,str): # is a string
         factor = 1.0    # divide each part by a factor to convert it to degrees
         negFlag = False # set to true if negative declination
         tmp = map(float,value.split(':'))
         if tmp[0] < 0:
            negFlag = True
         for i in range(len(tmp)): # convert each part to degrees
            tmp[i] = abs(tmp[i])/factor
            factor = factor*60.0

         degrees = sum(tmp)
         if coord == 'ra':
            if len(tmp) == 1: # if only 1 value input, assume it was already in
               return degrees # degrees
            else:
               return 15*degrees # convert hours to degrees
         elif coord == 'dec':
            if negFlag:
               return -1*degrees
            else:
               return degrees
      elif isinstance(value,(float,int)): # assume it is already in degrees
         return value

