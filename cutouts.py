import numpy as np
from astropy.io import fits as pyfits
from astropy import wcs
import os, sys

def getflux(fname):
   hdulist = pyfits.open(fname)
   table = hdulist[1].data
   flux = table['FLUX']
   hdulist.close()
   return flux

def cutout(kic, q, sfilepath, cluster, ra, dec, cutoutdims):

   # for opening files
   suffixes = {0:2009131105131, 1:2009166043257, 2:2009259160929, 3:2009350155506, 4:2010078095331, 5:2010174085026, 6:2010265121752, 7:2010355172524, 8:2011073133259, 9:2011177032512, 10:2011271113734, 11:2012004120508, 12:2012088054726, 13:2012179063303, 14:2012277125453, 15:2013011073258, 16:2013098041711, 17:2013131215648}
   
   workingdir = os.getcwd()
   os.chdir(sfilepath) # directory containing superstamps

   # get cutout
   serialnumber = 0
   target_x = 0
   target_y = 0
   skip = False
   if cluster == 6791:
      searchrange = np.arange(25,45)
   elif cluster == 6819:
      searchrange = np.arange(45,65)
   for j in searchrange:
      tempfname = f'kplr1000009{j}-{suffixes[q]}_lpd-targ.fits'
      hdulist_temp = pyfits.open(tempfname)
      hd1 = hdulist_temp[1].header
      hd2 = hdulist_temp[2].header
      wtemp = wcs.WCS(hd1, keysel=['binary'])
      x_max = hd2['NAXIS1']-0.5
      y_max = hd2['NAXIS2']-0.5
      target_x, target_y = wtemp.wcs_world2pix(ra, dec, 0)
      if target_x > -0.5 and target_x <= x_max and target_y > -0.5 and target_y <= y_max:
         serialnumber = j
         hdulist_temp.close()
         break
      else:
         continue
   if serialnumber == 0:
      skip = True
      # print(f'KIC {kic} is not available in quarter {q}')

   if skip == False:
      fname = f'kplr1000009{serialnumber}-{suffixes[q]}_lpd-targ.fits'
      hdulist1 = pyfits.open(fname)
      parameters = hdulist1[0].header
      channel = parameters['CHANNEL']
      table = hdulist1[1].data
      flux_full = table['FLUX']
      time1 = table['TIME']
      cadence = table['CADENCENO']
      hd1 = hdulist1[1].header
      w = wcs.WCS(hd1, keysel=['binary'])
      hd2 = hdulist1[2].header
      xdim = hd2['NAXIS1']
      ydim = hd2['NAXIS2']
      hdulist1.close()

      if (channel%2) == 0:
         eo = 0
      else:
         eo = 1

      fulldim = cutoutdims*2-1
      flux1 = np.zeros((len(time1),fulldim,fulldim))
      pointmask = np.zeros((fulldim,fulldim))
      x = int(np.round(target_x))
      y = int(np.round(target_y))
      r = 0
      l = 0
      b = 0
      t = 0
      bolster_r = 0
      bolster_l = 0
      bolster_b = 0
      bolster_t = 0

      # check what adjustments need doing and whether they can be done

      if cluster == 6819:
         serialnumber -= 20

      if xdim-0.5-target_x < cutoutdims: # right edge
         r += 1
         if serialnumber >= 25 and serialnumber <= 34: # right edge, no can do
            missing = int(cutoutdims - (xdim - x))
            # print(f'postage stamp deformed, missing {missing} pixels in the x direction')
            r -= 1
            bolster_r += missing
      # elif np.abs(-0.5-target_x) < cutoutdims-1: # left edge
      if target_x - cutoutdims < -0.5:
         l += 1
         if serialnumber >= 35 and serialnumber <= 44: # left edge, no can do
            missing = int(cutoutdims - np.abs(-1 - x))
            # print(f'postage stamp deformed, missing {missing} pixels in the x direction')
            l -= 1
            bolster_l += missing
      if ydim-0.5-target_y < cutoutdims: # bottom edge
         b += 1
         if serialnumber == 34 or serialnumber == 44: # bottom edge, no can do
            missing = int(cutoutdims - (ydim - y))
            # print(f'postage stamp deformed, missing {missing} pixels in the y direction')
            b -= 1
            bolster_b += missing
      if np.abs(-0.5-target_y) < cutoutdims-1: # top edge
         t += 1
         if serialnumber == 35 or serialnumber == 25: # top edge, no can do
            missing = int(cutoutdims - np.abs(-1 - y))
            # print(f'postage stamp deformed, missing {missing} pixels in the y direction')
            t -=1
            bolster_t += missing

      # then, do whatever adjustments can be done
      if cluster == 6819:
         serialnumber += 20

      if r == 1 and b == 1:
         missingx = int(cutoutdims - (xdim - x))
         missingy = int(cutoutdims - (ydim - y))
         # print('postage stamp deformed, missing', missingx, 'pixels in the x direction')
         # print('postage stamp deformed, missing', missingy, 'pixels in the y direction')
         flux_r = getflux(f'kplr1000009{serialnumber-10}-{suffixes[q]}_lpd-targ.fits')
         flux_b = getflux(f'kplr1000009{serialnumber+1}-{suffixes[q]}_lpd-targ.fits')
         flux_c = getflux(f'kplr1000009{serialnumber-9}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, :(cutoutdims*2-1)-missingy, :(cutoutdims*2-1)-missingx] = flux_full[:, ydim-((cutoutdims*2-1)-missingy):ydim, xdim-((cutoutdims*2-1)-missingx):xdim]
         flux1[:, :(cutoutdims*2-1)-missingy, (cutoutdims*2-1)-missingx:] = flux_r[:, ydim-((cutoutdims*2-1)-missingy):ydim, 0:missingx]
         flux1[:, (cutoutdims*2-1)-missingy:, :(cutoutdims*2-1)-missingx] = flux_b[:, 0:missingy, xdim-((cutoutdims*2-1)-missingx):xdim]
         flux1[:, (cutoutdims*2-1)-missingy:, (cutoutdims*2-1)-missingx:] = flux_c[:, 0:missingy, 0:missingx]
      elif l == 1 and b == 1:
         missingx = int(cutoutdims - np.abs(-1 - x))
         missingy = int(cutoutdims - (ydim - y))
         # print('postage stamp deformed, missing', missingx, 'pixels in the x direction')
         # print('postage stamp deformed, missing', missingy, 'pixels in the y direction')
         flux_l = getflux(f'kplr1000009{serialnumber+10}-{suffixes[q]}_lpd-targ.fits')
         flux_b = getflux(f'kplr1000009{serialnumber+1}-{suffixes[q]}_lpd-targ.fits')
         flux_c = getflux(f'kplr1000009{serialnumber+11}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, :(cutoutdims*2-1)-missingy, missingx:] = flux_full[:, ydim-((cutoutdims*2-1)-missingy):ydim, 0:(cutoutdims*2-1)-missingx]
         flux1[:, :(cutoutdims*2-1)-missingy, :missingx] = flux_l[:, ydim-((cutoutdims*2-1)-missingy):ydim, xdim-missingx:xdim]
         flux1[:, (cutoutdims*2-1)-missingy:, :missingx] = flux_b[:, 0:missingy, 0:missingx]
         flux1[:, (cutoutdims*2-1)-missingy:, missingx:] = flux_c[:, 0:missingy, xdim-((cutoutdims*2-1)-missingx):xdim]
      elif r == 1 and t == 1:
         missingx = int(cutoutdims - (xdim - x))
         missingy = int(cutoutdims - np.abs(-1 - y))
         # print('postage stamp deformed, missing', missingx, 'pixels in the x direction')
         # print('postage stamp deformed, missing', missingy, 'pixels in the y direction')
         flux_r = getflux(f'kplr1000009{serialnumber-10}-{suffixes[q]}_lpd-targ.fits')
         flux_t = getflux(f'kplr1000009{serialnumber-1}-{suffixes[q]}_lpd-targ.fits')
         flux_c = getflux(f'kplr1000009{serialnumber-11}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, missingy:, :(cutoutdims*2-1)-missingx] = flux_full[:, ydim-((cutoutdims*2-1)-missingy):ydim, xdim-((cutoutdims*2-1)-missingx):xdim]
         flux1[:, missingy:, (cutoutdims*2-1)-missingx:] = flux_r[:, ydim-((cutoutdims*2-1)-missingy):ydim, 0:missingx]
         flux1[:, :missingy, :(cutoutdims*2-1)-missingx] = flux_t[:, ydim-missingy:ydim, xdim-((cutoutdims*2-1)-missingx):xdim]
         flux1[:, :missingy, (cutoutdims*2-1)-missingx:] = flux_c[:, ydim-missingy:ydim, 0:missingx]
      elif l == 1 and t == 1:
         missingx = int(cutoutdims - np.abs(-1 - x))
         missingy = int(cutoutdims - np.abs(-1 - y))
         # print('postage stamp deformed, missing', missingx, 'pixels in the x direction')
         # print('postage stamp deformed, missing', missingy, 'pixels in the y direction')
         flux_l = getflux(f'kplr1000009{serialnumber+10}-{suffixes[q]}_lpd-targ.fits')
         flux_t = getflux(f'kplr1000009{serialnumber-1}-{suffixes[q]}_lpd-targ.fits')
         flux_c = getflux(f'kplr1000009{serialnumber+9}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, missingy:, missingx:] = flux_full[:, ydim-((cutoutdims*2-1)-missingy):ydim, xdim-((cutoutdims*2-1)-missingx):xdim] # 4 main
         flux1[:, missingy:, :missingx] = flux_l[:, ydim-((cutoutdims*2-1)-missingy):ydim, xdim-missingx:xdim]
         flux1[:, :missingy, missingx:] = flux_t[:, ydim-missingy:ydim, xdim-((cutoutdims*2-1)-missingx):xdim]
         flux1[:, :missingy, :missingx] = flux_c[:, ydim-missingy:ydim, xdim-missingx:xdim]
      elif r == 1: # right edge
         missing = int(cutoutdims - (xdim - x))
         # print('postage stamp deformed, missing', missing, 'pixels in the x direction')
         flux_r = getflux(f'kplr1000009{serialnumber-10}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, bolster_t:(cutoutdims*2-1)-bolster_b, :(cutoutdims*2-1)-missing] = flux_full[:, y-cutoutdims+1+bolster_t:y+cutoutdims-bolster_b, xdim-((cutoutdims*2-1)-missing):xdim]
         flux1[:, bolster_t:(cutoutdims*2-1)-bolster_b, (cutoutdims*2-1)-missing:] = flux_r[:, y-cutoutdims+1+bolster_t:y+cutoutdims-bolster_b, 0:missing]
      elif l == 1: # left edge
         missing = int(cutoutdims - np.abs(-1 - x))
         # print('postage stamp deformed, missing', missing, 'pixels in the x direction')
         flux_l = getflux(f'kplr1000009{serialnumber+10}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, bolster_t:(cutoutdims*2-1)-bolster_b, missing:] = flux_full[:, y-cutoutdims+1+bolster_t:y+cutoutdims-bolster_b, 0:(cutoutdims*2-1)-missing]
         flux1[:, bolster_t:(cutoutdims*2-1)-bolster_b, :missing] = flux_l[:, y-cutoutdims+1+bolster_t:y+cutoutdims-bolster_b, xdim-missing:xdim]
      elif b == 1: # bottom edge
         missing = int(cutoutdims - (ydim - y))
         # print('postage stamp deformed, missing', missing, 'pixels in the y direction')
         flux_b = getflux(f'kplr1000009{serialnumber+1}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, :(cutoutdims*2-1)-missing, bolster_l:(cutoutdims*2-1)-bolster_r] = flux_full[:, ydim-((cutoutdims*2-1)-missing):ydim, x-cutoutdims+1+bolster_l:x+cutoutdims-bolster_r]
         flux1[:, (cutoutdims*2-1)-missing:, bolster_l:(cutoutdims*2-1)-bolster_r] = flux_b[:, 0:missing, x-cutoutdims+1+bolster_l:x+cutoutdims-bolster_r]
      elif t == 1: # top edge
         missing = int(cutoutdims - np.abs(-1 - y))
         # print('postage stamp deformed, missing', missing, 'pixels in the y direction')
         flux_t = getflux(f'kplr1000009{serialnumber-1}-{suffixes[q]}_lpd-targ.fits')
         flux1[:, missing:, bolster_l:(cutoutdims*2-1)-bolster_r] = flux_full[:, 0:(cutoutdims*2-1)-missing, x-cutoutdims+1+bolster_l:x+cutoutdims-bolster_r] # here flux_3 is main
         flux1[:, :missing, bolster_l:(cutoutdims*2-1)-bolster_r] = flux_t[:, ydim-missing:ydim, x-cutoutdims+1+bolster_l:x+cutoutdims-bolster_r]
      else:
         flux1[:, bolster_t:(cutoutdims*2-1)-bolster_b, bolster_l:(cutoutdims*2-1)-bolster_r] = flux_full[:, y-cutoutdims+1+bolster_t:y+cutoutdims-bolster_b, x-cutoutdims+1+bolster_l:x+cutoutdims-bolster_r]

   else:
      flux1 = 0
      time1 = 0
      eo = 0
      w = 0
      cadence = 0

   os.chdir(workingdir)

   return skip, flux1, time1, eo, w, cadence
