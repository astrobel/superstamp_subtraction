import numpy as np
import scipy.interpolate as spi
# import scipy.misc as spm
import regrid as rg
import nancleaner as nc
import centroid as ct
import os, sys

# def circle(radius, factor, y, x): # radius is fraction of aperture taken up by circle
#    mask_rg = np.zeros((y*factor, x*factor))
#    for index, val in np.ndenumerate(mask_rg):
#       if factor % 2 == 1:
#          centre_x = int(np.floor((x*factor)/2))
#          centre_y = int(np.floor((y*factor)/2))
#          dist = np.sqrt((centre_y - index[0])**2 + (centre_x - index[1])**2)
#          if dist < (x*factor)/radius and dist < (y*factor)/radius:
#             mask_rg[index] += 1 # centre
#       else:
#          centre_a = int(np.round(newdim/2))
#          centre_b = int(np.round(newdim/2))
#          dist_a = np.sqrt((centre_a - index[0])**2 + (centre_b - index[1])**2)
#          # dist_b = np.sqrt((centre_b - index[0])**2 + (centre_b - index[1])**2)
#          if dist_a < newdim/radius:# and dist_b < newdim/radius:
#             mask_rg[index] += 1 # centre

#    return mask_rg

def gaussmask(factor, ydim, xdim, isolation): # ydim/xdim should be oldx and oldy

   sigmax = factor*isolation
   sigmay = factor*isolation
   x, y = np.meshgrid(np.arange(0,xdim*factor,1), np.arange(0,ydim*factor,1))
   distx = np.sqrt((x-np.round((xdim/2)*factor))**2)
   disty = np.sqrt((y-np.round((ydim/2)*factor))**2)
   mask = np.exp(-(distx**2/(2*sigmax**2) + disty**2/(2*sigmay**2)))

   return mask

def zerocutter(flux, mask, factor): # flux should be time x 5 x 5, mask 5 x 5

   locs = np.argwhere(np.nanmean(flux, axis=0) == 0)

   x_cut = []
   y_cut = []

   for j in range(len(locs)):
      x_cut.append(locs[j][1])
      y_cut.append(locs[j][0])

   x_cut_u, x_counts = np.unique(x_cut, return_counts=True)
   y_cut_u, y_counts = np.unique(y_cut, return_counts=True)

   x_cut_final = x_cut_u[(x_counts == flux.shape[2])]
   y_cut_final = y_cut_u[(y_counts == flux.shape[1])]
   x_cut_mask = []
   y_cut_mask = []
   
   for i in x_cut_final:
      for j in range(factor):
         x_cut_mask.append((i*factor)+j)
   for i in y_cut_final:
      for j in range(factor):
         y_cut_mask.append((i*factor)+j)

   flux = np.delete(flux, x_cut_final, axis=2)
   flux = np.delete(flux, y_cut_final, axis=1)
   mask = np.delete(mask, x_cut_mask, axis=1)
   mask = np.delete(mask, y_cut_mask, axis=0)

   return flux, mask

def subtract(flux1, time1, cadence, q, factor, cfilepath, centroid_type, cutoutdims, cluster, isolation): # , mask

   workingdir = os.getcwd()

   if centroid_type == 'e':
      os.chdir(cfilepath)
      centroids = np.loadtxt(f'{cluster}q{q}centroid.dat', delimiter=',')
      os.chdir(workingdir)
      x_shifts = centroids[:,0]
      y_shifts = centroids[:,1]

   # make mask first before cutting down flux
   mask = gaussmask(factor, flux1.shape[1], flux1.shape[2], isolation)

   # cut down both mask and flux, i.e. get rid of zeros 
   flux1, mask = zerocutter(flux1, mask, factor)
   
   # prepare arrays/get dimensions etc
   if centroid_type == 'e':
      fluxnew, timenew, x_shifts, y_shifts = nc.nancleaner3d_c(flux1, time1, x_shifts, y_shifts)
   elif centroid_type == 't':
      fluxnew, timenew = nc.nancleaner3d(flux1, time1)
   dump, newcadence = nc.nancleaner3d(flux1, cadence)
   timespan = fluxnew.shape[0]

   # regridding
   oldy, oldx = fluxnew.shape[1], fluxnew.shape[2] # reset dims
   newy, newx = fluxnew.shape[1]*factor, fluxnew.shape[2]*factor
   re_flux = np.zeros((timespan, newy, newx))
   shifted = np.zeros((timespan, newy, newx))

   nanmask = np.asarray(np.where(np.isnan(fluxnew) == True))
   rg_nanmask = np.array((nanmask[0], nanmask[1]*factor, nanmask[2]*factor))

   for i in range(timespan):
      interpolant = spi.RectBivariateSpline(np.arange(oldy), np.arange(oldx), np.nan_to_num(fluxnew[i]), kx=1, ky=1)
      re_flux[i] = interpolant(np.linspace(0, oldy-1, newy), np.linspace(0, oldx-1, newx)) #/ (factor*factor)

   # regridding base cutout for plotting (oy vey)
   avg_original = np.nanmean(fluxnew, axis=0)
   regridded_base = rg.regrid(avg_original, factor)

   # centroiding
   originalmask = np.zeros(flux1[0].shape)#+3
   originalmask[cutoutdims-2:cutoutdims+1,cutoutdims-2:cutoutdims+1] += 3 # 3x3 box for checking
   mask_rg = rg.regrid_slow(originalmask,factor)
   if centroid_type == 't':
      x_cent, y_cent = ct.centroid(timespan, re_flux, mask_rg)
      midx, midy = (newx-1)/2, (newy-1)/2
      x_shifts = (x_cent - midx)
      y_shifts = (y_cent - midy)
   elif centroid_type == 'e':
      x_cent = 0
      y_cent = 0

   for i in range(timespan):
      interpolant = spi.RectBivariateSpline(np.arange(newy), np.arange(newx), re_flux[i,:,:], kx=1, ky=1)
      # interpolant = spi.interp2d(np.arange(newy), np.arange(newx), re_flux[i,:,:], kind='cubic')
      shifted[i,:,:] = interpolant(np.linspace(0-y_shifts[i], newy-y_shifts[i], newy), np.linspace(0-x_shifts[i], newx-x_shifts[i], newx))

   avgshift = np.nanmean(shifted, axis=0)

   # mask making
   maskrange = np.nanmax(avgshift) - np.nanmin(avgshift)
   # maskmax = np.nanmax(avgshift) - 0.5 * maskrange
   maskmin =  np.nanmin(avgshift) + 0.05 * maskrange # cut for background

   # defining six different standard masks
   # maskarr = np.zeros((6,newy,newx))
   # maskarr[0][(cutoutdims-1)*factor:cutoutdims*factor,(cutoutdims-1)*factor:cutoutdims*factor] +=1
   # maskarr[1][(cutoutdims-2)*factor:(cutoutdims+1)*factor,(cutoutdims-1)*factor:cutoutdims*factor] += 1
   # maskarr[1][(cutoutdims-1)*factor:cutoutdims*factor,(cutoutdims-2)*factor:(cutoutdims-1)*factor] += 1
   # maskarr[1][(cutoutdims-1)*factor:cutoutdims*factor,cutoutdims*factor:(cutoutdims+1)*factor] += 1
   # maskarr[2] += 1
   # maskarr[3] += circle(9, factor, oldy, oldx)
   # maskarr[4] += circle(5, factor, oldy, oldx)
   # maskarr[5] += circle(3, factor, oldy, oldx)
   # if mask != None:
   #    maskarr = np.delete(maskarr, np.delete(np.arange(6), mask), axis=0)

   # annulus
   mask_bg = np.zeros((newy, newx))
   counter = 0
   for index, val in np.ndenumerate(avgshift):
      if avgshift[index] < maskmin: # and np.isnan(avgshift[index]) == False:
         mask_bg[index] = -1
         counter += 1

   # calculating background
   lclength = len(timenew)
   for j in range(lclength):
      bg = np.nansum(shifted[j][np.where(mask_bg==-1)])
      bg /= counter
      shifted[j] -= bg

   avgflux = np.nanmean(shifted, axis=0) # new average image after 
   output = np.c_[newcadence, timenew]

   new_flux = np.zeros(lclength)

   ap_av = np.nansum(avgflux*mask)
   for j in range(lclength):
      placeholder = shifted[j] - avgflux
      aperture = np.nansum(placeholder*mask)
      # annulus = sum(placeholder[np.where(mask==-1)])
      # annulus /= counter
      new_flux[j] = aperture + ap_av #- annulus

   output = np.c_[output, new_flux]

   # for deprecated version with multiple masks:
   # for i in range(maskarr.shape[0]):
   #    new_flux = np.zeros(lclength)

   #    ap_av = np.nansum(avgflux[np.where(maskarr[i]==1)])
   #    for j in range(lclength):
   #       placeholder = shifted[j] - avgflux
   #       aperture = np.nansum(placeholder[np.where(maskarr[i]==1)])
   #       # annulus = sum(placeholder[np.where(mask==-1)])
   #       # annulus /= counter
   #       new_flux[j] = aperture + ap_av #- annulus

   #    output = np.c_[output, new_flux]
   
   return output, mask, avgflux, regridded_base, x_cent, y_cent # prev maskarr
