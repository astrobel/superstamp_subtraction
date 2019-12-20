import numpy as np
import scipy.interpolate as spi
import regrid as rg
import nancleaner as nc
import centroid as ct
import os, sys

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

def sapper(flux1, time1, cadence, q, factor, cutoutdims, isolation):

   workingdir = os.getcwd()

   # regridding
   fluxnew, timenew = nc.nancleaner3d(flux1, time1)

   # make mask first before cutting down flux
   mask = gaussmask(factor, flux1.shape[1], flux1.shape[2], isolation)

   # cut down both mask and flux, i.e. get rid of zeros 
   flux1, mask = zerocutter(flux1, mask, factor)

   dump, newcadence = nc.nancleaner3d(flux1, cadence)
   timespan = fluxnew.shape[0]
   oldy, oldx = fluxnew.shape[1], fluxnew.shape[2]
   newy, newx = fluxnew.shape[1]*factor, fluxnew.shape[2]*factor
   re_flux = np.zeros((timespan, newy, newx))

   nanmask = np.asarray(np.where(np.isnan(fluxnew) == True))
   rg_nanmask = np.array((nanmask[0], nanmask[1]*factor, nanmask[2]*factor))

   for i in range(timespan):
      re_flux[i] = rg.regrid(fluxnew[i], factor)

   # regridding base cutout for plotting (oy vey)
   avg_original = np.mean(fluxnew, axis=0)
   regridded_base = rg.regrid(avg_original, factor)
   avgflux = np.mean(re_flux, axis=0)

   # mask making
   maskrange = np.nanmax(avgflux) - np.nanmin(avgflux)
   maskmin =  np.nanmin(avgflux) + 0.05 * maskrange # cut for background
   
   # annulus
   mask_bg = np.zeros((newy, newx))
   counter = 0
   for index, val in np.ndenumerate(avgflux):
      if avgflux[index] < maskmin: # and np.isnan(avgshift[index]) == False:
         mask_bg[index] = -1
         counter += 1

   # calculating background
   lclength = len(timenew)
   for j in range(lclength):
      bg = np.nansum(re_flux[j][np.where(mask_bg==-1)])
      bg /= counter
      re_flux[j] -= bg

   avgflux = np.mean(re_flux, axis=0) # new average image after 
   output = np.c_[newcadence, timenew]

   ap_av = np.nansum(avgflux*mask)

   new_flux = np.zeros(lclength)

   for j in range(lclength):
      placeholder = re_flux[j] - avgflux
      aperture = np.nansum(placeholder*mask)
      # annulus = sum(placeholder[np.where(mask==-1)])
      # annulus /= counter
      new_flux[j] = aperture + ap_av #- annulus

   output = np.c_[output, new_flux]

   return output, mask, avgflux, regridded_base
