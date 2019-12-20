import numpy as np
import matplotlib.pyplot as plt
import outliers as ol
from astropy.convolution import convolve, Box1DKernel, Gaussian1DKernel

def gausssmooth(time, flux, kern):
   """Takes a time series and applies a Gaussian filter,
         removing edge effects.
      time: ndarray, 1D
      flux: ndarray, 1D
      kern: gaussian smoothing
   """

   gauss_flux = convolve(flux, Gaussian1DKernel(kern), boundary='extend')

   flux2 = flux - gauss_flux

   return flux2, gauss_flux # correct, fit

def fitandclip(time, flux, kernel=100, sigma=3, iters=3):

   for i in range(iters):
      
      # fit

      flux_temp, smth_flux = gausssmooth(time, flux, kernel)

      time_temp, flux_temp = ol.clip(time, flux_temp, sigma)
      # return only the points that were counted as not outliers and repeat

      if i == iters-1:
         time = time_temp
         flux = flux_temp
      else:
         flux = [flux[j] for j in range(len(time)) if time[j] in time_temp]
         time = [time[j] for j in range(len(time)) if time[j] in time_temp]
         kernel /= 2

   return np.asarray(time), np.asarray(flux)
