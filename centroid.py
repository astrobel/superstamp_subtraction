import numpy as np
import matplotlib.pyplot as plt
import sys

def centroid(timespan, fluxarr, mask):

   # inputs: length of flux array, flux array (3 dimensional), aperture mask for centroiding over only target

   x_cent = np.zeros(timespan)
   y_cent = np.zeros(timespan)

   for i in range(timespan):
      xfsum = 0
      yfsum = 0
      fsum = 0
      temp = fluxarr[i,:]
      
      for index, val in np.ndenumerate(temp):
      
         if mask[index] == 3:
            xfsum += index[1] * temp[index]
            yfsum += index[0] * temp[index]
            fsum += temp[index]
         else:
            pass

      x_cent[i] = xfsum / fsum
      y_cent[i] = yfsum / fsum

   return x_cent, y_cent


def centroid_nomask(timespan, fluxarr):

   # inputs: length of flux array, flux array (3 dimensional), aperture mask for centroiding over only target

   x_cent = np.zeros(timespan)
   y_cent = np.zeros(timespan)

   for i in range(timespan):
      xfsum = 0
      yfsum = 0
      fsum = 0
      temp = fluxarr[i,:]
      for index, val in np.ndenumerate(temp):

         xfsum += index[1] * temp[index]
         yfsum += index[0] * temp[index]
         fsum += temp[index]

      x_cent[i] = xfsum / fsum
      y_cent[i] = yfsum / fsum

   return x_cent, y_cent
