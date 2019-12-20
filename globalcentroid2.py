# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages as pdf
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from cycler import cycler
from astropy.io import fits as pyfits
from astropy import wcs
from astropy.utils.exceptions import AstropyWarning
import nancleaner as nc
import centroid as ct
import smoothing, os, sys, warnings, argparse

warnings.simplefilter('ignore', category=AstropyWarning)

parser = argparse.ArgumentParser(description='Perform image subtraction on NGC 6791 and NGC 6819 targets.')
parser.add_argument('-s', '--sfilepath', required=False, default='./superstamps/', type=str, help='File location of superstamps')
parser.add_argument('-n', '--ngc', required=True, choices=[6791, 6819], type=int, help='Which cluster?')
parser.add_argument('-p', '--plots', dest='show', default=False, type=bool, help='Show plots?')

params = parser.parse_args()

cluster = params.ngc
suffixes = {0:2009131105131, 1:2009166043257, 2:2009259160929, 3:2009350155506, 4:2010078095331, 5:2010174085026, 6:2010265121752, 7:2010355172524, 8:2011073133259, 9:2011177032512, 10:2011271113734, 11:2012004120508, 12:2012088054726, 13:2012179063303, 14:2012277125453, 15:2013011073258, 16:2013098041711, 17:2013131215648}
if cluster == 6791:
   qs = np.arange(1,18)
   searchrange = np.arange(25,45)
elif cluster == 6819:
   qs = [1,2,3,4,5,7,8,9,11,12,13,15,16,17]
   searchrange = np.arange(45,65)

kiclist = np.loadtxt(f'{cluster}bright.dat', delimiter=',')
kics = kiclist[:,0]
ras = kiclist[:,1]
decs = kiclist[:,2]
kps = kiclist[:,3]
centroiddims = 1 # i.e. one pixel either side of target x,y

for q in qs:

   workingdir = os.getcwd()

   globalx = []
   globaly = []

   skipcount = 0

   for k in range(len(kiclist)):

      kic = kics[k]
      ra = ras[k]
      dec = decs[k]

      x = 0
      y = 0
      stamp = 0

      os.chdir(params.sfilepath) # directory containing superstamps
      for j in searchrange:
         fname = f'kplr1000009{j}-{suffixes[q]}_lpd-targ.fits'
         hdulist = pyfits.open(fname)
         hd1 = hdulist[1].header
         table = hdulist[1].data
         hd2 = hdulist[2].header
         w = wcs.WCS(hd1, keysel=['binary'])
         x_max = hd2['NAXIS1']
         y_max = hd2['NAXIS2']
         target_x, target_y = w.wcs_world2pix(ra, dec, 0)
         if target_x > 0 and target_x <= x_max and target_y > 0 and target_y <= y_max:
            x = int(np.round(target_x))
            y = int(np.round(target_y))
            flux = table['FLUX']
            time = table['TIME']
            stamp = j
            hdulist.close()
            break
         else:
            hdulist.close()
            continue

      os.chdir(workingdir)

      try:
         cutout = flux[:,int(y-centroiddims):int(y+centroiddims+1), int(x-centroiddims):int(x+centroiddims+1)]

         x_cent, y_cent = ct.centroid_nomask(len(time), cutout)
         xmed = x_cent-np.nanmedian(x_cent)
         ymed = y_cent-np.nanmedian(y_cent)

         if k == 0:
            globalx += list(xmed)
            globaly += list(ymed)
         else:
            globalx = [sum(i) for i in zip(xmed, globalx)]
            globaly = [sum(i) for i in zip(ymed, globaly)]
      except:
         skipcount += 1

   globalx = [i/(len(kiclist)-skipcount) for i in globalx]
   globaly = [i/(len(kiclist)-skipcount) for i in globaly]
   
   os.chdir('./centroids/')
   np.savetxt(f'{cluster}q{q}centroid.dat', np.c_[globalx, globaly], delimiter=',')
   os.chdir(workingdir)
