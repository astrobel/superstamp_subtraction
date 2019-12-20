import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gs
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
from astropy.stats import LombScargle
from astropy.utils.exceptions import AstropyWarning
import cutouts as cut
import sap_photometry as sap
import eventfixer as fx
import fc2 as fc
import os, sys, warnings, argparse

C = ['#e6e7e8','#e5e6e7','#e4e6e7','#e4e5e6','#e3e5e5','#e3e4e5','#e2e4e4','#e2e3e4','#e1e2e3','#e1e2e3','#e0e1e2','#e0e1e1','#dfe0e1','#dfdfe0','#dedfe0','#dededf','#dddedf','#dcddde','#dcdddd','#dbdcdd','#dbdbdc','#dadbdc','#dadadb','#d9dada','#d9d9da','#d8d9d9','#d8d8d9','#d7d7d8','#d7d7d8','#d6d6d7','#d6d6d6','#d5d5d6','#d4d4d5','#d4d4d5','#d3d3d4','#d3d3d3','#d2d2d3','#d2d2d2','#d1d1d2','#d1d0d1','#d0d0d1','#d0cfd0','#cfcfcf','#cfcecf','#cecdce','#cecdce','#cdcccd','#cccccc','#cccbcc','#cbcbcb','#cbcacb','#cac9ca','#cac9ca','#c9c8c9','#c9c8c8','#c8c7c8','#c8c6c7','#c7c6c7','#c7c5c6','#c6c5c5','#c6c4c5','#c5c4c4','#c4c3c4','#c4c2c3','#c3c2c3','#c3c1c2','#c2c0c1','#c1c0c1','#c1bfc0','#c0bebf','#bfbebe','#bfbdbe','#bebcbd','#bdbcbc','#bdbbbc','#bcbabb','#bbbaba','#bbb9ba','#bab8b9','#b9b8b8','#b8b7b8','#b8b6b7','#b7b6b6','#b6b5b6','#b6b4b5','#b5b4b4','#b4b3b4','#b4b2b3','#b3b1b2','#b2b1b2','#b2b0b1','#b1afb0','#b0afb0','#b0aeaf','#afadae','#aeadae','#aeacad','#adabac','#acabab','#acaaab','#aba9aa','#aaa9a9','#aaa8a9','#a9a7a8','#a8a7a7','#a8a6a7','#a7a5a6','#a6a5a5','#a5a4a5','#a5a3a4','#a4a3a3','#a3a2a3','#a3a1a2','#a2a1a1','#a1a0a1','#a19fa0','#a09e9f','#9f9e9f','#9f9d9e','#9e9c9d','#9d9c9d','#9d9b9c','#9c9a9b','#9b9a9b','#9b999a','#9a9899','#999898','#999798','#989697','#979697','#979596','#969596','#969495','#959494','#959394','#949393','#939293','#939292','#929192','#929091','#919090','#918f90','#908f8f','#908e8f','#8f8e8e','#8e8d8e','#8e8d8d','#8d8c8d','#8d8c8c','#8c8b8b','#8c8a8b','#8b8a8a','#8b898a','#8a8989','#898889','#898888','#888788','#888787','#878686','#878586','#868585','#868485','#858484','#848384','#848383','#838283','#838282','#828181','#828181','#818080','#818080','#807f7f','#7f7e7f','#7f7e7e','#7e7d7e','#7e7d7d','#7d7c7c','#7d7c7c','#7c7b7b','#7c7b7b','#7b7a7a','#7a7a7a','#7a7979','#797879','#797878','#787777','#787777','#777676','#767676','#767575','#757575','#757474','#747474','#747373','#747373','#737373','#737272','#737272','#727272','#727171','#727171','#717171','#717070','#717070','#707070','#706f6f','#6f6f6f','#6f6f6f','#6f6e6e','#6e6e6e','#6e6e6e','#6e6d6d','#6d6d6d','#6d6d6d','#6d6c6c','#6c6c6c','#6c6c6c','#6c6b6b','#6b6b6b','#6b6a6a','#6b6a6a','#6a6a6a','#6a6969','#6a6969','#696969','#696868','#686868','#686868','#686767','#676767','#676767','#676666','#666666','#666666','#666565','#656565','#656565','#656464','#646464','#646464','#646363','#636363','#636363','#626262','#626262','#626262','#616161','#616161','#616161','#606060','#606060','#606060','#5f5f5f','#5f5f5f','#5f5f5f','#5e5e5e']
middling_grey = mpl.colors.ListedColormap(C)
quarterlines = [120, 131, 169, 259, 351, 442, 538, 629, 734, 808, 906, 1000, 1098, 1182, 1273, 1372, 1471, 1558]

warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=UserWarning)

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
      r'\usepackage{helvet}',
      r'\usepackage[EULERGREEK]{sansmath}',
      r'\sansmath'
]
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True

class MidpointNormalize(colors.Normalize):
   # from https://matplotlib.org/users/colormapnorms.html#custom-normalization-two-linear-ranges
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

makefigs = True # 0: no plots, 1: yes plots

parser = argparse.ArgumentParser(description='Perform image subtraction on NGC 6791 and NGC 6819 targets.')
parser.add_argument('-k', '--kic', required=True, type=int, help='KIC ID')
parser.add_argument('-s', '--stampfilepath', required=False, default='./superstamps/', type=str, help='File location of superstamps')
parser.add_argument('-r', '--regridfactor', required=False, default=2, type=int, help='Regridding factor')
parser.add_argument('-m', '--missing', required=False, default=None, nargs='+', choices=np.arange(1,18), type=int, help='Any quarters to cut?')

params = parser.parse_args()

kic = params.kic
factor = params.regridfactor
missing = params.missing

cluster = 0
kic = str(kic)
if kic.startswith('2'):
   cluster = 6791
   quarterlist = list(np.arange(1, 18))
elif kic.startswith('5') or kic.startswith('4'):
   cluster = 6819
   quarterlist = [1,2,3,4,5,7,8,9,11,12,13,15,16,17]
if missing != None:
   quarterlist = [i for i in quarterlist if i not in missing]
kic = int(kic)
emptyquarters = []

members = pd.read_csv(f'kiclist{cluster}.csv')
b = members['phot_bp_mean_mag']
r = members['phot_rp_mean_mag']
g = members['phot_g_mean_mag']
br = members['bp_rp']
meanprob = members['meanprob']
kics = members['kic_kepler_id']
ras = members['ra_2000']
decs = members['dec_2000']
teff = members['teff_val']
Kp = members['kic_kepmag']
crowding = members['crowding_mag']

target_ra = ras[kics ==  kic].values[0]
target_dec = decs[kics == kic].values[0]
cutoutdims = 3
residuals = crowding/Kp
isolation = residuals[kics == kic].values[0]

try:
   os.mkdir(str(kic))
except FileExistsError:
   pass

kern = 100 # for smoothing
sigma = (factor**2) * 2 # for clipping

time = []
x_cent = []
y_cent = []

raw = np.zeros((1,2))
masked = np.zeros((1,2))
means = []

for i, q in enumerate(quarterlist):

   skip, flux1, time1, eo, w, cadence = cut.cutout(kic, q, params.stampfilepath, cluster, target_ra, target_dec, cutoutdims)
   if skip == False:
      output, maskarr, avgflux, regridded_base = sap.sapper(flux1, time1, cadence, q, factor, cutoutdims, isolation)
      tosave = pd.DataFrame(data=output)
      fileheader = ['cadence', 'time (BJD-2454833)', 'flux']
      os.chdir(f'./{kic}/')
      tosave.to_csv(f'kic{kic}_q{q}_{factor}_sap.dat', sep=' ', float_format='%.5f', header=fileheader, index=False)
      os.chdir('..')

      # for column in range(nummasks):
      flux = output[:,2]
      time = output[:,1]
      means.append(np.nanmean(flux))
      time, flux = fx.fixer(q, time, flux, cluster)
      temp = np.c_[time, flux]
      raw = np.r_[raw, temp]
      time, flux = fc.fitandclip(time, flux)
      temp = np.c_[time, flux]
      masked = np.r_[masked, temp]

      avgflux = np.flipud(avgflux)
      if eo == 0:
         avgflux = np.fliplr(avgflux)
   elif skip == True:
      emptyquarters.append(q)

for q in emptyquarters:
   quarterlist.remove(q)

raw = np.delete(raw, 0, axis=0)
masked = np.delete(masked, 0, axis=0)
masked[:,1] = masked[:,1]/np.nanmean(means)

### DATA HANDLING ###

# for amplitude summary plot
medbins = int(280 / 30)
bins = np.zeros(medbins)
maxforplotting = 0

newdim = regridded_base.shape[0]

# for i in range(nummasks):

# fitted flux
time = masked[:,0]
flux = masked[:,1]

# amplitude spectrum
ofac = int(5)
hifac = int(1)
freqs, ampls = LombScargle(np.asarray(time), np.asarray(flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
hifac *= (283/11.57)/max(freqs)
freqs, ampls = LombScargle(np.asarray(time), np.asarray(flux)).autopower(method='fast', normalization='psd', samples_per_peak=ofac, nyquist_factor=hifac)
ampls = ampls * 4. / len(time)
ampls = np.sqrt(ampls)
ampls *= 1e6
freqs *= 11.57
# exec("freqs%d = freqs" % i)
# exec("ampls%d = ampls" % i)
maxforplotting = np.max(ampls)

# for amplitudes summary plot
for j in range(medbins):
   bins[j] = np.median(ampls[(freqs > j*30) & (freqs <= (j+1)*30)])

os.chdir(f'./{kic}/')
stardata = open(f'{kic}_sap.dat', 'w')
stardata.write(f'Teff = {teff[kics==kic].values[0]:.0f} K\n')
stardata.write(f'G = {g[kics==kic].values[0]:.2f}\n')
stardata.write(f'B - R = {br[kics==kic].values[0]:.2f}\n')
stardata.write(f'membership: mean gaussian posterior prob = {meanprob[kics==kic].values[0]:.2f}\n')
stardata.write(f'high freq noise = {bins[-1]:.2f} ppm\n')
stardata.write(f'S/N = {max(ampls)/bins[-1]:.2f}\n')
stardata.close()
os.chdir('..')

cmap = pl.cm.Greys_r
maskcmap = cmap(np.arange(cmap.N))
maskcmap[:,-1] = np.linspace(1, 0, cmap.N)
maskcmap = ListedColormap(maskcmap)

if makefigs == True:

   plt.figure(figsize=(15,10))
   grid = gs.GridSpec(8,5, width_ratios=[1,2,1.5,1,1])

   img = plt.subplot(grid[0:3,0])
   amp = plt.subplot(grid[0:4,1])
   cmd = plt.subplot(grid[4:,0:2])
   first = plt.subplot(grid[0:2,2:])
   smth = plt.subplot(grid[2:4,2:])
   zoom = plt.subplot(grid[4:6,2:])
   ft = plt.subplot(grid[6:,2:])
   
   # exec("freqs = freqs%d" % i)
   # exec("ampls = ampls%d" % i)

   # star and mask
   img.set_title(f'{kic}')
   star = img.imshow(regridded_base, cmap='OrRd')
   img.imshow(maskarr, cmap=maskcmap, alpha=0.5) # for one gaussian mask only
   star.axes.get_xaxis().set_ticklabels([])
   star.axes.get_yaxis().set_ticklabels([])
   star.axes.get_xaxis().set_ticks([])
   star.axes.get_yaxis().set_ticks([])
   if np.isnan(teff[kics==kic].values[0]) == False:
      img.text(0.5, newdim+factor, r'T$_{eff}$'+f' = {teff[kics==kic].values[0]:.0f} K')
   img.text(0.5, newdim+(factor*1.5), f'G = {g[kics==kic].values[0]:.2f}')
   img.text(0.5, newdim+(factor*2), f'B - R = {br[kics==kic].values[0]:.2f}')
   img.text(0.5, newdim+(factor*2.5), 'membership: gaussian')
   img.text(0.5+(0.5*factor), newdim+(factor*3), f'posterior prob = {meanprob[kics==kic].values[0]:.2f}')
   img.text(0.5, newdim+(factor*3.5), f'high freq noise = {bins[-1]:.2f} ppm')
   img.text(0.5, newdim+(factor*4), f'S/N = {max(ampls)/bins[-1]:.2f}')

   # CMD
   cmdplot = cmd.scatter(br, g, c=meanprob, cmap=middling_grey, s=3)
   cmd.plot(br[kics==kic].values[0], g[kics==kic].values[0], 'rX', fillstyle='none', ms=9)
   if br[kics==kic].values[0] > 0.5 and br[kics==kic].values[0] < 2.2:
      cmd.set_xlim(0.5, 2.2)
      cmd.set_ylim(ymin=12.5)
   elif br[kics==kic].values[0] < 0.5:
      cmd.set_xlim(xmax=2.2)
      cmd.set_ylim(ymin=12.5)
   elif br[kics==kic].values[0] > 2.2:
      cmd.set_xlim(xmin=0.5)
   else:
      pass
   cmd.invert_yaxis()
   cmd.set_xlabel('B - R')
   cmd.set_ylabel('G')
   cbar = plt.colorbar(cmdplot, ax=cmd)
   cbar.set_label('posterior prob')

   # amplitude summary
   amp.xaxis.set_major_locator(MaxNLocator(integer=True))
   main = amp.scatter(range(medbins), bins, c='r', linewidths=0.5, edgecolors='k')
   amp.set_xlabel('bin number')
   amp.set_ylabel('median amplitude by bin (ppm)')

   # raw flux
   time = raw[:,0]
   flux = raw[:,1]
   for q in range(len(quarterlines)):
      first.axvline(x=quarterlines[q], color='red', linestyle='--', alpha=0.7)
   first.plot(time, flux, 'k.', ms=3)
   first.set_xlim(quarterlines[0], 1600)
   first.set_xlabel('time (d)')
   first.set_ylabel('flux')

   # fitted flux
   time = masked[:,0]
   flux = masked[:,1]
   smth.plot(time, flux, 'k.', ms=3)
   smth.set_xlim(quarterlines[0], 1600)
   smth.set_xlabel('time (d)')
   smth.set_ylabel('normalised flux')

   # fitted flux: zoomed
   if 9 in emptyquarters: # default nice quarter to plot is 9, otherwise pick highest quarter to zoom in on
      qa = quarterlines[max(quarterlist)-1]
      qb = quarterlines[max(quarterlist)]
   else:
      qa = quarterlines[9]
      qb = quarterlines[10]
   zoom.plot(time[(time >= qa) & (time <= qb)], flux[(time >= qa) & (time <= qb)], 'k.', ms=3)
   zoom.set_xlim(qa, qb)
   zoom.set_xlabel('time (d)')
   zoom.set_ylabel('normalised flux')

   # amplitude spectrum
   ft.loglog(freqs, ampls, 'k-', lw=0.5)
   ft.set_xlim(0.1, 283)
   ft.set_ylim(0.1, np.nanmax(maxforplotting)*1.1)
   ft.set_xlabel(r'frequency ($\mu$Hz)')
   ft.set_ylabel('amplitude (ppm)')

   plt.tight_layout()
   os.chdir(f'./{kic}/')
   plt.savefig(f'kic{kic}_factor{factor}_mask{i}_sap_summary.png')
   os.chdir('..')
   plt.close()

else:
   pass

print(f'{kic} done')
