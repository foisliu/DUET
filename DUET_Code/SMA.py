#!/usr/bin/env python3
import numpy as np
from astropy.io import fits
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy import constants as const
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import aplpy
import os
import matplotlib.pyplot as plt
from skimage import measure
from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetex=True)
#rcParams['text.latex.preamble'] = r'\usepackage{qunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
rcParams.update({'font.size': 24})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

hdul = fits.open('/Users/liuxihe/Desktop/我的/sgrb1.cont.image.fits')[0] 
data = hdul.data * 1e3 
dist = 8.1 # kpc

rms = 4.5

levels = np.arange(3, 28, 5) * rms
header = hdul.header
w = WCS(hdul.header)

crval1 = header['CRVAL1']
crval2 = header['CRVAL2']
cdelt1 = header['CDELT1']
cdelt2 = header['CDELT2']
ra = (np.arange(hdul.shape[1]) - header['CRPIX1']) * cdelt1 + crval1
dec = (np.arange(hdul.shape[0]) - header['CRPIX2'])* cdelt2 + crval2

mask = data > 25 * rms
y, x = np.where(mask)
coords = np.column_stack((x, y))
values = data[y, x]

fig = plt.figure()

ax = fig.add_subplot(111, projection=w)

ax.imshow(data, cmap='inferno', vmin=-16.229334392352256, vmax=95.32722747628739,origin='lower')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('SMA image with rms contours')
ax.contour(data, levels=levels,linewidths=2,colors='lime')
ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
# 找到大于等于阈值的闭合区域
contours = measure.find_contours(data, levels[-1])
contour_max_coord_list = []  
for contour in contours:
    if data[int(contour[0, 0]), int(contour[0, 1])] >= levels[-1]:
        # 找到区域内值最大位置的像素坐标
        contour_mask = np.zeros_like(data, dtype=bool)
        contour_mask[np.round(contour[:, 0]).astype(int), np.round(contour[:, 1]).astype(int)] = True
        contour_values = data[contour_mask]
        contour_max_idx = np.argmax(contour_values)
        contour_max_coords = np.argwhere(contour_mask)[contour_max_idx]
        wcs_coords = w.all_pix2world(contour_max_coords.reshape(1, -1), 0)
        contour_max_coord = SkyCoord(wcs_coords[0][0], wcs_coords[0][1], unit='deg', frame='fk5')
        contour_max_coord_list.append(contour_max_coord)
        coord_pix = w.all_world2pix([[contour_max_coord.ra.deg, contour_max_coord.dec.deg]], 0)[0]
        coord_x = coord_pix[0]
        coord_y = coord_pix[1]
        coord_x, coord_y = np.flip([coord_x, coord_y])
        contour_max_coord_fk5 = contour_max_coord.transform_to('fk5')
        print(f"（RA,Dec）坐标系：RA={contour_max_coord_fk5.ra}, Dec={contour_max_coord_fk5.dec}")
        ax.plot(coord_x, coord_y, marker='*', markersize=20, color='blue') 
cbar = fig.colorbar(ax.images[0], ax=ax, shrink=0.8, pad=0.05, aspect=30,orientation='horizontal')
#cbar.set_label('mJy/beam')
cbar.ax.set_position([0.12, 0.15, 0.77, 0.03])
cbar.mappable.set_clim(vmin=-20, vmax=100)
#ax.set_xlim(125, 175)
#ax.set_ylim(75,125)



plt.savefig('/Users/liuxihe/Desktop/SMA_rms_part.pdf', dpi=1000,bbox_inches='tight')


#x, y = wcs.all_world2pix(ra, dec, 0)
#skycoord = SkyCoord.from_pixel(np.arange(header['NAXIS1']), np.arange(header['NAXIS2']), wcs=wcs)

#data = np.flipud(data)
#data = np.fliplr(data)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection=wcs)
#ax.imshow(data, cmap='inferno')
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_title('SMA data image with rms contours')
#ax.contour(data, levels=levels,linewidths=0.5,colors='lime')
#ax.set_xlabel('Galactic Longitude (deg)')
#ax.set_ylabel('Galactic Latitude (deg)')
#ax.invert_xaxis()
#ax.invert_yaxis()
#plt.savefig('/Users/liuxihe/Desktop/SMA_rms.pdf', dpi=500,bbox_inches='tight')