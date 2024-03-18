#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib import rc, rcParams
from skimage import measure
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

rcParams.update({'font.size': 10})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

hdul = fits.open('/Users/liuxihe/Desktop/我的/final_fits/band6.image.new1.fits')[0]
data = hdul.data * 1e3 

rms = 4.5e-2

levels = np.arange(3.5, 28.5, 5) * rms

header = hdul.header
wcs = WCS(header)

# 获取原始图像的坐标信息
crval1 = header['CRVAL1']
crval2 = header['CRVAL2']
cdelt1 = header['CDELT1']
cdelt2 = header['CDELT2']

# 计算每个像素的原始坐标
ra = (np.arange(hdul.shape[1]) - header['CRPIX1']) * cdelt1 + crval1
dec = (np.arange(hdul.shape[0]) - header['CRPIX2']) * cdelt2 + crval2

mask = data > 25 * rms
y, x = np.where(mask)
coords = np.column_stack((x, y))
values = data[y, x]

fig = plt.figure()
ax = fig.add_subplot(111, projection=wcs)
ax.imshow(data, cmap='inferno', vmin=-0.1623136558741939, vmax=1.7233099994687848,origin='lower')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('ALMA band6 data image with rms contours')
ax.contour(data, levels=levels,linewidths=1,colors='lime')
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
        wcs_coords = wcs.all_pix2world(contour_max_coords.reshape(1, -1), 0)
        contour_max_coord = SkyCoord(wcs_coords[0][0], wcs_coords[0][1], unit='deg', frame='fk5')
        contour_max_coord_list.append(contour_max_coord)
        coord_pix = wcs.all_world2pix([[contour_max_coord.ra.deg, contour_max_coord.dec.deg]], 0)[0]
        coord_x = coord_pix[0]
        coord_y = coord_pix[1]
        coord_x, coord_y = np.flip([coord_x, coord_y])
        contour_max_coord_fk5 = contour_max_coord.transform_to('fk5')
        print(f"（RA,Dec）坐标系：RA={contour_max_coord_fk5.ra.to_string(unit=u.hourangle, sep=':', pad=True)}, Dec={contour_max_coord_fk5.dec.to_string(unit=u.degree, sep=':', pad=True)}")
        ax.plot(coord_x, coord_y, marker='*', markersize=10, color='blue') 
cbar = fig.colorbar(ax.images[0], ax=ax, shrink=0.8, pad=0.05, aspect=30,orientation='horizontal')
#cbar.set_label('mJy/beam')
cbar.ax.set_position([0.12, 0.15, 0.77, 0.03])
cbar.mappable.set_clim(vmin=-1, vmax=2)
ax.set_xlim(475, 1675)
ax.set_ylim(500,1275)
plt.show()
#plt.savefig('/Users/liuxihe/Desktop/ALMA_band6_rms_part.pdf', dpi=1000,bbox_inches='tight')