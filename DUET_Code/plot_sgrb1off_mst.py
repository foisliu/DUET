#!/usr/bin/env python3
import numpy as np
from astropy.io import fits
import aplpy

import matplotlib.pyplot as plt
#import matplotlib.colors as colors

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetex=True)
#rcParams['text.latex.preamble'] = [r'\usepackage{siunitx}',r'\sisetup{detect-all}',r'\usepackage{helvet}',r'\usepackage{sansmath}',r'\sansmath']
rcParams.update({'font.size': 15})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

####
hdu = fits.open('/Users/liuxihe/Desktop/DUET_figure_2/sgrb1off_south_cont_selfcal.commonbeam.image.fits')[0]#使用SgrB1off的band6image
hdu.data = hdu.data * 1e3 # Convert to mJy/beam

dist = 8.178 # kpc

fig = plt.figure(figsize=(10, 5))

mc = [266.6907, -28.532]
f1 = aplpy.FITSFigure(hdu,figure=fig,subplot=[0.185,0.03,0.80,0.92])
f1.show_colorscale(cmap='gist_heat_r', vmax=28, vmin=5e-2, aspect='equal', stretch='log')
f1.recenter(mc[0], mc[1], width=0.029, height=0.018)
# Show primary beams
f1.show_contour('/Users/liuxihe/Desktop/DUET_figure_2/sgrb1off_south_cont_selfcal.pb.fits',levels=[0.3,0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=10)
#f1.show_contour('../../SgrB1off_c1p1/CONT/sgrb1off_c1p1_cont.pb.fits',levels=[0.3,0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=10)
#f1.show_contour('../../SgrB1off_w1/CONT/sgrb1off_w1_cont.pb.fits',levels=[0.3,0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=10,slices='0')

f1.tick_labels.set_xformat('hh:mm:ss')
f1.tick_labels.set_yformat('dd:mm:ss')
f1.ticks.set_xspacing(2./3600.*15.)
f1.ticks.set_yspacing(30./3600.)
f1.axis_labels.set_xpad(0.4)
f1.axis_labels.set_ypad(-0.4)
f1.set_nan_color('white')
f1.ticks.set_color('black')
f1.ticks.set_minor_frequency(4)
f1.ticks.set_length(10)

# Scale bar
f1.add_scalebar(0.1/dist/1e3/np.pi*180,label=r'0.1 pc',corner='bottom right')
f1.scalebar.set_font_size('large')
# Beam
f1.add_beam()
f1.beam.set_frame(True)
f1.beam.set_color('black')
f1.beam.set_linewidth(0.1)
# Color bar
f1.add_colorbar()
f1.colorbar.set_location('right')
#f1.colorbar.set_tick_pad(-34)
f1.colorbar.set_pad(0.05)
f1.colorbar.set_width(0.15)
f1.colorbar.set_ticks([0.05,0.1,0.2,0.5,1,2,5,10,20])
f1.colorbar.set_axis_label_pad(8)
#f1.colorbar.set_axis_label_text('$I_\mathrm{1.3\,mm}$ (Jy/beam)')

f1.add_label(0.05, 0.92, r'Cloud ‘e’ Cores & MST', relative=True, weight='black', size='large', ha='left', zorder=50)

# Plot condensations and MST
ra_arr,dec_arr = np.loadtxt('/Users/liuxihe/Desktop/band6_data_core_v2.dat',usecols=(1,2),unpack=True)#only band6
f1.show_markers(ra_arr,dec_arr,marker='+',s=50 ,ec='blue',fc='blue',linewidths=1.25,alpha=0.75,zorder=10)

edges,ras,decs,dras,ddecs = np.loadtxt('/Users/liuxihe/Desktop/only_band6_v2.dat', usecols=(0,1,2,3,4), unpack=True)
edges = edges * 3600.
#idx_pruned = np.where(edges <= 7.)
#ras_pruned = ras[idx_pruned]
#decs_pruned = decs[idx_pruned]
#dras_pruned = dras[idx_pruned]
#ddecs_pruned = ddecs[idx_pruned]
#f1.show_arrows(ras_pruned,decs_pruned,dras_pruned,ddecs_pruned,head_length=0.,head_width=0.,wid th=4.0,alpha=0.5,color='k')
f1.show_arrows(ras,decs,dras,ddecs,head_length=0.,head_width=0.,width=1.0,alpha=0.5,color='k',zorder=20)

# Show zoom-in regions
w = 0.006
mcz = [266.6948, -28.5334]
#f1.show_rectangles(mcz[0],mcz[1],w,w,color='lime',linewidth=2,linestyle='dashed',zorder=30)
#f1.add_label(mcz[0],mcz[1]+0.0036,'A',size='large',color='green',zorder=31)

f1.save('sgrb1off_MST.pdf', dpi=300)
