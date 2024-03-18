#!/usr/bin/env python3
import numpy as np
from astrodendro import Dendrogram, pp_catalog
from astropy.io import fits
from astropy.modeling.models import BlackBody
from astropy import constants as const
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import aplpy
import os
import matplotlib.pyplot as plt

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetex=True)
#rcParams['text.latex.preamble'] = r'\usepackage{qunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
rcParams.update({'font.size': 24})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'


hdul = fits.open('/Users/liuxihe/Desktop/DUET_figure_2/SgrB1off_band3_cont_rebin30.image.tt0.regrid2b6.new.fits')[0] 
data = hdul.data * 1e3 

dist = 8.1 # kpc
header = hdul.header
#wcs = WCS(header)

#crval1 = header['CRVAL1']
#crval2 = header['CRVAL2']
#cdelt1 = header['CDELT1']
#cdelt2 = header['CDELT2']
#ra = (np.arange(hdul.shape[1]) - header['CRPIX1']) * cdelt1 + crval1
#dec = (np.arange(hdul.shape[0]) - header['CRPIX2'])* cdelt2 + crval2


fig = plt.figure(figsize=(10, 6))
mc = [hdul.header['CRVAL1'], hdul.header['CRVAL2']]
f1 = aplpy.FITSFigure(hdul,figure=fig,subplot=[0.185,0.03,0.80,0.92])
f1.show_colorscale(cmap='inferno', vmax=0.000041257113084711234, vmin=-0.00002673653611230131, aspect='equal')
f1.recenter(mc[0], mc[1], width=0.033, height=0.019)

f1.tick_labels.set_xformat('hh:mm:ss')
f1.tick_labels.set_yformat('dd:mm:ss')
f1.ticks.set_xspacing(1./3600.*15.)
f1.ticks.set_yspacing(10./3600.)
f1.axis_labels.set_xpad(0.4)
f1.axis_labels.set_ypad(-0.4)
f1.set_nan_color('white')
f1.ticks.set_color('black')
f1.ticks.set_minor_frequency(4,2)
f1.ticks.set_length(6)

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
#f1.colorbar.set_ticklabels([0.02,0.05,0.1,0.2,0.5,1])
f1.colorbar.set_ticks([-0.00003,-0.00002,-0.00001,0,0.00001,0.00002,0.00003,0.00004])
f1.colorbar.set_axis_label_pad(15)

f1.set_title('ALMA_band3 after smoothing')

f1.save('/Users/liuxihe/Desktop/ALMA_band3.pdf',dpi=150)
