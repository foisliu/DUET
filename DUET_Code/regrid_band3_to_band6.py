#!/usr/bin/env python3
import numpy as np
from astropy.io import fits
#from astropy.wcs import WCS
from reproject import reproject_interp
#import matplotlib.pyplot as plt

b3im = '/Users/liuxihe/Desktop/DUET_figure_2/SgrB1off_band3_cont_rebin30v2.pb.tt0.fits'
b6im = '/Users/liuxihe/Desktop/DUET_figure_2/sgrb1off_south_cont_selfcal.commonbeam.image.pbcor.fits'

map0 = fits.open(b6im)[0]
i0 = map0.data
h0 = map0.header

map1 = fits.open(b3im)[0]
i1 = map1.data
h1 = map1.header

#h0['bmaj'] = h1['bmaj']
#h0['bmin'] = h1['bmin']
#h0['bpa'] = h1['bpa']

# Reproject Parsons data to CMZPol2 data
newi1, footprint = reproject_interp(map1, h0)

fits.writeto('SgrB1off_band3_cont_rebin30.image.tt0.regrid2b6.new.pb.fits', newi1, h0, overwrite=True)
