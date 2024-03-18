#!/usr/bin/env python3
import numpy as np
from astrodendro import Dendrogram, pp_catalog
from astropy.io import fits
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy import constants as const
from astropy import wcs
from astropy.coordinates import SkyCoord
import aplpy
import os
import matplotlib.pyplot as plt

from matplotlib import rc, rcParams
from matplotlib.font_manager import fontManager, FontProperties
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetex=True)
#rcParams['text.latex.preamble'] = r'\usepackage{qunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
rcParams.update({'font.size':24})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

# whether to overwrite existing leave catalogs
owcat = True
# whether to 'prune' leaves (e.g., whose area above 5sigma is smaller than one beam)
flagleaf = True
# whether to label the leave indices in the plot
labelleaf = True

####
hdu = fits.open('/Users/liuxihe/Desktop/DUET/DUET_figure_2/sgrb1off_south_cont_selfcal.commonbeam.image.fits')[0]
hdu_pbcor = fits.open('/Users/liuxihe/Desktop/DUET/DUET_figure_2/sgrb1off_south_cont_selfcal.commonbeam.image.pbcor.fits')[0]
hdu.data = hdu.data * 1e3 # Convert to mJy/beam
hdu_pbcor.data = hdu_pbcor.data * 1e3 # Convert to mJy/beam
hdu_pb=fits.open('/Users/liuxihe/Desktop/DUET/DUET_figure_2/sgrb1off_south_cont_selfcal.pb.fits')[0].data

dist = 8.1 # kpc

#### Whether to compute dendrogram & core properties
if (not os.path.isfile('sgrb1off_dendro_cores_deg.dat')) or owcat:
    # min_value/min_delta = 4/1 rms
    # min_npix = pixels in one beam, pi * bmaj * bmin / 4 / (pixel)^2
    rms = 4e-2
    beam = 44.64091467190921
    d = Dendrogram.compute(hdu.data, min_value=4*rms, min_delta=1*rms, min_npix=beam)
    
    if flagleaf:
        leaves_repro = []
        for leaf in d.leaves:
            indices = leaf.indices()
            N_pix = (hdu.data[indices]>5*rms).sum()
            if N_pix > beam:
                leaves_repro.append(leaf)
    else:
        leaves_repro = d.leaves
    
    # Drop leaves outside of the PB
    leaves_repro = []
    for leaf in d.leaves:
        peak, value = leaf.get_peak()
        if hdu_pb[peak] >= 0.3:
            leaves_repro.append(leaf)
    
    print('Number of leaves: {0}'.format(np.size(leaves_repro)))
    
    metadata = {}
    metadata['data_unit'] = u.Jy / u.beam
    metadata['spatial_scale'] =  0.04 * u.arcsec
    metadata['beam_major'] =  0.306 * u.arcsec
    metadata['beam_minor'] =  0.206 * u.arcsec
    cat = pp_catalog(leaves_repro, metadata)

    fluxes = []
    for leaf in leaves_repro:
        indices = leaf.indices()
        flux = hdu_pbcor.data[:,:][indices].sum() / beam
        fluxes.append(flux)
    fluxes = np.array(fluxes)
    area   = cat['area_exact'].data # in arcsec^2
    area_beam = area / (np.pi * metadata['beam_major'].value * metadata['beam_minor'].value / 4. / np.log(2))
    fluxes_err = rms * np.sqrt(area_beam)

    
    # Coordinates of leaves
    coords = []
    radecs = []
    radec_deg = []
    # Get the peak coordinates of the cores
    peaks = []
    for leaf in leaves_repro:
        peak = leaf.get_peak()[0]
        peaks.append(peak)
    print(peaks)
    w = wcs.WCS(hdu.header)
    for x,y in zip(list(cat['x_cen']),list(cat['y_cen'])):
        coord = w.all_pix2world(x, y, 1)
        coords.append(coord)
        radec_deg.append(str(coord[0])+' '+str(coord[1]))
        c = SkyCoord(coord[0], coord[1], frame='fk5', unit='deg')
        #radecs.append(c.to_string('hmsdms'))
        rah = int(c.ra.hms.h); ram = int(c.ra.hms.m); ras = c.ra.hms.s
        decd = int(c.dec.dms.d); decm = int(c.dec.dms.m); decs = c.dec.dms.s
        if decd <= 0:
            #decd = '$-$'+str(abs(decd)).zfill(2)
            decd = '-'+str(abs(decd)).zfill(2)
            decm = str(abs(decm)).zfill(2)
            decs = '{:.2f}'.format(abs(decs)).zfill(5)
        newstr = str(rah)+':'+str(ram).zfill(2)+':'+'{:.2f}'.format(ras).zfill(5)+'  '+decd+':'+decm+':'+decs
        radecs.append(newstr)
    
    # Convert pixel coordinates to world coordinates (RA and Dec)
    wcs_coords = w.all_pix2world(peaks, 0)
    sky_coords = SkyCoord(ra=wcs_coords[:, 0], dec=wcs_coords[:, 1], unit=(u.deg, u.deg), frame='fk5')

    # Calculate masses
    mjy = u.mJy.cgs.scale
    pc = const.pc.cgs.value
    msun = const.M_sun.cgs.value
    mp = const.m_p.cgs.value
    distance = dist * 1e3 * pc
    # Effective radius
    ra = np.sqrt(area / np.pi)/3600./180.*np.pi * distance / pc # in pc
    ##
    freq = 225.987 # GHz 
    tdust = 20. # K
    kappa = 0.899 # cm2 g-1
    bb = BlackBody(temperature = tdust*u.K)
    BT = bb(freq*u.GHz).cgs.value
    masses = 1e2 * fluxes * mjy * distance**2 / BT / kappa / msun
    densities = masses*msun / (4./3. * np.pi * (ra * pc)**3) / (2.8 * mp)
    ## The leaves are ordered from bottom to top of the map
    # This is for papers, coordinates are LaTeX formatted
    with open('sgrb1off_dendro_cores.dat', 'w') as newfile:
        for row in list(zip(np.arange(len(radecs)),radecs,np.round(ra*206265,decimals=-1),fluxes,fluxes_err,masses,densities)):
            print("{:>4} {: >26} {:.0f} {:.2f} {:.4f} {:.2f} {:.1e}".format(*row), file=newfile)
    # This is for codes, coordinates are in degrees
    with open('sgrb1off_dendro_cores_deg.dat', 'w') as newfile:
        for row in list(zip(list(cat['_idx']),radec_deg,ra*206265,fluxes,fluxes_err,masses,densities)):
            print("{: >3} {: >30} {:.0f} {:.2f} {:.4f} {:.2f} {:.2e}".format(*row), file=newfile)
    # Read the existing data from the file
    data = []
    with open('sgrb1off_dendro_cores.dat', 'r') as f:
        for line in f:
            data.append(line.strip().split())

    # Add the peaks to the data
    for i, row in enumerate(data):
        row.append(str(sky_coords[i].to_string('hmsdms')))

    # Write the updated data to the file
    with open('sgrb1off_dendro_cores.dat', 'w') as f:
        for i, leaf in enumerate(leaves_repro):
            #f.write(f"{i+1}\t{leaf.get_peak()[0]}\t{sky_coords[i].to_string('hmsdms')}\n")
            f.write('\t'.join(data[i]) + '\n')
else:
    print('Core catalogs already exist!')

####################################################
# Create empty mask. For each leaf we do an 'or' operation with the mask so
# that any pixel corresponding to a leaf is set to True.
if (not os.path.isfile('sgrb1off_dendro_mask.fits')) or owcat:
    mask = np.zeros(hdu.data.shape, dtype=bool)
    for leaf in leaves_repro:
        mask = mask | leaf.get_mask()
    
    newmask = np.array(mask, dtype=int)
    fits.writeto('sgrb1off_dendro_mask.fits', newmask, hdu.header, overwrite=True)

    index_only6 = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/only_band6.dat', usecols=(0,), dtype=int)
    leaves_only6 = []
    for i, leaf in enumerate(leaves_repro):
        if i in index_only6:
            leaves_only6.append(leaf)
    mask_only6 = np.zeros(hdu.data.shape, dtype=bool)
    for leaf in leaves_only6:
        mask_only6 = mask_only6 | leaf.get_mask()
    newmask_only6 = np.array(mask_only6, dtype=int)
    fits.writeto('/Users/liuxihe/Desktop/sgrb1off_dendro_mask_only_band6.fits', newmask_only6, hdu.header, overwrite=True)

    index_band6 = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/band3_band6.dat', usecols=(0,), dtype=int)
    leaves_band6 = []
    for i, leaf in enumerate(leaves_repro):
        if i in index_band6:
            leaves_band6.append(leaf)
    mask_band6 = np.zeros(hdu.data.shape, dtype=bool)
    for leaf in leaves_band6:
        mask_band6 = mask_band6 | leaf.get_mask()
    newmask_band6 = np.array(mask_band6, dtype=int)
    fits.writeto('/Users/liuxihe/Desktop/sgrb1off_dendro_mask_band3_band6.fits', newmask_band6, hdu.header, overwrite=True)

    # Now we create a FITS HDU object to contain this, with the correct header
    mask_hdu = fits.PrimaryHDU(mask.astype('short'), hdu.header)
else:
    mask_hdu = fits.open('sgrb1off_dendro_mask.fits')[0]
    print('Core mask file already exists!')

####################################################
# We then use APLpy to make the final plot
fig = plt.figure(figsize=(10, 6))
#mc = [266.6907, -28.5187]
mc = [hdu.header['CRVAL1'], hdu.header['CRVAL2']]
f1 = aplpy.FITSFigure(hdu,figure=fig,subplot=[0.185,0.03,0.80,0.92])
f1.show_colorscale(cmap='cubehelix_r', vmax=37, vmin=4e-1, aspect='equal', stretch='log')
f1.recenter(mc[0], mc[1], width=0.029, height=0.016)
# Show primary beams
#f1.show_contour('sgrb1off_south_cont_selfcal.pb.fits',levels=[0.3,0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=30)
f1.show_contour(hdu_pb,levels=[0.3,0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=30)

mask_hdu_only3 = fits.open('/Users/liuxihe/Desktop/sgrb1off_dendro_mask_only_band3.fits')[0]
mask_hdu_only6 = fits.open('/Users/liuxihe/Desktop/sgrb1off_dendro_mask_only_band6.fits')[0]
mask_hdu_band6 = fits.open('/Users/liuxihe/Desktop/sgrb1off_dendro_mask_band3_band6.fits')[0]
f1.show_contour(mask_hdu_only3, colors='red', linewidths=1)
f1.show_contour(mask_hdu_only6, colors='cyan', linewidths=1)
f1.show_contour(mask_hdu_band6, colors='blue', linewidths=1)

f1.tick_labels.set_xformat('hh:mm:ss')
f1.tick_labels.set_yformat('dd:mm:ss')
f1.ticks.set_xspacing(2./3600.*15.)
f1.ticks.set_yspacing(20./3600.)
f1.axis_labels.set_xpad(0.4)
f1.axis_labels.set_ypad(0)
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
#f1.colorbar.set_ticklabels([0.1,0.2,0.5,1,2,5,10,20])
f1.colorbar.set_ticks([0.4,0.8,1.6,6.4,32,])
f1.colorbar.set_axis_label_pad(10)
#f1.colorbar.set_axis_label_text('$I_\mathrm{1.3\,mm}$ (mJy/beam)')

f1.add_label(0.03, 0.93, r'Cloud e Cores', relative=True, weight='black', size='large', ha='left', zorder=50)

f1.set_title('$I_\mathrm{1.3\,mm}$ (mJy/beam)')

#if labelleaf:
    #ids, ras, decs = np.loadtxt('sgrb1off_dendro_cores_deg.dat', usecols=(0,1,2), unpack=True)
    #rcParams.update({'font.size': 20})
    #for idx, ra, dec in zip(ids, ras, decs):
        #f1.add_label(ra, dec, str(int(idx)), weight='black', size='xx-small', ha='left', zorder=30)
#plt.show()
f1.save('sgrb1off_cores.pdf',dpi=150)
