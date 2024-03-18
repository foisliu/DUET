#!/usr/bin/env python3
import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u

cube = SpectralCube.read('/Users/liuxihe/Desktop/band3.pbcor.fits')
beam = radio_beam.Beam(major=0.309*u.arcsec, minor=0.207*u.arcsec, pa=-85.4661*u.deg)
new_cube = cube.convolve_to(beam)
new_cube.write('/Users/liuxihe/Desktop/band3.pbcor.new.fits')