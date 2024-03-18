#!/usr/bin/env python3
import numpy as np

rms1 = 0.306
rms2 = 0.206
Pixel_increment = 0.04
beam = ((np.pi / np.log(2))*(0.25*(rms1 * rms2))) / (Pixel_increment**2)

print(beam)