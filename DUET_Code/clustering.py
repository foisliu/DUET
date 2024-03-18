#!/usr/bin/env python3
import numpy as np

def dist(coord1, coord2):
    ra1,dec1,ra2,dec2 = coord1[0]*np.pi/180,coord1[1]*np.pi/180,coord2[0]*np.pi/180,coord2[1]*np.pi/180
    return np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2))/np.pi*180.

distance = 8.1 # kpc

# read in condensation coordinates
ras,decs = np.loadtxt('/Users/liuxihe/Desktop/only_band6_v2.dat',usecols=(1,2),unpack=True)
pos = np.column_stack((ras, decs))

# Edges are in degree
edges = np.loadtxt('/Users/liuxihe/Desktop/only_band6_v2.dat',usecols=(0),unpack=True)

# Derive the normalized mean edge length of the MST, m_bar
Nc = len(ras)

mean_ra  = (max(ras) + min(ras)) / 2.
mean_dec = (max(decs) + min(decs)) / 2.
meancore = [mean_ra, mean_dec]

dist_cores = [dist(meancore, corei) for corei in pos]

Rc = max(dist_cores)
area = np.pi * Rc**2

mbar = np.sum(edges) / np.sqrt(Nc * area)

# Derive the normalized correlation length, s_bar
all_dist = []
for idx in np.arange(Nc):
    for idy in np.arange(idx+1,Nc):
        onedist = dist(pos[idx], pos[idy])
        all_dist.append(onedist)

#print(len(all_dist))
#print(Nc*(Nc-1)/2.)

sbar = np.mean(all_dist) / Rc

print("m_bar is:", mbar)
print("s_bar is:", sbar)
print("Q factor is:", mbar / sbar)
