#!/usr/bin/env python3
import numpy as np
import FragMent

#target = 'sgrb1off'

# read in condensation coordinates
ras,decs = np.loadtxt('/Users/liuxihe/Desktop/band6_data_core_v2.dat',usecols=(1,2),unpack=True)
#ras,decs = np.loadtxt('../dendro/'+target+'_dendro_cores_deg_3sigma.dat',usecols=(1,2),unpack=True)
#ras,decs = np.loadtxt('../dendro/'+target+'_dendro_cores_deg_sig2sigma.dat',usecols=(1,2),unpack=True)

pos = np.column_stack((ras, decs))

mst_seps, mst = FragMent.MST(pos)

idx, idy = np.where(mst > 0)

posarray = np.hstack((pos[idx],pos[idy]))

delta_ras = [pair[2]-pair[0] for pair in posarray]
delta_decs = [pair[3]-pair[1] for pair in posarray]

newras = ras[idx]
newdecs = decs[idx]

with open('only_band6_v2.dat','w') as fout:
#with open(target+'_mst_log_3sigma.dat','w') as fout:
#with open(target+'_mst_log_sig2sigma.dat','w') as fout:
    for row in list(zip(mst_seps,newras,newdecs,delta_ras,delta_decs)):
        print("{:.15f} {:.10f} {:.10f} {:.10f} {:.10f}".format(*row), file=fout)