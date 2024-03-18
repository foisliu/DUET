#!/usr/bin/env python3
def dms_to_dd(d, m, s):
    if d < 0:
        dd = d - m/60 - s/3600
    else:
        dd = d + m/60 + s/3600
    return dd

def hms_to_dd(h, m, s):
    hh = h*15 + m*15/60 + s*15/3600
    return hh

with open ('/Users/liuxihe/Desktop/DUET_data_2/sgrb1off_dendro_cores_band6.dat', 'r') as f:
    data = f.readlines () 

new_data = []

for line in data:
    columns = line.split ()
    ra = columns [1]
    dec = columns [2]
    ra_h, ra_m, ra_s = ra.split(':')
    dec_d, dec_m, dec_s = dec.split(':')
    ra_hh = hms_to_dd (float (ra_h), float (ra_m), float (ra_s))
    dec_dd = dms_to_dd (float (dec_d), float (dec_m), float (dec_s))
    new_data.append (' '.join ([columns [0], str (ra_hh), str (dec_dd)] + columns [3:]))

with open ('band6_data_core_v2.dat', 'w') as f:
    for line in new_data:
        f.write (line + '\n') 


with open ('band6_data_core_v2.dat', 'r') as f:
    print (f.read ())
