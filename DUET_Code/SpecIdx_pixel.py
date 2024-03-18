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
rcParams.update({'font.size': 24})
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

dist = 8.1 # kpc

# 读取图像数据，假设它们存储在numpy数组中
hdu1 = fits.open('/Users/liuxihe/Desktop/DUET_figure_2/SgrB1off_band3_cont_rebin30.image.tt0.regrid2b6.new.pbcor.fits')[0] # 频率为nu1的图像
hdu1.data = hdu1.data * 1e3 # 获取img1中的数据数组
hdu2 = fits.open('/Users/liuxihe/Desktop/DUET_figure_2/sgrb1off_south_cont_selfcal.commonbeam.image.pbcor.fits')[0] # 频率为nu2的图像
hdu2.data = hdu2.data * 1e3 # 获取img2中的数据数组
hdu_pb=fits.open('/Users/liuxihe/Desktop/DUET_figure_2/sgrb1off_south_cont_selfcal.pb.fits')[0].data

# 获取fits文件的头文件信息
header1 = hdu1.header
header2 = hdu2.header

# Define the frequencies of the two bands
nu1 = 93598.4 # Frequency of image1 in MHz
nu2 = 225987 # Frequency of image2 in MHz

# 计算图像的rms参数
rms1 = 4.5e-2
rms2 = 9.7e-3

levels = np.arange(3, 33, 5) * rms1

# 设置rms阈值
threshold1 = 2*rms1
threshold2 = 2*rms2

# 过滤低于阈值的像素点
hdu1.data[hdu1.data < threshold1] = 0
hdu2.data[hdu2.data < threshold2] = 0


# 对图像数据进行对数变换
log_img1 = np.log10(hdu1.data) # 对数据进行对数变换
log_img2 = np.log10(hdu2.data) # 对数据进行对数变换

# 水平翻转和竖直翻转图像数据
#log_img1 = np.flipud(log_img1)
#log_img2 = np.flipud(log_img2)

# 计算每个像素点处的光谱指数 
spectral_index = (log_img1 - log_img2) / np.log10(nu1 / nu2)#nu1和nu2分别为两个波段的频率

####################################################
#Draw the spectral index figure
#ax = plt.gca()
#ax.invert_xaxis()
#ax.invert_yaxis()
fig = plt.figure(figsize=(10, 6))
mc = [hdu1.header['CRVAL1'], hdu1.header['CRVAL2']]
f1 = aplpy.FITSFigure(hdu1,figure=fig,subplot=[0.255,0.09,0.70,0.82])
f1.recenter(mc[0]+ 0.005, mc[1]-0.001, width=0.007, height=0.0006)


f1.show_contour(hdu_pb,levels=[0.3,0.5],colors='gold',linewidths=2.0,linestyles='dashed',zorder=30)
f1.tick_labels.set_xformat('hh:mm:ss')
f1.tick_labels.set_yformat('dd:mm:ss')
f1.ticks.set_xspacing(2./3600.*15.)
f1.ticks.set_yspacing(20./3600.)
f1.axis_labels.set_xpad(0.4)
f1.axis_labels.set_ypad(1.8)
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

f1.show_contour(hdu1.data, levels=levels, linewidths=0.5, colors='black')

plt.imshow(spectral_index, cmap="jet", vmin=2, vmax=4.5)#extent为坐标轴范围，vmin和vmax为谱指数的最小和最大值
plt.xlabel("Right Ascension (deg)")
plt.ylabel("Declination (deg)")
plt.title("Spectral index map at ALMA band 3 and 6")
#plt.xlim([266.700, 266.690])
#plt.ylim([-28.5375, -28.5300])
cbar = plt.colorbar(location='right', shrink=0.8, aspect=30)
#cbar.ax.set_ylabel('Spectral index')
plt.savefig('/Users/liuxihe/Desktop/SpecIdx_pixel.pdf', dpi=500)
#plt.show()