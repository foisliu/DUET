#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plfit
from numpy.random import rand,seed
from plfit import plfit

MyPL = plfit(mydata)
MyPL.plotpdf(log=True)

# 读取文件
df = pd.read_csv('/Users/liuxihe/Desktop/DUET_data_2/sgrb1off_dendro_cores_band6.dat', sep='\s+', header=None, usecols=[5], names=['Condensation Masses'])

# 将数据转换为 numpy 数组
mass_data = np.array(df['Condensation Masses'])

# 生成区间边界
bins = np.logspace(np.log10(np.min(mass_data)), np.log10(np.max(mass_data)), 31)

# 定义常量
MASS_SENSITIVITY_FACTOR = 0.27 # x轴上的坐标值

# 筛选出大于等于云核质量下限值的云核质量
mass_data = mass_data[np.where(mass_data >= MASS_SENSITIVITY_FACTOR )]

# 计算每个区间内的云核数量
counts, edges = np.histogram(mass_data, bins=bins)

# 计算区间宽度
widths = edges[1:] - edges[:-1]

# 计算\frac{\mathrm{d} N}{\mathrm{d} \log_{}{M} }
dn_dlog10m = counts / widths

# 计算每个柱子的中心位置
centers = edges[:-1]

# 绘制柱形图
plt.bar(centers, dn_dlog10m, width=widths, align='edge', color='grey', edgecolor='black')

# 设置图形参数
plt.title('Condensation Masses')
plt.xlabel('Condensation Masses ($M_{\odot } $)')
plt.ylabel(r'$\mathrm{d} N/\mathrm{d} \log_{}{M} $')
plt.xscale('log')
plt.yscale('log')
plt.grid()

# 在图中画一条垂直的蓝线
plt.axvline(MASS_SENSITIVITY_FACTOR , color='blue', linewidth=2, linestyle='--', label='Mass sensitivity')

# 在横轴上画一条red sticks
plt.stem(mass_data, 0.01 * np.ones_like(mass_data), linefmt='r-', markerfmt='r-', basefmt='r-', label='Data points')


#显示图形
plt.savefig('/Users/liuxihe/Desktop/core_mass.pdf', dpi=300)