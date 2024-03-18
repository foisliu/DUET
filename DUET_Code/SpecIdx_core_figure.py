#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches  # 导入patches模块

data1 = np.loadtxt("/Users/liuxihe/Desktop/DUET/DUET_data_2/spectral_index_band36.dat", usecols=(1,))
SpecIdx1 = data1
bins = np.linspace(1.5, 4.5, 40)  # 生成区间边界
counts, edges = np.histogram(SpecIdx1, bins)
plt.bar(edges[:-1], counts, width=0.05, color="grey", align="edge")

data2 = np.loadtxt("/Users/liuxihe/Desktop/DUET/DUET_data_2/spectral_index_band6.dat", usecols=(1,))
SpecIdx2 = data2
bins = np.linspace(1.5, 4.5, 40)  # 生成区间边界
counts, edges = np.histogram(SpecIdx2, bins)
plt.bar(edges[:-1], counts, width=0.05, color="purple", align="edge", alpha = 0.5)

# 创建图例
grey_patch = mpatches.Patch(color='grey', label='Band 3&6')
purple_patch = mpatches.Patch(color='purple', label='Band 6')
plt.legend(handles=[grey_patch, purple_patch], loc='upper right')

plt.xlabel("Cloud core SpecIdx")
plt.ylabel("The number of cloud cores")
plt.title("Statistical Histogram of Core SpecIdx")
plt.savefig('/Users/liuxihe/Desktop/SpecIdx_core.pdf', dpi=300)