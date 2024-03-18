#!/usr/bin/env python3
import numpy as np

# 定义band3和band6的频率，单位为GHz
band3_freq = 93.5984# 你可以根据实际情况修改这个值
band6_freq = 225.987 # 你可以根据实际情况修改这个值

# 打开aaa.dat和bbb.dat文件，分别读取第四列数据，存储为numpy数组
band3_data = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/band6_band3.dat', usecols=4)
#band3_data = 4.3e-2 # min value of rms
band6_data = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/band3_band6.dat', usecols=4)
#band6_data = 2e-1 # min value of rms

band3_data_error = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/band6_band3.dat', usecols=5)
band6_data_error = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/band3_band6.dat', usecols=5)

# 计算spectral index，使用公式S_nu = nu^alpha
spectral_index = (np.log10(band6_data / band3_data)) / (np.log10(band6_freq / band3_freq))

spectral_index_error = np.log10((band6_data / band3_data) * (np.sqrt((band6_data_error / band6_data) ** 2 + (band3_data_error / band3_data) ** 2)))

spectral_index_error_error = abs(spectral_index_error/((band6_data / band3_data) * np.log(10)))
# 创建一个新的numpy数组，第一列为序号，第二列为spectral index
output_data = np.column_stack((np.arange(1, len(spectral_index) + 1), spectral_index, spectral_index_error_error))

# 输出一个名为fff.dat的文件，使用空格分隔每一列
np.savetxt('spectral_index_band3.dat', output_data, fmt='%d %.4f %.4f', delimiter=' ')
