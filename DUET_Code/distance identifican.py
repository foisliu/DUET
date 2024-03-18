#!/usr/bin/env python3
# 导入所需的库
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord,Angle


# 读取aaa.dat和bbb.dat两个文件，分别存储为numpy数组
band3 = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/sgrb1off_dendro_cores_band3.dat',dtype=np.str_,encoding='utf-8')
band6 = np.loadtxt('/Users/liuxihe/Desktop/DUET/DUET_data_2/sgrb1off_dendro_cores_band6.dat',dtype=np.str_,encoding='utf-8')
# 获取两个文件的行数和列数
n3, m3 = band3.shape
n6, m6 = band6.shape

# 定义一个函数，用来计算两个坐标之间的距离
def distance(ra1, dec1, ra2, dec2):
    # 创建SkyCoord对象，使用ICRS坐标系
    c1 = SkyCoord(ra=ra1*u.deg, dec=dec1*u.deg, frame='icrs')
    c2 = SkyCoord(ra=ra2*u.deg, dec=dec2*u.deg, frame='icrs')
    # 计算两个坐标之间的角距离，并转换为弧秒
    dist = c1.separation(c2).to(u.arcsec)
    # 返回距离值
    return dist.value

# 创建一个空列表，用来存储合并后的云核数据
merged = []

# 遍历band3的每一行
for i in range(n3):
    # 获取band3的第i行的最后一列，作为peak坐标
    ra1 = band3[i, -2]
    dec1 = band3[i, -1]
    #print(ra1,dec1)
    ra1 = Angle(ra1).degree
    dec1 = Angle(dec1).degree
    # 创建一个布尔变量，用来标记是否找到与band3匹配的band6云核
    found = False
    # 遍历band6的每一行
    for j in range(n6):
        # 获取band6的第j行的最后一列，作为peak坐标
        ra2 = band6[j, -2]
        dec2 = band6[j, -1]
        ra2 = Angle(ra2).degree
        dec2 = Angle(dec2).degree
        #print(ra1, dec1, ra2, dec2)
        # 计算两个peak坐标之间的距离
        dist = distance(ra1, dec1, ra2, dec2)
        # 如果距离小于0.253弧秒，认为是同一个云核
        if dist < 0.251:
            # 将band6的第j行数据添加到合并列表中，并在最后一列添加"band3_band6"标记
            merged.append(np.append(band6[j], "band3_band6"))
            merged.append(np.append(band3[i], "band6_band3"))
            # 将found变量设为True，表示找到了匹配的云核
            found = True
            # 跳出内层循环，继续下一个band3云核的匹配
            break
    # 如果没有找到匹配的云核，将band3的第i行数据添加到合并列表中，并在最后一列添加"only_band3"标记
    if not found: 
        merged.append(np.append(band3[i],"only_band3"))

# 遍历band6的每一行
for j in range(n6):
    # 获取band6的第j行的最后一列，作为peak坐标
    ra2 = band6[j, -2]
    dec2 = band6[j, -1]
    #print(ra1,dec1)
    ra2 = Angle(ra2).degree
    dec2 = Angle(dec2).degree
    # 创建一个布尔变量，用来标记是否找到与band3匹配的band6云核
    found = False
    # 遍历band6的每一行
    for k in range(n3):
        # 获取band6的第j行的最后一列，作为peak坐标
        ra1 = band3[k, -2]
        dec1 = band3[k, -1]
        ra1 = Angle(ra1).degree
        dec1 = Angle(dec1).degree
        #print(ra1, dec1, ra2, dec2)
        # 计算两个peak坐标之间的距离
        dist = distance(ra1, dec1, ra2, dec2)
        # 如果距离小于0.253弧秒，认为是同一个云核
        if dist < 0.251:
            # 将band6的第j行数据添加到合并列表中，并在最后一列添加"band3_band6"标记
            #merged.append(np.append(band6[j], "band3_band6"))
            #merged.append(np.append(band3[i], "band6_band3"))
            # 将found变量设为True，表示找到了匹配的云核
            found = True
            # 跳出内层循环，继续下一个band3云核的匹配
            break
    # 如果没有找到匹配的云核，将band6的第j行数据添加到合并列表中，并在最后一列添加"only_band6"标记
    if not found: 
        merged.append(np.append(band6[j],"only_band6"))
# 将合并列表转换为numpy数组
merged = np.array(merged)

# 保存合并数组到ccc.dat文件中
np.savetxt('ccc.dat', merged, fmt='%s')

# Load data from ccc.dat file
data = np.loadtxt("ccc.dat", dtype=str)

# Filter rows where last column contains "band3_band6"
filtered_data = data[data[:, -1] == "band3_band6"]

# Save filtered data to ddd.dat file
np.savetxt("band3_band6.dat", filtered_data, fmt="%s")

# Load data from ccc.dat file
data = np.loadtxt("ccc.dat", dtype=str)

# Filter rows where last column contains "band3_band6"
filtered_data = data[data[:, -1] == "band6_band3"]

# Save filtered data to ddd.dat file
np.savetxt("band6_band3.dat", filtered_data, fmt="%s")

# Load data from ccc.dat file
data = np.loadtxt("ccc.dat", dtype=str)

# Filter rows where last column contains "band3_band6"
filtered_data = data[data[:, -1] == "only_band3"]

# Save filtered data to ddd.dat file
np.savetxt("only_band3.dat", filtered_data, fmt="%s")

# Load data from ccc.dat file
data = np.loadtxt("ccc.dat", dtype=str)

# Filter rows where last column contains "band3_band6"
filtered_data = data[data[:, -1] == "only_band6"]

# Save filtered data to ddd.dat file
np.savetxt("only_band6.dat", filtered_data, fmt="%s")
