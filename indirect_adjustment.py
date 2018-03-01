# -*- coding: utf-8 -*-

import numpy as np

# 基线数据
line = [[0 for i in range(5)] for j in range(24)]
line[0] = ['4', '1', -1632.834, -445.751, -1.070]     #L1
line[1] = ['1', '2', -519.547, 2730.001, 1.087]       #L2
line[2] = ['2', '5', 1744.190, 37.123, 1.888]         #L3
line[3] = ['2', '6', 1484.372, 2952.909, 1.329]       #L4
line[4] = ['2', '3', -371.070, 3287.184, 0.901]       #L5
line[5] = ['6', '3', -1855.484, 334.284, -0.437]
line[6] = ['3', '7', 1391.588, 2302.878, -0.120]
line[7] = ['8', '4', -1410.208, 761.462, 0.528]
line[8] = ['4', '5', -408.232, 2321.422, 1.917]
line[9] = ['5', '9', 1966.853, -408.542, -0.458]
line[10] = ['6', '5', 259.748, -2915.705, 0.559]
line[11] = ['10', '6', -1410.776, 993.282, -0.609]
line[12] = ['6', '11', 1410.785, 1006.710, 2.208]
line[13] = ['7', '6', 463.895, -2637.183, 0.527]
line[14] = ['7', '11', 1874.660, -1630.448, 2.757]
line[15] = ['12', '8', -1280.313, -1652.853, -0.558]
line[16] = ['9', '8', -148.436, -2674.318, -2.021]
line[17] = ['9', '12', 1131.864, -1021.458, -1.436]
line[18] = ['13', '9', -1057.654, -1040.031, 1.566]
line[19] = ['9', '10', -815.796, 2331.041, 0.521]
line[20] = ['10', '14', 1687.875, 343.308, -0.982]
line[21] = ['11', '14', 1687.916, -1656.719, -2.569]
line[22] = ['13', '12', 74.239, -2061.476, 0.139]
line[23] = ['14', '13', 185.542, -1634.288, -1.129]

# 点位数据，将P1(0, 0, 0)置为已知点
point = [[0 for i in range(3)] for j in range(14)]
point[0] = [0, 0, 0]
point[1] = [x + y for x, y in zip(point[0], line[1][2:])]   #P2
point[2] = [x + y for x, y in zip(point[1], line[4][2:])]   #P3
point[3] = [x - y for x, y in zip(point[0], line[0][2:])]   #P4
point[4] = [x + y for x, y in zip(point[1], line[2][2:])]   #P5
point[5] = [x + y for x, y in zip(point[1], line[3][2:])]
point[6] = [x + y for x, y in zip(point[2], line[6][2:])]
point[7] = [x - y for x, y in zip(point[3], line[7][2:])]
point[8] = [x + y for x, y in zip(point[4], line[9][2:])]
point[9] = [x + y for x, y in zip(point[8], line[19][2:])]
point[10] = [x + y for x, y in zip(point[5], line[12][2:])]
point[11] = [x + y for x, y in zip(point[8], line[17][2:])]
point[12] = [x - y for x, y in zip(point[8], line[18][2:])]
point[13] = [x + y for x, y in zip(point[9], line[20][2:])]

# 从基线数据中提取系数构造B矩阵
matrix_B = [[0 for i in range(14)] for j in range(24)]

# 将基线中起点位置对应的系数置为-1，终点位置对应的系数置为1
for i in range(0, 24):
    matrix_B[i][int(line[i][1]) - 1] = 1
    matrix_B[i][int(line[i][0]) - 1] = -1

# 提取基线数据中的x, y, z增量，方便后续计算
line_data = [line[iter][2:5] for iter in range(0,24)]

# 计算误差方程常数项
l = line_data - np.dot(np.array(matrix_B), np.array(point))

# 等权观测，构造单位对角阵
matrix_P = np.eye(24)

# 计算法方程的解x
Nbb = np.dot(np.dot(np.transpose(matrix_B), matrix_P), matrix_B)
W = np.dot(np.dot(np.transpose(matrix_B), matrix_P), l)
x = np.dot(np.linalg.inv(Nbb), W)

# 将x带入误差方程，求得改正数
V = np.dot(matrix_B, x) - l

# 由基线观测值与改正数算得平差值
result = line_data + V

# 输出改正数矩阵
print '\n改正数：'
print V

# 输出平差值
print '\n平差值：'
print result