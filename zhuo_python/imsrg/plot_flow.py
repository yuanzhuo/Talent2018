#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:48:23 2018

@author: zhuo
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from matplotlib import axes
from mpl_toolkits.mplot3d import Axes3D



i_list = [0,1,2,3,4,5,6,7,8,9,10]#,10,15,20]
i_list = [0,1,2,3,4,5,10,20,40,60,80]#,10,15,20]
i_list = [0,1,2,3,4,5,10,20,30,40]
#i_list = [0,1]

# for i in i_list:
#     F_1 = "data/srg_conv_" + str(i) + ".txt"
#     data_F = np.loadtxt(F_1)
#     x_len = len(data_F)
#     x = np.arange(0,x_len)
#     y = x
#     plt.contourf(x, y, data_F[x][y], 8, alpha = 0.1, cmap = plt.cm.hot)
#     out_f = "data/temp" + str(i) +".eps"
#     plt.savefig(out_f,format='eps')
#
# for i in i_list:
#     F_1 = "data/srg_conv_f_" + str(i) + ".txt"
#     data_F = np.loadtxt(F_1)
#     x_len = len(data_F)
#     x = np.arange(0,x_len)
#     y = x
#     plt.contourf(x, y, data_F[x][y], 8, alpha = 0.1, cmap = plt.cm.hot)
#     out_f = "data/temp_f" + str(i) +".eps"
#     plt.savefig(out_f,format='eps')

# for i in i_list:
#     F_1 = "data/srg_conv_G_" + str(i) + ".txt"
#     data_F = np.loadtxt(F_1)
#     x_len = len(data_F)
#     x = np.arange(0,x_len)
#     y = x
#     plt.contour(x, y, data_F[x][y], 12, alpha = 1, cmap = plt.cm.Blues)
#     out_g = "data/temp_G_Test" + str(i) +".eps"
#     plt.savefig(out_g,format='eps')
#     print(out_g)

figure = plt.figure()

ax = Axes3D(figure)

X = np.arange(0, 256)

Y = X
F_1 = "data/srg_conv_G_0.txt"
data_F = np.loadtxt(F_1)
ax.plot_surface(X, Y, data_F, rstride=1, cstride=1, cmap='rainbow')

plt.show()


