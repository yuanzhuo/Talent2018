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
import math


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
# figure = plt.figure()
#
# ax = Axes3D(figure)
#
# X = np.arange(0, 256)
#
# Y = X
# F_1 = "data/srg_conv_G_0.txt"
# data_F = np.loadtxt(F_1)
# ax.plot_surface(X, Y, data_F, rstride=1, cstride=1, cmap='rainbow')
#
# plt.show()

def draw_heatmap(data,xlabels,ylabels,index):

    #cmap = cm.Blues
    cmap = cm.binary
    
    #cmap = cm.get_cmap('rainbow',10)


    figure=plt.figure(facecolor='w')

    ax=figure.add_subplot(2,1,1,position=[0.1,0.15,0.8,0.8])

    ax.set_yticks(range(len(ylabels)))

    ax.set_yticklabels(ylabels)

    ax.set_xticks(range(len(xlabels)))

    ax.set_xticklabels(xlabels)
    ax.xaxis.tick_top()

    vmax=data[0][0]

    vmin=data[0][0]

    for i in data:

        for j in i:

            if j>vmax:

                vmax=j

            if j<vmin:

                vmin=j

    map=ax.imshow(data,interpolation='nearest',cmap=cmap,aspect='auto',vmin=vmin,vmax=vmax)

    cb=plt.colorbar(mappable=map,cax=None,ax=None,shrink=0.5)
    
    out_f = "data/H_" + str(index) +".eps"
    plt.savefig(out_f,format='eps')
    #plt.show()

for i in i_list:#[0:1]:
    F_1 = "data/srg_conv_H_" + str(i) + ".txt"
    data_F = np.loadtxt(F_1)
    # x_len = len(data_F)
    # x = np.arange(0,x_len)
    # y = x
    # plt.contour(x, y, data_F[x][y], 12, alpha = 1, cmap = plt.cm.Blues)
    # out_g = "data/temp_G_Test" + str(i) +".eps"
    # plt.savefig(out_g,format='eps')
    # print(out_g)
    xlabels=['$0p0h$','$1p1h$','$2p2h$','$3p3h$','$4p4h$']

    #ylabels=['a','b','c','d','e','f','g','h','i','j']
    
    len_data = len(data_F)
    #print(data_F)
    #data_F=np.random.rand(6,6)
    for x in range(0,len_data):
        for y in range(0,len_data):
            if(x==y and data_F[x][y]>1):
                print(x,y)
                data_F[x][y] = data_F[x][y]/10.0
            # if(data_F[x][y]>1):
            #     data_F[x][y] = math.log(data_F[x][y])
    
    print(data_F)
    print("len_data : ",len_data)
    draw_heatmap(data_F,xlabels[0:len_data],xlabels[0:len_data],i)









    
