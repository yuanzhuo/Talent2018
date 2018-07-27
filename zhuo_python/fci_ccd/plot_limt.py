#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:48:23 2018

@author: zhuo
"""

import numpy as np
import matplotlib.pyplot as plt



F_1 = "data/ccsd_conv_four_sp12.txt"

data_F = np.loadtxt(F_1)

g_num =len(data_F)
print(len(data_F))

data_F_t = data_F[0:10:1]

F_fig_g = ["g = -1","g = -0.8","g = -0.6","g = -0.4","g = -0.2","g = -0.0","g = 0.2","g = 0.4","g = 0.6","g = 0.8","g = 1"]
f = 0
for data_t in data_F_t:
    index = np.where(data_t<-0.000001)[0]
    #print(index)
    #print(len(index))

    data = np.zeros((len(index)))
    j=0
    for i in index:
        data[j] = data_t[i]
        j+=1
    if(len(data) == 0):
        continue
    data_final = data[len(data)-1]

    # j=0
    # for i in index:
    #     data[j] = abs(data[j] - data_final)
    #     j+=1

    #print(data)
    #plt.plot(np.log2(data))
    print("f : ",f)
    plt.plot(np.log2(data), linewidth=3,label=F_fig_g[f])
    #plt.plot(data, linewidth=3,label=F_fig_g[f])

    f +=1
    #plt.plot(np.log2(-1*data), linewidth=3)
plt.xlim(1,200)
plt.legend(loc='uper right')
plt.text(65, -2.5, r'Particle Num : 4 SP: 12')
plt.ylabel("ln(E_diff)")
#plt.ylabel("E_corr")
plt.xlabel("loop number")
plt.savefig("data/temp.eps",format='eps')
