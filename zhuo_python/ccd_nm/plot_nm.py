#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:48:23 2018

@author: zhuo
"""

import numpy as np
import matplotlib.pyplot as plt

F_1 = "data/Nmax_1.txt"
data_F_1 = np.loadtxt(F_1)
len_dat = len(data_F_1)
x_1 = np.zeros((len_dat))
y_1 = np.zeros((len_dat))
z_1 = np.zeros((len_dat))

for i in np.arange(0,len_dat):
    x_1[i] = data_F_1[i][0]
    y_1[i] = data_F_1[i][1]/14.0
    z_1[i] = data_F_1[i][2]/14.0

F_2 = "data/Nmax_2.txt"
data_F_2 = np.loadtxt(F_2)
len_dat = len(data_F_2)
x_2 = np.zeros((len_dat))
y_2 = np.zeros((len_dat))
z_2 = np.zeros((len_dat))

for i in np.arange(0,len_dat):
    x_2[i] = data_F_2[i][0]
    y_2[i] = data_F_2[i][1]/14.0
    z_2[i] = data_F_2[i][2]/14.0

F_3 = "data/Nmax_3.txt"
data_F_3 = np.loadtxt(F_3)
len_dat = len(data_F_3)
x_3 = np.zeros((len_dat))
y_3 = np.zeros((len_dat))
z_3 = np.zeros((len_dat))

for i in np.arange(0,len_dat):
    x_3[i] = data_F_3[i][0]
    y_3[i] = data_F_3[i][1]/14.0
    z_3[i] = data_F_3[i][2]/14.0


lab_fig = ["Nmax = 1","Nmax = 2","Nmax = 3"]
plt.plot(x_1,z_1, linewidth=3,label=lab_fig[0])
plt.plot(x_2,z_2, linewidth=3,label=lab_fig[1])
plt.plot(x_3,z_3,"H", lw=8,label=lab_fig[2])

plt.ylabel("E_corr/A")
plt.xlabel("density")
plt.legend(loc='upper right')

plt.savefig("data/temp.eps")

