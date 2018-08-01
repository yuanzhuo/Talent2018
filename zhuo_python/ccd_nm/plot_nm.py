#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 15:48:23 2018

@author: zhuo
"""

import numpy as np
import matplotlib.pyplot as plt

F_1 = "data/Nmax_1.txt"
data_F = np.loadtxt(F_1)
len_dat = len(data_F)
x = np.zeros((len_dat))
y = np.zeros((len_dat))

for i in np.arange(0,len_dat):
    x[i] = data_F[i][0]
    y[i] = data_F[i][1]/14.0

data = y
print(data)
plt.plot(x,y,'*', linewidth=3)
plt.savefig("data/temp.eps")