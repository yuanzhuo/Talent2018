#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 14:22:15 2018

@author: zhuo
"""
import numpy as np
import matplotlib.pyplot as plt



F_1 = "four_sp4.txt"
F_2 = "four_sp8.txt"
F_3 = "four_sp12.txt"

E_1 = "eight_sp10.txt"
E_2 = "eight_sp12.txt"
E_3 = "eight_sp14.txt"

F = [F_1,F_2,F_3]
E = [E_1,E_2,E_3]

data_four =[]
for file in E:
    data_four_t = np.loadtxt(file)
    data = np.zeros((2,len(data_four_t)))
    i=0
    while i<len(data_four_t):
        data[0][i]=data_four_t[i][0]
        data[1][i]=data_four_t[i][1]
        i+=1
    plt.plot(data[1], linewidth=3)

data_four_t = np.loadtxt(file)
data = np.zeros((2,len(data_four_t)))
i=0
while i<len(data_four_t):
    data[0][i]=data_four_t[i][0]
    data[1][i]=data_four_t[i][1]
    i+=1
plt.plot(data[1], linewidth=3)

data_four_t = np.loadtxt(F_1)
data_four.append(data_four_t)
plt.savefig("temp.eps")


    
#print(data)
#plt.plot(data[1],marker="*")