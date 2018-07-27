#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 14:22:15 2018

@author: zhuo
"""
import numpy as np
import matplotlib.pyplot as plt



F_1 = "data/four_sp4.txt"
F_2 = "data/four_sp8.txt"
F_3 = "data/four_sp12.txt"

E_1 = "data/eight_sp10.txt"
E_2 = "data/eight_sp12.txt"
E_3 = "data/eight_sp14.txt"

F = [F_1,F_2,F_3]
E = [E_1,E_2,E_3]

F_fig = ["FCI SP orbits: 4","FCI SP orbits: 8","FCI SP orbits: 12"]
E_fig = ["FCI SP orbits: 10","FCI SP orbits: 12","FCI SP orbits: 14"]

ccsd_F_1 = "data/ccsd_four_sp4.txt"
ccsd_F_2 = "data/ccsd_four_sp8.txt"
ccsd_F_3 = "data/ccsd_four_sp12.txt"
F_ccsd = [ccsd_F_1,ccsd_F_2,ccsd_F_3]

ccsd_E_1 = "data/ccsd_eight_sp10.txt"
ccsd_E_2 = "data/ccsd_eight_sp12.txt"
ccsd_E_3 = "data/ccsd_eight_sp14.txt"
E_ccsd = [ccsd_E_1,ccsd_E_2,ccsd_E_3]

F_fig_C = ["CCD SP orbits: 4","CCD SP orbits: 8","CCD SP orbits: 12"]
E_fig_C = ["CCD SP orbits: 10","CCD SP orbits: 12","CCD SP orbits: 14"]






data_four =[]
f=0
for file in F:
    data_four_t = np.loadtxt(file)
    data = np.zeros((2,len(data_four_t)))
    i=0
    while i<len(data_four_t):
        data[0][i]=data_four_t[i][0]
        data[1][i]=data_four_t[i][1]
        i+=1
    plt.plot(data[0],data[1], linewidth=3, label=F_fig[f])

    standard = np.zeros(len(data[0]))

    plt.plot(data[0],standard, '--')
    f+=1
f=0
for file in F_ccsd:
    data_four_t = np.loadtxt(file)
    data = np.zeros((2,len(data_four_t)))
    i=0
    while i<len(data_four_t):
        data[0][i]=data_four_t[i][0]
        data[1][i]=data_four_t[i][1]
        i+=1
    plt.plot(data[0],data[1],'*', linewidth=3, label=F_fig_C[f])


    f+=1
    

plt.xlim(-1.05,1.05)
plt.legend(loc='lower left')
plt.text(60, .025, r'particle num : 8')
# data_four_t = np.loadtxt(file)
# data = np.zeros((2,len(data_four_t)))
# i=0
# while i<len(data_four_t):
#     data[0][i]=data_four_t[i][0]
#     data[1][i]=data_four_t[i][1]
#     i+=1
# plt.plot(data[1], linewidth=3)

data_four_t = np.loadtxt(F_1)
data_four.append(data_four_t)
plt.savefig("data/temp.eps")



#print(data)
#plt.plot(data[1],marker="*")