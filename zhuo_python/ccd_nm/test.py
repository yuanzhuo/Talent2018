#!/usr/bin/python
import numpy as np
import itertools
import sys
import matplotlib.pyplot as plt
from basis_k import *
from minnesota_me import *
from ccsd_nm import *

#from ccsd_me import *



# particle_num = 14
#
# N_max = 1
# sps_k=Basis_SP_K(N_max,particle_num)
# sps_k.cal_nL(0,0.08)
# sps_k.build()
#
# tb_k = Basis_TB_k(sps_k)
# tb_k.build()
#tb_k.print_ch(0)
#quit()
# minn_me = Minn_ME(sps_k,tb_k)
# a_xyz = np.array([0,0,0])
# b_xyz = np.array([-1.1233123,0,0])
# p_ch_xyz = np.array([0,0,0])
#
# s2_1 = -1
# s2_2 = 1
# k12_ab = a_xyz - b_xyz
# S12_2 = s2_1 + s2_2
# tb_f = TB_k(0,0,p_ch_xyz,k12_ab,s2_1,s2_2,S12_2)
#
# c_xyz = np.array([-1.1233123,0,0])
# d_xyz = np.array([0,0,0])
#
# k12_cd = c_xyz - d_xyz
#
# tb_i = TB_k(1,0,p_ch_xyz,k12_cd,s2_1,s2_2,S12_2)
#
# val_0 = minn_me.cal_V_sub(0,tb_f,tb_i)
# val_1 = minn_me.cal_V_sub(1,tb_f,tb_i)
# print(a_xyz,b_xyz,c_xyz,d_xyz)
# print(val_0,val_1," \t ",val_0+val_1)
#print("pre factor : ",minn_me.pre_factor)
#print("Kp_RST * 4 : ",minn_me.Kp_RST*4)

#quit()

# ==============++++++++++==============#

N_max = 1
particle_num = 14


rho_list = np.arange(0.025,0.4,0.05)
E_tot_vec = np.zeros((len(rho_list),3))
ccsd_limt = np.zeros((len(rho_list),3))
ccsd_conv = np.zeros((len(rho_list),200))
i=0
#rho_list = [0.08]

for rho in rho_list:
    sps_k=Basis_SP_K(N_max,particle_num)
    sps_k.cal_nL(0,rho)
    sps_k.build()
    #sps_k.print_state()
    print("-------!!!! -------test_point 1 ")
    #quit()
    tb_k = Basis_TB_k(sps_k)
    print("-------!!!! -------test_point 2 ")
    tb_k.build_2()
    print("-------!!!! -------test_point 3 ")
    #quit()
    #tb_k.print_ch(38)

    minn_me = Minn_ME(sps_k,tb_k)
    print("-------!!!! -------test_point 4 ")
    ccsd_nm = CCSD_NM_ME(sps_k,tb_k,minn_me)
    print("-------!!!! -------test_point 5 ")
    #ccsd_nm.build_V_2()
    print("-------!!!! -------test_point 6 ")
    #quit()
    Iter_Info = ccsd_nm.iter()
    E_tot = ccsd_nm.Etot_cal()
    E_corr = ccsd_nm.Ec_cal()
    #ccsd_nm.E_HF_test()
    j=0
    for it in Iter_Info:
         ccsd_conv[i][j]=it
         j+=1
    #quit()
    #print("Pre-Fac",minn_me.pre_factor)

    print("** rho :",rho,"\t *** E_c : ",ccsd_nm.E_corr," *** E_ref/14 : ",ccsd_nm.E_ref/14," *** E_tot : ",ccsd_nm.E_tot)
    E_tot_vec[i][0] = rho
    E_tot_vec[i][2] = ccsd_nm.E_corr
    E_tot_vec[i][1] = ccsd_nm.E_tot

    i+=1
#np.savetxt('data/Nmax_2.txt',E_tot_vec)
np.savetxt("data/ccsd_conv_nm.txt",ccsd_conv)

# ==============++++++++++==============#





