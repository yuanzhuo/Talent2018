#!/usr/bin/python
import numpy as np
import itertools
import sys
import matplotlib.pyplot as plt
from basis import *
from slater import *
from pairing_me import *
from hamiltonian import *
from ccsd_me import *


# sp_obits_num_p = 8
# particle_num = 4
#
# pairing_d = 1
# pairing_g = -1*0.5
#
# pairing_g = pairing_g*2
#
# sps_t=Basis_SP(sp_obits_num_p,pairing_d,particle_num)
# sps_t.build()
# sps_t.print_state()

# ==============++++++++++==============#
# basis_exc = PH_Combain(sps_t,particle_num)
# ph_type_vec=[2,4]
# basis_exc.build(ph_type_vec,1)
# # basis_exc.print_BE()
# # basis_exc.print_pair()
#
# slater = Build_Slater(sps_t.sp_num,particle_num)
# slater.build_ph(basis_exc)
# slater.print_ph()
#
# me = Pairing_ME(sps_t)
# g_vec = np.arange(-1.0,1.1,0.1)
# res_vec = np.zeros((len(g_vec),2))
# ccsd_vec = np.zeros((len(g_vec),2))
# i=0
# for g_v in g_vec:
#     pairing_g = g_v
#     me.build(1,pairing_g)
#     h = H_system(sps_t,slater,me)
#     h.build_me()
#     E_corr = h.diag()[1]
#     #E_corr = -1
#     print(E_corr)
#     #res_vec[i]=E_corr
#     res_vec[i][0]=g_v
#     res_vec[i][1]=E_corr
#     i+=1
#
# i=0
# for res in res_vec:
#     print(res[0],'\t',res[1])
#     i+=1
#
# np.savetxt('eight_sp14.txt',res_vec)

# print(res_vec[0])
# #plt.plot(res_vec[0],res_vec[1])
#


# ==============++++++++++==============#
sp_obits_num_p = 12
particle_num = 4

sps_t=Basis_SP(sp_obits_num_p,1,particle_num)
sps_t.build()
sps_t.print_state()
pairing_me = Pairing_ME(sps_t)
#
i=0
#g_vec = [-0.8]
g_vec = np.arange(-1,1,0.2)
#ccsd_limt = np.zeros((len(g_vec),3))
#ccsd_conv = np.zeros((len(g_vec),2))
ccsd_conv = np.zeros((len(g_vec),200))
#g_vec = [-1]
for g_v in g_vec:
    print("g_v : ", g_v)
    pairing_g = g_v
    pairing_me.build(1,pairing_g)
    ccsd_me = CCSD_ME(sps_t,pairing_me)
    H_bar = ccsd_me.build_Hbar()
    ccsd_me.cal_New_T2_pphh(H_bar)
    Iter_Info = ccsd_me.iter()
    Ec = ccsd_me.Ec_cal()
    loop_time = ccsd_me.loop_time
    j = 0
    for it in Iter_Info:
        ccsd_conv[i][j]=it
        j+=1
    # ccsd_limt[i][0]=g_v
    # ccsd_limt[i][1]=loop_time
    # ccsd_limt[i][2]=Ec
    i+=1
    print("== g_v == : ",g_v,"\t loop_time : ",ccsd_me.loop_time,"\t Ec : ",Ec)
#
np.savetxt("data/ccsd_conv_four_sp12.txt",ccsd_conv)
# conf = str(10)+"-"+str(1)+"-"+str(10)+"-"+str(1)
# pairing_me.print_me()
# print(conf,pairing_me.me_dic[conf])








