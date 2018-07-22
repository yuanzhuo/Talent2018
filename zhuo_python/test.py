#!/usr/bin/python
import numpy as np
import itertools
import sys
import matplotlib.pyplot as plt
from basis import *
from slater import *
from pairing_me import *
from hamiltonian import *


sp_obits_num_p = 14
particle_num = 8

pairing_d = 1
pairing_g = -1*0.5

pairing_g = pairing_g*2

sps_t=Basis_SP(sp_obits_num_p,pairing_d)
sps_t.build()
#sps_t.print_state()
basis_exc = PH_Combain(sps_t,particle_num)
#test_vec = basis_exc.cal_combination(0,3,2)
#ph_type_vec=[2,4]
ph_type_vec=[2,4]

basis_exc.build(ph_type_vec,1)
basis_exc.print_BE()
basis_exc.print_pair()

slater = Build_Slater(sps_t.sp_num,particle_num)
slater.build_ph(basis_exc)
slater.print_ph()
#slater.countBit(15)
#h = [0]
#p = [6]
#sla_e=slater.cal_sla_bin(h,p,3)
#print(format(sla_e,'0b'))

#S_bin =

pairing_g = -1*0.5
pairing_g = pairing_g*2

g_vec = np.arange(-1,1,0.1)
res_vec = np.zeros((len(g_vec),2))
i=0
for g_v in g_vec:
    pairing_g = g_v*2
    me = Pairing_ME(sps_t)
    me.build(1,pairing_g)
    #me.print_me()
    h = H_system(sps_t,slater,me)
    #A = h.position2Bit(460,4)
    #print(A)
    #print(h.cal_me(496,496))
    h.build_me()
    #h.print_me()
    E_corr = h.diag()[1]
    print(E_corr)
    #res_vec[i]=E_corr
    res_vec[i][0]=g_v
    res_vec[i][1]=E_corr
    i+=1

i=0
for res in res_vec:
    print(res[0],'\t',res[1])
    i+=1

np.savetxt('eight_sp14.txt',res_vec)
#print(res_vec[0])
#plt.plot(res_vec[0],res_vec[1])











