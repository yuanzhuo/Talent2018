#!/usr/bin/python
import numpy as np
import itertools
import sys
from basis import *
from slater import *
from pairing_me import *
from hamiltonian import *


sp_obits_num_p = 12
particle_num = 4

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

me = Pairing_ME(sps_t)
me.build(1,pairing_g)
#me.print_me()
h = H_system(sps_t,slater,me)
A = h.position2Bit(460,4)
print(A)
print(h.cal_me(496,496))
h.build_me()
h.print_me()
h.diag()










