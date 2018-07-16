#!/usr/bin/python
import numpy as np

class SP:
    'single particle state'
    #sp_num=0;
    def __init__(self,p,spin,energy,index):
        self.p=p
        self.spin=spin
        self.energy=energy
        self.index=index
        #Basis_SP.sp_num+=1

class Basis_SP:
    'A collection of Single Particle Basis'
    states = []
    def __init__(self,sp_num,energy_gap):
        self.sp_num = sp_num
        self.energy_gap=energy_gap

    def build(self):
        p=0
        index = 0
        while(p<self.sp_num):
            spin = -1
            energy = p * self.energy_gap
            sp_t=SP(p,spin,energy,index)
            Basis_SP.states.append(sp_t)
            index+=1
            spin = +1
            sp_t=SP(p,spin,energy,index)
            Basis_SP.states.append(sp_t)
            index+=1
            p+=1

    def print_states(self):
        print ("index \t p \t spin \t energy")
        i=0
        for sp_t in Basis_SP.states:
            print(i,sp_t.p,sp_t.spin,sp_t.energy)
            i+=1

class State_Excit:
    index_a=-1
    index_b=-1
    Flag=-1  # 0: 0p0h  1: 1p1h  2:2p2h ....
    def __init__(self,flag,index_a,index_b):
        self.Flag = flag
        self.index_a = index_a
        self.index_b = index_b

#class Basis_Excit:
#    'A collection particle-hole excition'




sps_t=Basis_SP(2,1)
sps_t.build();
sps_t.print_states();






