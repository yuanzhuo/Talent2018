#!/usr/bin/python
import numpy as np
import itertools
import sys

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
    def __init__(self,sp_num,energy_gap,particle_num):
        self.sp_num_p = sp_num
        self.sp_num = sp_num * 2
        self.energy_gap=energy_gap
        self.particle_num = particle_num
        self.state = []
        self.state_size = 0

    def build(self):
        p=0
        index = 0
        while(p<self.sp_num_p):
            spin = -1
            energy = p * self.energy_gap
            sp_t=SP(p,spin,energy,index)
            self.state.append(sp_t)
            index+=1
            spin = +1
            sp_t=SP(p,spin,energy,index)
            self.state.append(sp_t)
            index+=1
            p+=1
        self.state_size = len(self.state)
        self.fermi_index = self.find_fermi()
        print("SP fermi_index : ",self.fermi_index)

    def find_fermi(self):
        fermi_index = -1;
        #print("###",self.particle_num)
        if((self.particle_num% 2) == 0):
            fermi_index=self.particle_num-1
        else:
            print('odd particle can`t find fermi surface')
            fermi_index=self.particle_num-1
            sys.exit()
        return fermi_index

    def print_state(self):
        print ("index \t p \t spin \t energy")
        i=0
        for sp_t in self.state:
            print(i,"\t ",sp_t.p,"\t ",sp_t.spin,"\t ",sp_t.energy)
            i+=1

