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
    def __init__(self,sp_num,energy_gap):
        self.sp_num_p = sp_num
        self.sp_num = sp_num * 2
        self.energy_gap=energy_gap
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

    def print_state(self):
        print ("index \t p \t spin \t energy")
        i=0
        for sp_t in self.state:
            print(i,sp_t.p,sp_t.spin,sp_t.energy)
            i+=1

