#!/usr/bin/python
import numpy as np
import itertools
import sys
from basis import *

class Pairing_ME:
    def __init__(self,sp):
        i=0
        self.sp_state = sp.state
        self.me_dic = {}

    def build(self,d,g):
        self.d = d
        self.g = g
        print(len(self.sp_state),10)
        for sp_a in self.sp_state:
            for sp_b in self.sp_state:
                if(sp_a.p != sp_b.p or sp_a.index == sp_b.index ):
                    continue
                ab_p = sp_a.p
                for sp_c in self.sp_state:
                    for sp_d in self.sp_state:
                        if(sp_c.p != sp_d.p or sp_c.index  == sp_d.index ):
                            continue
                        cd_p = sp_c.p

                        #print(sp_c.p,sp_d.p)
                        conf = str(sp_a.index) +"-" + str(sp_b.index) +"-"+ str(sp_c.index) +"-"+str(sp_d.index)

                        #conf = [sp_a.p, sp_b.p, sp_c.p, sp_d.p]
                        #if(ap_p == cd_p):
                            #self.me_dic[conf] = -1*g
                        #else:
                        phase = 1
                        if(sp_a.index%2 != sp_c.index%2):
                            phase = -1
                        if(ab_p == cd_p):
                            self.me_dic[conf] = -0.5*phase*g
                        else:
                            self.me_dic[conf] = -0.5*phase*g

                        #print(conf,-1)

    def print_me(self):
        print(self.me_dic)





