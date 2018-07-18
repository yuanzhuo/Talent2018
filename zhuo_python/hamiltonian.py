#!/usr/bin/python
import numpy as np
import itertools
import sys
from slater import *
from pairing_me import *
from basis import *

class H_system:
    def __init__(self,sp,slater,me):
        self.sp_state = sp.state
        self.sla_0 = slater.sla_0
        self.sla_ph = slater.sla_ph
        self.me_dic = me.me_dic
        self.d = me.d
        self.g = me.g
        self.particle_num = slater.sla_particle_num
        self.h_tot = np.zeros((slater.sla_num,slater.sla_num))

    def clearBit(self, int_type, offset):
        mask = ~(1 << offset)
        return(int_type & mask)
    def countBit(self,a):
        '整数的二进制表达里有多少个1，复杂度仅为1的个数'
        num = 0
        while a != 0:
            #print(bin(a),bin(a-1))
            a = a & (a - 1)  # 就是这个操作，需要掌握，它的本质含义是抹去了0不考虑
            num += 1
        return num

    def position2Bit(self,a,num):
        '找到前 num 个1的位置'
        i = 0
        j = 0
        pos = []
        #print(format(a,'0b'),a)
        while(a!=1 and i < 100 and j<num):
            #print(format(a,'0b'),a)
            if((a&1)==1):
                pos.append(i)
                j+=1
            a=a>>1
            i+=1
        return pos


    def cal_me(self,sla_f,sla_i):
        'interaction part and with e_i single particle energy'
        #print(format(sla_f,'0b'),sla_f)
        sla_i_t = self.clearBit(sla_i,len(self.sp_state))
        #print(format(sla_i_t,'0b'),sla_i_t,format(sla_i,'0b'))
        diff = sla_f^sla_i_t
        #print(format(diff,'0b'),diff)
        num = self.countBit(diff)
        if(num > 5):
            return 0.0
        else:
            pos = self.position2Bit(diff,4)

        len_pos = len(pos)
        #print('pos',pos,'len_pos',len_pos)
        if(len_pos == 0):
            pos2 = self.position2Bit(sla_i,self.particle_num)
            pos3 =[]
            for i in pos2:
                pos3.append(len(self.sp_state)-i-1)

            e_i = 0
            print('pos3',pos3)
            for index in pos3:
                e_i += self.sp_state[index].p * self.d
            print('\t e_i ',e_i)
            return e_i-1*self.d
        elif(len_pos == 4):
            conf = str(pos[0]) + str(pos[1]) +str(pos[2]) +str(pos[3])
            #print('conf', conf)
            if conf in self.me_dic:
                return self.me_dic[conf]
        else:
             return 0.0

    def build_me(self):
        i = -1
        for sl_f in self.sla_ph:
            i+=1
            j = -1
            for sl_i in self.sla_ph:
                j+=1
                #print(sl_f,sl_i,'#')
                #print('\t ~~~~~~~~~ ',i,j,format(sl_f,'0b'),format(sl_i,'0b'))
                self.h_tot[i][j] = self.cal_me(sl_f,sl_i)

    def diag(self):
        a,b=np.linalg.eig(self.h_tot)
        print(a)
        print(b)

    def print_me(self):
        for i in self.h_tot:
            print(i)






