#!/usr/bin/python
import numpy as np
import itertools
import sys
from slater import *
from pairing_me import *
from basis import *

class H_system:
    def __init__(self,sp,slater,v1b,v2b):
        self.sp_state = sp.state
        self.sla_0 = slater.sla_0
        self.sla_ph = slater.sla_ph
        #self.me_dic = me.me_dic
        self.v1b = v1b
        self.v2b = v2b
        #self.d = me.d
        #self.g = me.g
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

    def cal_me_v12b(self,sla_f,sla_i):

        'interaction part and with e_i single particle energy'
        #print(format(sla_f,'0b'),sla_f)
        sla_i_t = self.clearBit(sla_i,len(self.sp_state))
        #print(format(sla_i_t,'0b'),sla_i_t,format(sla_i,'0b'))
        diff = sla_f^sla_i_t
        same = sla_f&sla_i_t
        diff_f = sla_f^same
        diff_i = sla_i^same

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
            #print('pos3',pos3)
            e_i = 0
            e_ij = 0
            for index_i in pos3:
                e_i += self.v1b[index_i][index_i]
                for index_j in pos3:
                    e_ij += 0.5*self.v2b[index_i][index_j][index_i][index_j]

            res = e_i + e_ij
            #print('\t e_i ',e_i)
            return res
        elif(len_pos == 2):
            pos_f = self.position2Bit(diff_f,1)
            pos_i = self.position2Bit(diff_i,1)
            pos2 = self.position2Bit(sla_i,self.particle_num-2)
            pos3 =[]
            for i in pos2:
                pos3.append(len(self.sp_state)-i-1)
            e_i = self.v1b[pos_f[0]][pos_i[0]]
            e_ij = 0
            for index_i in pos3:
                e_ij += self.v2b[pos_f[0]][index_i][pos_i[0]][index_i]

            res = e_i + e_ij
            #print('\t e_i ',e_i)
            return res
        elif(len_pos == 4):
            pos_f = self.position2Bit(diff_f,2)
            pos_i = self.position2Bit(diff_i,2)
            index_f = pos_f[0]
            index_fp = pos_f[1]
            index_i = pos_i[0]
            index_ip = pos_i[1]

            res = self.v2b[index_f][index_fp][index_i][index_ip]
            return res
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
                self.h_tot[i][j] = self.cal_me_v12b(sl_f,sl_i)

    def diag(self):
        a,b=np.linalg.eig(self.h_tot)
        #print(min(a), min(a)-self.h_tot[0][0])
        #print(b)
        return min(a), min(a)-self.h_tot[0][0]

    def print_me(self):
        for i in self.h_tot:
            print(np.round(i,2))






