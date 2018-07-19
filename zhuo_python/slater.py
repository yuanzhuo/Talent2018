#!/usr/bin/python
import numpy as np
import itertools
import sys
from basis import *

class State_Excit:
    # index_p & _h : particle & hole single particle index
    def __init__(self,tup_h,tup_p):
        self.tup_h = tup_h
        self.tup_p = tup_p

class PH_Combain:
    'A collection particle-hole excition'
    # 0: 0p0h  1: 1p1h  2:2p2h ....
    def __init__(self, basis_sp, particle_num):
        if ( particle_num > len(basis_sp.state) ):
            print('particle_num > len(basis_sp.state)')
            sys.exit()
        self.sp = basis_sp
        self.particle_num = particle_num
        self.state_ph = []
        #self.sp.print_state()
        self.fermi_index = self.find_fermi()

    def find_fermi(self):
        fermi_index = -1;
        print("###",self.particle_num)
        if((self.particle_num% 2) == 0):
            fermi_index=self.particle_num-1
        else:
            print('odd particle can`t find fermi surface')
            fermi_index=self.particle_num-1
            sys.exit()
        return fermi_index

    def cal_combination(self, beg, end, num):
        if ( (beg>end) or (end-beg+1 )< num  ):
            print('Wrong happened at cal_combination')
            sys.exit()
        list_t = np.arange(beg,end+1)
        combin_tup = itertools.combinations(list_t, num)
        return combin_tup


    def cal_ph(self, ph_type):
        # ph_state [i](tup_h,tup_p)
        ph_state_t = []
        beg = 0
        end = self.fermi_index

        hole_combin = self.cal_combination(beg,end,ph_type)
        beg = self.fermi_index+1
        end = len(self.sp.state)-1
        particle_combin = self.cal_combination(beg,end,ph_type)
        for h in hole_combin:
            for p in particle_combin:
                h_t = (h)
                p_t = (p)
                SE_t =State_Excit(h_t,p_t)
                ph_state_t.append(SE_t)

        return ph_state_t

    def cal_ph_pair_2p2h(self):
        ph_type = 2
        ph_state_pair = []
        beg_x = 0
        end_x = len(self.sp.state)-1
        beg = self.sp.state[beg_x].p
        end = self.sp.state[end_x].p
        combin_pair = self.cal_combination(beg, end, ph_type  )
        fermi_index_p = self.sp.state[self.fermi_index].p

        for pair in combin_pair:
            size_pair = len(pair)
            h = pair[0]
            p = pair[1]
            #print(h,p,fermi_index_p)
            if(h<=fermi_index_p and p>fermi_index_p):

                h_tup = (2*h,2*h+1)
                p_tup = (2*p,2*p+1)
                SE_t =State_Excit(h_tup,p_tup)
                ph_state_pair.append(SE_t)
        return ph_state_pair

    def cal_ph_pair_4p4h(self):

        ph_state_pair = []
        beg_x = 0
        end_x = len(self.sp.state)-1
        beg = self.sp.state[beg_x].p
        end = self.sp.state[end_x].p
        combin_pair = self.cal_combination( beg, end, 4 )
        fermi_index_p = self.sp.state[self.fermi_index].p

        for pair in combin_pair:
            size_pair = len(pair)
            h_list = (pair[0],pair[1])
            p_list = (pair[2],pair[3])
            h_list_max = max(h_list)
            h_list_min = min(p_list)
            #print(h,p,fermi_index_p)
            if(h_list_max<=fermi_index_p and h_list_min>fermi_index_p):
                h_tup = (2*pair[0],2*pair[0]+1,2*pair[1],2*pair[1]+1)
                p_tup = (2*pair[2],2*pair[2]+1,2*pair[3],2*pair[3]+1)
                SE_t =State_Excit(h_tup,p_tup)
                ph_state_pair.append(SE_t)
        return ph_state_pair

    def build(self,ph_type_vec,pair_flag):
        # ph_state [ph_type][i](tup_h,tup_p)
        # pair_flag = 1 noly build pair_states
        self.ph_state = []
        if(pair_flag != 1):
            for ph_type in ph_type_vec:
                ph_state_t = self.cal_ph(ph_type)
                self.ph_state.append(ph_state_t)
        else:
            ph_cal_dic = {0:self.cal_ph(0),2:self.cal_ph_pair_2p2h(),4:self.cal_ph_pair_4p4h()}
            for ph_type in ph_type_vec:
                ph_state_t = ph_cal_dic[ph_type]
                self.ph_state.append(ph_state_t)

    def print_BE(self):
        i = 0
        ii = self.sp.state_size-1
        print('self.fermi_index',self.fermi_index)
        for sp_t in reversed(self.sp.state):
            print(ii,'\t',sp_t.p,'\t',sp_t.spin,'\t',sp_t.energy)
            if(ii == self.fermi_index+1):
                print('---- fermi -- surface ---- ')
            i+=1
            ii-=1

    def print_pair(self):
        for i in self.ph_state:
            index = 0
            for j in i:
                h = j.tup_h
                p = j.tup_p
                print(index,"\t h :",h,"\t p :",p)
                index += 1

class Build_Slater:
    'build binary slater determination'
    def __init__(self,sp_num, particle_num):
        sla_int_0 = 0
        num_list = np.arange(0,particle_num+1)
        for i in num_list:
            sla_int_0 += 2**i
        bit_str = "0" + str(sp_num)+ "b"
        #sla_b_0 = format(sla_int_0,bit_str)
        #sla_0 = np.array([sla_b_0])
        sla_b_0 = sla_int_0<<(sp_num-particle_num)
        #sla_b_0 = bin(sla_int_0)
        #sla_b_0 = 0b1111
        print(format(sla_b_0,'0b'))
        #sla_b_0 = sla_b_0 << 8
        #sla_b_0 = bin(sla_b_0)
        #print(format(sla_b_0,'0b'))
        self.sla_0 = sla_b_0
        self.sla_particle_num = self.countBit(sla_int_0)-1
        self.sp_num = sp_num
        #print(format(self.sla_0,'0b'))

    def setBit(self, int_type, offset):
        mask = 1 << offset
        return(int_type | mask)

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

    def cal_sla_bin(self,h_index,p_index,fermi_index):
        if(self.sla_particle_num-1 != fermi_index):
            print('Wrong happened at cal_sla_bin fermi_index', fermi_index, 'self.sla_particle_num', self.sla_particle_num-1)
        sla_e = self.sla_0
        #print(format(sla_e,'0b'),'##')
        for h in h_index:
            sla_e = self.clearBit(sla_e,self.sp_num-h-1)
            #print('set_hole\t',h,format(sla_e,'0b'),sla_e)
        for p in p_index:
            sla_e = self.setBit(sla_e,self.sp_num-p-1)
            #print('set_particle\t',p,format(sla_e,'0b'),sla_e)

        return sla_e


    def build_ph(self, ph_comb):
        #self.sla_ph[ph_type][i]
        size_conf_type = len(ph_comb.ph_state)
        sla_ph = []
        sla_ph_block = []
        sla_0_a = [self.sla_0]
        sla_ph_block.append(sla_0_a)
        sla_ph.append(self.sla_0)
        sla_ph_num = 1
        for ph_state in ph_comb.ph_state:
            size_t = len(ph_state)
            print('ph_size : ', size_t)
            sla_0_a = []
            i=0
            for pair in ph_state:
                h = pair.tup_h
                p = pair.tup_p
                print(p,h)
                sla = self.cal_sla_bin(h,p,ph_comb.fermi_index)
                #print(sla)
                print(format(sla,'0b'),sla)
                sla_0_a.append(sla)
                sla_ph.append(sla)
                sla_ph_num += 1
                i+=1
                #print(format(sla_0_a[i],'0b'))
            sla_ph_block.append(sla_0_a)
        self.sla_ph_block = sla_ph_block
        self.sla_ph = sla_ph
        self.sla_num = sla_ph_num


    def print_ph(self):
        #print the ph_type : int & binary form
        i=0
        print('sla num : ', self.sla_num)
        for sla_ph in self.sla_ph_block:
            print("ph_type : ",i)
            i+=1
            index = 0
            #len_sla = len(sla_ph)
            if(i==1):
                print("\t i : ", index,"\t",sla_ph[0],"\t",format(sla_ph[0],'0b'))
            else:
                for sla in sla_ph:
                    print(sla)
                    print("\t i : ", index,"\t",sla,"\t",format(sla,'0b'))
                    index += 1

