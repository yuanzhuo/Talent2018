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
        self.sp_num = sp_num
        self.energy_gap=energy_gap
        self.state = []
        self.state_size = 0

    def build(self):
        p=0
        index = 0
        while(p<self.sp_num):
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
        beg = 0
        end = self.fermi_index

        particle_combin = self.cal_combination(beg,end,ph_type)
        for h in hole_combin:
            for p in particle_combin:
                SE_t =State_Excit(h,p)
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
                ph_state_t = self.cal_ph_pair_4p4h()
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
        num_list = np.arange(0,particle_num)
        for i in num_list:
            sla_int_0 += 2**i
        bit_str = "0" + str(sp_num)+ "b"
        #sla_b_0 = format(sla_int_0,bit_str)
        #sla_0 = np.array([sla_b_0])
        sla_b_0 = sla_int_0<<(sp_num-particle_num)
        #sla_b_0 = bin(sla_int_0)
        #sla_b_0 = 0b1111
        print(format(sla_b_0,'0b'))
        sla_b_0 = sla_b_0 << 8
        #sla_b_0 = bin(sla_b_0)
        print(format(sla_b_0,'0b'))
        self.sla_0 = sla_int_0
        self.sla_particle_num = self.countBit(sla_int_0)
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

    def cal_sla_bin(h_index,p_index,fermi_index):
        if(self.sla_particle_num != fermi_index):
            print('Wrong happened at cal_sla_bin')
        sla_e = self.sla_0
        for h in h_index:
            self.clearBit(sla_e,h)
        for p in p_index:
            self.setBit(sla_e,h)


    def build_ph(self, ph_comb):
        size_conf_type = len(ph_comb.ph_state)
        self.sla_ph = np.zeros((size_conf_type,size_conf_type))
        for ph_state in ph_comb.ph_state:
            for pair in ph_state:
                i=0




sps_t=Basis_SP(4,1)
sps_t.build()
#sps_t.print_state()
basis_exc = PH_Combain(sps_t,4)
#test_vec = basis_exc.cal_combination(0,3,2)
ph_type_vec=[2,4]

basis_exc.build(ph_type_vec,1)
basis_exc.print_BE()
basis_exc.print_pair()

slater = Build_Slater(8,4)
slater.countBit(15)
#S_bin =








