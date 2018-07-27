#!/usr/bin/python
import numpy as np
import itertools
import sys

class SP_K:
    'single particle state'
    #sp_num=0;
    def __init__(self,k_xyz,s2_z,k2_tot):
        self.k_xyz=np.array(k_xyz)
        self.k2_tot = k2_tot
        self.s2_z = s2_z
        self.index=-1
        #Basis_SP.sp_num+=1

class Basis_SP_K:
    'A collection of Single Particle Basis'
    def __init__(self,N_max,particle_num):
        self.N_max = N_max
        self.particle_num = particle_num
        self.state = []
        self.state_size = 0
        self.cal_nL(0, 0.08)

    def cal_nL(self,theta,rho):
        self.L_k = self.particle_num/(1.0*rho)
        self.n_k = np.arange(-1*self.N_max,self.N_max+1)
        self.k_step = (2*np.pi/self.L_k)
        self.k_vec = self.n_k * self.k_step

    def take_k2_tot(self,elem):
        return elem.k2_tot

    def build(self):
        print("self.n_k : ",self.n_k)
        print("L_k = ",self.L_k,"\t k_step = ",self.k_step)
        print("self.k_vec : ", self.k_vec)
        index = 0
        #sz2_list = np.arange(-1,1,2)
        for k_x in self.k_vec:
            for k_y in self.k_vec:
                for k_z in self.k_vec:
                    k2_tot = k_x*k_x + k_y*k_y + k_z*k_z
                    sp_kt = SP_K([k_x,k_y,k_z],-1,k2_tot)
                    index += 1
                    self.state.append(sp_kt)
                    sp_kt = SP_K([k_x,k_y,k_z],1,k2_tot)
                    self.state.append(sp_kt)
                    index += 1
        self.state.sort(key=self.take_k2_tot)
        i = 0
        while i<len(self.state):
            self.state[i].index = i
            i+=1

        self.find_fermi()
        print("SP fermi_index : ",self.fermi_index)
        print("state size : ",len(self.state))


    def find_fermi(self):
        self.fermi_index = self.particle_num - 1


    def print_state(self):
        print("index\t k_x\t k_y\t k_z\t S2_z\t k^2_tot")
        for sp in self.state:
            print(sp.index, "\t ",round(sp.k_xyz[0],4),"\t ",round(sp.k_xyz[1],4),"\t ",round(sp.k_xyz[2],4),"\t ",round(sp.s2_z,4),"\t ",round(sp.k2_tot,4))


class TB_k:
    'two body state'
    def __init__(self,index_1,index_2,K_tot,k12,S12_2):
        # K_tot [0,1,2] x,y,z
        # k12 [0,1,2] kx1-kx2,ky1-ky2,kz1-kz2
        self.index_1 = index_1
        self.index_2 = index_2
        self.K_tot = K_tot
        self.k12 = k12
        self.S12_2 = S12_2

class TB_P_Channel:
    'channel for two body state'
    def __init__(self,p12_xyz,p12_tot_2):
        self.P12_xyz = np.array(p12_xyz)
        self.P12_2tot = p12_tot_2
        self.index = -1

class Basis_TB_k:
    'two body bais'
    def __init__(self,sps_k):
        self.sps_k = sps_k
        self.P12_channel = []
        self.state_hh = []
        self.state_pp = []
        self.cal_P12_channel()
        self.sp_size = len(sps_k.state)
        self.nu_size = len(sps_k.state)-sps_k.fermi_index-1
        self.A_size = sps_k.fermi_index+1
        self.A_list = np.arange(0,self.A_size)
        self.nu_list = np.arange(0,self.nu_size)

    def cal_P12_channel(self):
        p12_vec = []
        for k1 in self.sps_k.k_vec:
            for k2 in self.sps_k.k_vec:
                #k12 = k1 - k2
                #print(k12)
                p12_vec.append(k1+k2)
        #print(p12_vec)

        for p12_x in p12_vec:
            for p12_y in p12_vec:
                for p12_z in p12_vec:
                    p12_tot_2 = p12_x*p12_x + p12_y*p12_y + p12_z*p12_z
                    P=TB_P_Channel([p12_x,p12_y,p12_z],p12_tot_2)
                    self.P12_channel.append(P)
        self.P12_channel.sort(key=self.take_k2_tot)
        i = 0
        while i<len(self.P12_channel):
            self.P12_channel[i].index = i
            i+=1
        print("P_channel_size : ",len(self.P12_channel))
        #self.print_P12_channel()

    def print_P12_channel(self):
        for p in self.P12_channel:
            print(p.index,"\t",round(p.P12_2tot,4),"\t",round(p.P12_xyz[0],4),"\t",round(p.P12_xyz[1],4),"\t",round(p.P12_xyz[2],4))

    def take_k2_tot(self,elem):
        return elem.P12_2tot

    def build(self):
        print("self.A_list",len(self.A_list),"\t",self.A_list)
        print("self.nu_list",len(self.nu_list),"\t",self.nu_list)
        # -------------------- build state_hh --------------------
        ch = 0
        for p_ch in self.P12_channel:
            #print("ch : ",ch)
            ch +=1
            p_ch_xyz = p_ch.P12_xyz
            #print(p_ch_xyz)
            tb_k_ch = []
            for a_i in self.A_list:
                a= self.sps_k.state[a_i]
                a_xyz = a.k_xyz
                for b_i in self.A_list:
                    b = self.sps_k.state[b_i]
                    b_xyz = b.k_xyz
                    ab_xyz = (a_xyz + b_xyz)
                    #print(ab_xyz,"\t",p_ch_xyz)

                    if(ab_xyz[0] != p_ch_xyz[0]):
                        continue
                    if(ab_xyz[1] != p_ch_xyz[1]):
                        continue
                    if(ab_xyz[2] != p_ch_xyz[2]):
                        continue
                    #print(a_xyz,"\t",b_xyz,"\t",ab_xyz,"\t",p_ch_xyz)
                    k12 = a_xyz - b_xyz
                    S12_2 = a.s2_z - b.s2_z
                    tb_t = TB_k(a.index,b.index,p_ch_xyz,k12,S12_2)
                    tb_k_ch.append(tb_t)
            #self.state.append(tb_k_ch)
            self.state_hh.append(tb_k_ch)
        # -------------------- build state_pp --------------------
        ch = 0
        for p_ch in self.P12_channel:
            #print("ch : ",ch)
            ch +=1
            p_ch_xyz = p_ch.P12_xyz
            #print(p_ch_xyz)
            tb_k_ch = []
            for a_i in self.nu_list:
                a= self.sps_k.state[a_i]
                a_xyz = a.k_xyz
                for b_i in self.nu_list:
                    b = self.sps_k.state[b_i]
                    b_xyz = b.k_xyz
                    ab_xyz = (a_xyz + b_xyz)
                    #print(ab_xyz,"\t",p_ch_xyz)

                    if(ab_xyz[0] != p_ch_xyz[0]):
                        continue
                    if(ab_xyz[1] != p_ch_xyz[1]):
                        continue
                    if(ab_xyz[2] != p_ch_xyz[2]):
                        continue
                    #print(a_xyz,"\t",b_xyz,"\t",ab_xyz,"\t",p_ch_xyz)
                    k12 = a_xyz - b_xyz
                    S12_2 = a.s2_z - b.s2_z
                    tb_t = TB_k(a.index,b.index,p_ch_xyz,k12,S12_2)
                    tb_k_ch.append(tb_t)
            #self.state.append(tb_k_ch)
            self.state_pp.append(tb_k_ch)
        print("tb_k.size_hh = ",len(self.state_hh))
        print("tb_k.size_pp = ",len(self.state_pp))
        # ch = 0
        # for tb_t in self.state_hh:
        #     if(len(tb_t)>0):
        #         print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state[ch].size : ",len(tb_t))
        #     ch += 1
        # ch = 0
        # for tb_t in self.state_pp:
        #     if(len(tb_t)>0):
        #         print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state[ch].size : ",len(tb_t))
        #     ch += 1








