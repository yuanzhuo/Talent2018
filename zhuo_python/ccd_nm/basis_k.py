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
        self.E = 0
        #Basis_SP.sp_num+=1

class Basis_SP_K:
    'A collection of Single Particle Basis'
    def __init__(self,N_max,particle_num):
        self.hbarc =  197.3269678792965;    # Mev.fm
        self.M_n = 938.918725 # Mev/c2
        self.hc2_2M = self.hbarc**2/(2.0*self.M_n)
        self.N_max = N_max
        self.particle_num = particle_num
        self.state = []
        self.state_size = 0
        #self.cal_nL(0, 0.08)

    def cal_nL(self,theta,rho):
        self.L_k = np.cbrt(self.particle_num/(1.0*rho))
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
        for k_x in self.n_k:
            for k_y in self.n_k:
                for k_z in self.n_k:
                    k2_tot = k_x*k_x + k_y*k_y + k_z*k_z
                    sp_kt = SP_K([k_x,k_y,k_z],-1,k2_tot)
                    index += 1
                    self.state.append(sp_kt)
                    sp_kt = SP_K([k_x,k_y,k_z],1,k2_tot)
                    self.state.append(sp_kt)
                    index += 1
        self.state.sort(key=self.take_k2_tot)
        self.state_size = len(self.state)
        i = 0
        while i<len(self.state):
            self.state[i].index = i
            i+=1

        self.find_fermi()
        self.build_Tk()
        print("SP fermi_index : ",self.fermi_index)
        print("state size : ",len(self.state))


    def find_fermi(self):
        self.fermi_index = self.particle_num - 1

    def build_Tk(self):

        i_list = np.arange(0,self.state_size)
        #print("self.state_size ?? : ",self.state_size,i_list)
        for i in i_list:
            kp_xyz = self.state[i].k_xyz
            k_xyz = kp_xyz
            val = self.cal_Tk(kp_xyz,k_xyz)
            #print("kinetic i ",i,val)
            self.state[i].E = val


    def cal_Tk(self,kp_xyz,k_xyz):
        # kinetic energy
        # kp_xyz[0,1,2] k_x,k_y,k_z
        # k_xyz[0,1,2]  k_x,k_y,k_z
        res = 0.0
        if(kp_xyz[0]!=k_xyz[0]):
            return res
        if(kp_xyz[1]!=k_xyz[1]):
            return res
        if(kp_xyz[2]!=k_xyz[2]):
            return res
        k_xyz_r = kp_xyz*self.k_step
        k2 = np.dot(k_xyz_r,k_xyz_r)
        res = self.hc2_2M*k2
        return res

    def print_state(self):
        print("index\t k_x\t k_y\t k_z\t S2_z\t k^2_tot\t E")
        for sp in self.state:
            print(sp.index, "\t ",round(sp.k_xyz[0],4),"\t ",round(sp.k_xyz[1],4),"\t ",round(sp.k_xyz[2],4),"\t ",round(sp.s2_z,4),"\t ",round(sp.k2_tot,4),round(sp.E,4))


class TB_k:
    'two body state'
    def __init__(self,index_1,index_2,K_tot,k12,s2_1,s2_2,S12_2):
        # K_tot [0,1,2] x,y,z
        # k12 [0,1,2] kx1-kx2,ky1-ky2,kz1-kz2
        self.index_1 = index_1
        self.index_2 = index_2
        self.s2_1 = s2_1
        self.s2_2 = s2_2
        self.K_tot = K_tot
        self.k12 = k12
        self.S12_2 = S12_2
        self.Channel = -1
        self.index = -1
        self.index_P = -1

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
        self.channel_size = len(self.P12_channel)
        self.sp_size = len(sps_k.state)
        self.nu_size = len(sps_k.state)-sps_k.fermi_index-1
        self.A_size = sps_k.fermi_index+1
        self.A_list = np.arange(0,self.A_size)
        self.nu_list = np.arange(0,self.nu_size)

    def cal_P12_channel(self):
        p12_vec = []
        print("self.sps_k.k_vec : ",self.sps_k.k_vec)
        k12_min = 2*min(self.sps_k.k_vec)
        k12_max = 2*max(self.sps_k.k_vec)
        p12_vec = np.arange(k12_min,k12_max+self.sps_k.k_step,self.sps_k.k_step)
        p12_vec = np.arange(-2*self.sps_k.N_max,2*self.sps_k.N_max+1)
        print(p12_vec)
        #sys.exit()
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
        self.print_P12_channel()

    def print_P12_channel(self):
        for p in self.P12_channel:
            print(p.index,"\t",round(p.P12_2tot,4),"\t",round(p.P12_xyz[0],4),"\t",round(p.P12_xyz[1],4),"\t",round(p.P12_xyz[2],4))

    def take_k2_tot(self,elem):
        return elem.P12_2tot

    # def cal_P12_channel_P_index(self):
    #     #ch = -1
    #     ch_list = np.arange(self.channel_size)
    #     for ch in ch_list:
    #
    #         tb_size_hh = self.state_hh[ch]
    #         tb_size_pp = self.state_pp[ch]
    #         self.P12_channel[ch].P_index_hh.reshape(tb_size_hh)
    #         self.P12_channel[ch].P_index_pp.reshape(tb_size_pp)


    def cal_P_index(self):
        ch_list = np.arange(self.channel_size)
        for ch in ch_list:
            len_hh = len(self.state_hh[ch])
            len_hh_l = np.arange(len_hh)
            for tb_i in len_hh_l:
                index_1 = self.state_hh[ch][tb_i].index_1
                index_2 = self.state_hh[ch][tb_i].index_2
                if(self.state_hh[ch][tb_i].index_P>=0):
                    continue
                for tb_j in len_hh_l:

                    if index_1 != self.state_hh[ch][tb_j].index_2:
                        continue
                    if index_2 != self.state_hh[ch][tb_j].index_1:
                        continue
                    self.state_hh[ch][tb_i].index_P = self.state_hh[ch][tb_j].index
                    self.state_hh[ch][tb_j].index_P = self.state_hh[ch][tb_i].index

                    break
            len_pp = len(self.state_pp[ch])
            len_pp_l = np.arange(len_pp)
            for tb_i in len_pp_l:
                index_1 = self.state_pp[ch][tb_i].index_1
                index_2 = self.state_pp[ch][tb_i].index_2
                if(self.state_pp[ch][tb_i].index_P>=0):
                    continue
                for tb_j in len_pp_l:
                    if index_1 != self.state_pp[ch][tb_j].index_2:
                        continue
                    if index_2 != self.state_pp[ch][tb_j].index_1:
                        continue
                    self.state_pp[ch][tb_i].index_P = self.state_pp[ch][tb_j].index
                    self.state_pp[ch][tb_j].index_P = self.state_pp[ch][tb_i].index

                    break


    def build(self):
        #print("self.A_list",len(self.A_list),"\t",self.A_list)
        #print("self.nu_list",len(self.nu_list),"\t",self.nu_list)
        # -------------------- build state_hh --------------------
        ch = -1
        for p_ch in self.P12_channel:
            #print("ch : ",ch)
            ch +=1
            p_ch_xyz = p_ch.P12_xyz
            #print(p_ch_xyz)
            tb_k_ch = []
            index = 0
            for a_i in self.A_list:
                a= self.sps_k.state[a_i]
                a_xyz = a.k_xyz
                for b_i in self.A_list:
                    b = self.sps_k.state[b_i]
                    if(a == b):
                        continue
                    b_xyz = b.k_xyz
                    ab_xyz = a_xyz + b_xyz
                    if(ab_xyz[0] != p_ch_xyz[0]):
                        continue
                    if(ab_xyz[1] != p_ch_xyz[1]):
                        continue
                    if(ab_xyz[2] != p_ch_xyz[2]):
                        continue
                    #print(a_xyz,"\t",b_xyz,"\t",ab_xyz,"\t",p_ch_xyz)
                    #print(ab_xyz,"\t",p_ch_xyz)
                    s2_1 = a.s2_z
                    s2_2 = b.s2_z
                    k12 = a_xyz - b_xyz
                    S12_2 = a.s2_z + b.s2_z
                    tb_t = TB_k(a.index,b.index,p_ch_xyz,k12,s2_1,s2_2,S12_2)
                    tb_t.index = index
                    tb_k_ch.append(tb_t)
                    #print(" ========= ")
                    index += 1
            #self.state.append(tb_k_ch)
            self.state_hh.append(tb_k_ch)
            #tb_k_ch = []
        # -------------------- build state_pp --------------------
        ch = -1
        for p_ch in self.P12_channel:
            #print("ch : ",ch)
            ch +=1
            p_ch_xyz = p_ch.P12_xyz
            #print(p_ch_xyz)
            tb_k_ch = []
            index = 0
            for a_i in self.nu_list:
                a_ii = a_i + self.sps_k.fermi_index + 1
                a= self.sps_k.state[a_ii]
                a_xyz = a.k_xyz
                for b_i in self.nu_list:
                    b_ii = b_i + self.sps_k.fermi_index + 1
                    b = self.sps_k.state[b_ii]
                    if(a == b):
                        continue
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
                    s2_1 = a.s2_z
                    s2_2 = b.s2_z
                    k12 = a_xyz - b_xyz
                    S12_2 = a.s2_z + b.s2_z
                    tb_t = TB_k(a.index,b.index,p_ch_xyz,k12,s2_1,s2_2,S12_2)
                    tb_t.index = index
                    tb_k_ch.append(tb_t)
                    index += 1
            #self.state.append(tb_k_ch)
            self.state_pp.append(tb_k_ch)
            tb_k_ch = []
        print("tb_k.size_hh = ",len(self.state_hh))
        print("tb_k.size_pp = ",len(self.state_pp))
        self.cal_P_index()
        ch = 0
        tb_tot_size_hh = 0
        tb_tot_size_pp = 0

        # for tb_t in self.state_pp[0]:
        #     print(tb_t.index,"\t i_1:",tb_t.index_1,"\t i_2:",tb_t.index_2,"\t i_P:",tb_t.index_P)

        for tb_t in self.state_hh:
            if(len(tb_t)>0):
                #print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state_hh[ch].size : ",len(tb_t))
                print("ch : ",ch,"tb_hh size : ",len(self.state_hh[ch]))
                tb_tot_size_hh += len(self.state_hh[ch])
            ch += 1
        ch = 0
        for tb_t in self.state_pp:
            if(len(tb_t)>0):
                #print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state_pp[ch].size : ",len(tb_t))
                print("ch : ",ch,"tb_pp size : ",len(self.state_pp[ch]))
                tb_tot_size_pp += len(self.state_pp[ch])
            ch += 1
        print("Tot tb_size_hh : ",tb_tot_size_hh)
        print("Tot tb_size_pp : ",tb_tot_size_pp)

    def find_channel(self,K_tot):
        for p_ch in self.P12_channel:
            p_ch_xyz = p_ch.P12_xyz
            if((K_tot[0] == p_ch_xyz[0])and(K_tot[1] == p_ch_xyz[1])and(K_tot[2] == p_ch_xyz[2])):
                return p_ch.index
        print("Wrong happened at find_channel")
        sys.exit()

    def build_2(self):
        #print("self.A_list",len(self.A_list),"\t",self.A_list)
        #print("self.nu_list",len(self.nu_list),"\t",self.nu_list)
        # -------------------- build state_hh --------------------
        len_ch = len(self.P12_channel)
        tb_l = []

        index_ch = np.zeros((len_ch))
        for a_i in self.A_list:
            a= self.sps_k.state[a_i]
            a_xyz = a.k_xyz
            for b_i in self.A_list:
                b = self.sps_k.state[b_i]
                if(a == b):
                    continue
                b_xyz = b.k_xyz
                ab_xyz = a_xyz + b_xyz

                #print(a_xyz,"\t",b_xyz,"\t",ab_xyz,"\t",p_ch_xyz)
                #print(ab_xyz,"\t",p_ch_xyz)
                s2_1 = a.s2_z
                s2_2 = b.s2_z
                k12 = a_xyz - b_xyz
                S12_2 = a.s2_z + b.s2_z
                tb_t = TB_k(a.index,b.index,ab_xyz,k12,s2_1,s2_2,S12_2)
                ch = self.find_channel(ab_xyz)
                #print("ch : ",ch)
                tb_t.Channel = ch
                tb_l.append(tb_t)
                #print(" ========= ")
                #index += 1
        #self.state.append(tb_k_ch)
        print("-------!!!! -------test_point 2 ----  1")
        tb_l.sort(key=lambda tb_l: tb_l.Channel)
        print("-------!!!! -------test_point 2 ----  2")


        ch_flag = 0

        for ch in np.arange(0,self.channel_size):
            tb_ch = []
            index = 0
            for tb in tb_l:
                if(tb.Channel == ch):
                    tb_t = tb
                    tb_t.index = index
                    tb_ch.append(tb)
                    index +=1

            self.state_hh.append(tb_ch)
        print("-------!!!! -------test_point 2 ----  3")

        # i=0
        # ch = 0
        # for tb_t in self.state_hh:
        #     if(len(tb_t)>0):
        #         print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state_hh[ch].size : ",len(tb_t))
        #         #print(i,tb_t.index_1,tb_t.index_2,tb_t.Channel)
        #     i+=1
        #     #tb_tot_size_hh += len(self.state_hh[ch])
        #     ch += 1
        # #quit()

        #tb_k_ch = []
        # -------------------- build state_pp --------------------
        tb_l = []
        for a_i in self.nu_list:
            a_ii = a_i + self.sps_k.fermi_index + 1
            a= self.sps_k.state[a_ii]
            a_xyz = a.k_xyz
            for b_i in self.nu_list:
                b_ii = b_i + self.sps_k.fermi_index + 1
                b = self.sps_k.state[b_ii]
                if(a == b):
                    continue
                b_xyz = b.k_xyz
                ab_xyz = (a_xyz + b_xyz)
                ch = self.find_channel(ab_xyz)
                if(len(self.state_hh[ch])==0):
                    continue
                #print(ab_xyz,"\t",p_ch_xyz)
                #print(a_xyz,"\t",b_xyz,"\t",ab_xyz,"\t",p_ch_xyz)
                s2_1 = a.s2_z
                s2_2 = b.s2_z
                k12 = a_xyz - b_xyz
                S12_2 = a.s2_z + b.s2_z
                tb_t = TB_k(a.index,b.index,ab_xyz,k12,s2_1,s2_2,S12_2)
                tb_t.index = index

                #print("ch : ",ch)
                tb_t.Channel = ch
                tb_l.append(tb_t)
        print("-------!!!! -------test_point 2 ----  4")

        tb_l.sort(key=lambda tb_l: tb_l.Channel)

        ch_flag = 0

        for ch in np.arange(0,self.channel_size):
            tb_ch = []
            index = 0
            for tb in tb_l:
                if(tb.Channel == ch):
                    tb_t = tb
                    tb_t.index = index
                    tb_ch.append(tb)
                    index +=1

            self.state_pp.append(tb_ch)

        print("tb_k.size_hh = ",len(self.state_hh))
        print("tb_k.size_pp = ",len(self.state_pp))
        self.cal_P_index()
        ch = 0
        tb_tot_size_hh = 0
        tb_tot_size_pp = 0

        # for tb_t in self.state_pp[0]:
        #     print(tb_t.index,"\t i_1:",tb_t.index_1,"\t i_2:",tb_t.index_2,"\t i_P:",tb_t.index_P)

        for tb_t in self.state_hh:
            if(len(tb_t)>0):
                #print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state_hh[ch].size : ",len(tb_t))
                #print("ch : ",ch,"tb_hh size : ",len(self.state_hh[ch]))
                tb_tot_size_hh += len(self.state_hh[ch])
            ch += 1
        ch = 0
        for tb_t in self.state_pp:
            if(len(tb_t)>0):
                #print("ch : ",ch,"\t",self.P12_channel[ch].P12_2tot,"\t",self.P12_channel[ch].P12_xyz,"\t state_pp[ch].size : ",len(tb_t))
                #print("ch : ",ch,"tb_pp size : ",len(self.state_pp[ch]))
                tb_tot_size_pp += len(self.state_pp[ch])
            ch += 1
        print("Tot tb_size_hh : ",tb_tot_size_hh)
        print("Tot tb_size_pp : ",tb_tot_size_pp)


    def print_ch(self,ch):
        tb_ch=self.state_hh[ch]
        print("ch : ",ch,"self.P12_channel.P12_xyz",self.P12_channel[ch].P12_xyz)
        for tb in tb_ch:
            a = tb.index_1
            b = tb.index_2
            print("------ ----- ")
            print("i : ",tb.index,"\t S = ",tb.S12_2,"\t P12: ",tb.k12/2)

            print("i : ",tb.index,"\t S = ",tb.S12_2,"\t a: ",a," s ",self.sps_k.state[a].s2_z," b: ",b,"\t sb",self.sps_k.state[b].s2_z)
            print(" \t",self.sps_k.state[a].k_xyz,"\t  ",self.sps_k.state[b].k_xyz)









