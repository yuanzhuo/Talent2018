#!/usr/bin/python

from cython.parallel import prange, parallel, threadid
from basis_k import *
import numpy as np
import sys


class CCSD_NM_ME:
    def __init__(self,sps_k,tb_k,Minn_me):
        self.sps_k = sps_k
        self.tb_k = tb_k
        self.Minn_me = Minn_me

        self.sp_size = len(sps_k.state)
        self.channel_size = tb_k.channel_size
        self.ch_list = np.arange(0,self.channel_size)
        # self.tb_A_size = len(tb_k.state_pp)
        # self.tb_nu_size = len(tb_k.state_pp)

        self.sp_A_size = sps_k.fermi_index+1
        self.sp_nu_size = len(sps_k.state)-sps_k.fermi_index-1


        self.T_pphh = []
        self.V_pphh = []
        self.V_hhhh = [] #np.zeros((0))
        self.V_pppp = []
        self.f_me_hh = np.zeros((self.sp_A_size,self.sp_A_size))
        self.f_me_pp = np.zeros((self.sp_nu_size,self.sp_nu_size))
        self.f_dom_abij = np.zeros((self.sp_nu_size,self.sp_nu_size,self.sp_A_size,self.sp_A_size))
        self.H_bar = np.zeros((self.channel_size,0,0))

        self.sp_A_list = np.arange(0,self.sp_A_size)
        self.sp_nu_list = np.arange(0,self.sp_nu_size)
        #print("sp_A_list : ",self.sp_A_list)
        #print("sp_nu_list : ",self.sp_nu_list)


        print("-------!!!! -------test_point 4 -- 1")
        self.build_f_pq()
        print("-------!!!! -------test_point 4 -- 2")
        self.build_V_2()
        print("-------!!!! -------test_point 4 -- 3")
        self.build_T0()


    def cal_f_pq(self,index_1,index_2):
        val = 0
        e_pq = 0.0
        p = self.sps_k.state[index_1]
        q = self.sps_k.state[index_2]
        p_xyz = p.k_xyz
        q_xyz = q.k_xyz


        if(index_1 == index_2):
            #e_pq = self.Minn_me.cal_Tk(p_xyz,q_xyz)
            e_pq = self.sps_k.state[index_2].E
            #return e_pq
            #print("index_1",index_1,"E : ",self.sp.state[index_1].energy)
        v_pq = 0.0
        for index_i in self.sp_A_list:
            i = self.sps_k.state[index_i]
            s2_p = p.s2_z
            s2_q = q.s2_z
            s2_i = i.s2_z

            k12_pi = p_xyz - i.k_xyz
            k12_qi = q_xyz - i.k_xyz

            Ktot_pi = p_xyz + i.k_xyz
            Ktot_qi = q_xyz + i.k_xyz

            S12_2_pi = s2_p + s2_i
            S12_2_qi = s2_q + s2_i

            tb_pi = TB_k(p.index,i.index,Ktot_pi,k12_pi,s2_p,s2_i,S12_2_pi)
            tb_qi = TB_k(q.index,i.index,Ktot_qi,k12_qi,s2_q,s2_i,S12_2_qi)

            v_pq_t = self.Minn_me.cal_V_neu_as(tb_pi,tb_qi)

            v_pq += v_pq_t

        val = e_pq + v_pq
        # if(abs(val)>0.1):
        #
        #     print("val : ",val,"\t e_pq : ",e_pq,"\t v_pq : ",v_pq)

        return val

    def build_f_pq(self):
        print("self.A_size : ",self.sp_A_size,"self.nu_size : ",self.sp_nu_size)
        print("self.sps_k.fermi_index : ",self.sps_k.fermi_index)

        # nn
        #print(self.sp_A_list)
        for i in self.sp_A_list:
            for j in self.sp_A_list:
                #print("ii , jj ",ii,jj)
                val=self.cal_f_pq(i,j)
                #print(i,j,val)
                self.f_me_hh[i][j]=val
        # pp
        #print(self.sp_nu_list)
        for i in self.sp_nu_list:
            for j in self.sp_nu_list:
                ii = i + self.sps_k.fermi_index+1
                jj = j + self.sps_k.fermi_index+1
                val=self.cal_f_pq(ii,jj)
                #if(abs(val) > 0):
                    #print(ii,jj,val)
                self.f_me_pp[i][j]=val

        # print(" f_me_hh : ")
        # print(self.f_me_hh.nonzero())
        # print(" f_me_pp : ")
        # print(self.f_me_pp.nonzero())

    def build_V(self):
        # V_me_pphh V_me_pppp V_me_hhhh
        # [channel][tb_i][tb_j]

        print("self.tb_k.state_hh : ",len(self.tb_k.state_hh))
        #print("ch_list : ",ch_list)
        # ch = 0
        # for tb_t in self.tb_k.state_hh:
        #     if(len(tb_t)>0):
        #         print("ch : ",ch,"\t state_hh[ch].size : ",len(tb_t))
        #     ch+=1
        V_hhhh_list = []
        V_pphh_list = []
        V_pppp_list = []


        for ch in self.ch_list:
            # hhhh
            #print("ch : ",ch)
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])
            # print("size_hh : ", size_hh)
            # print("size_pp : ", size_pp)

            #self.V_hhhh[ch]=self.V_hhhh[ch].reshape(size_hh)

            v_hhhh_t = np.zeros((size_hh,size_hh))
            nonzero_cont = 0
            for tb_hh_f in prange(self.tb_k.state_hh[ch]):#[0:4]:
                index_f = tb_hh_f.index
                for tb_hh_i in self.tb_k.state_hh[ch]:#[0:4]:
                    index_i = tb_hh_i.index
                    val = self.Minn_me.cal_V_neu_as(tb_hh_f,tb_hh_i)
                    # if val != 0.0 and ch == 0:
                    #     print("ch : ",ch,"\t index_f : ",index_f,"\t index_i : ",index_i,"\t val : ",val)
                    #     a_i = self.tb_k.state_hh[ch][index_f].index_1
                    #     b_i = self.tb_k.state_hh[ch][index_f].index_2
                    #     i_i = self.tb_k.state_hh[ch][index_i].index_1
                    #     j_i = self.tb_k.state_hh[ch][index_i].index_2
                    #     a_xyz = self.sps_k.state[a_i].k_xyz
                    #     b_xyz = self.sps_k.state[b_i].k_xyz
                    #     i_xyz = self.sps_k.state[i_i].k_xyz
                    #     j_xyz = self.sps_k.state[j_i].k_xyz
                    #     a_xyz[np.argwhere(a_xyz<0)]=-1
                    #     a_xyz[np.argwhere(a_xyz>0)]=+1
                    #     b_xyz[np.argwhere(b_xyz<0)]=-1
                    #     b_xyz[np.argwhere(b_xyz>0)]=+1
                    #     i_xyz[np.argwhere(i_xyz<0)]=-1
                    #     i_xyz[np.argwhere(i_xyz>0)]=+1
                    #     j_xyz[np.argwhere(j_xyz<0)]=-1
                    #     j_xyz[np.argwhere(j_xyz>0)]=+1
                    #     nonzero_cont +=1
                    #
                    #     print(a_xyz,b_xyz,i_xyz,j_xyz)
                    #print(index_f,index_i)
                    v_hhhh_t[index_f][index_i] = val
            V_hhhh_list.append(v_hhhh_t)
            #print("++++++ nonzero : ",nonzero_cont)

            #sys.exit()
            # pphh

            v_pphh_t = np.zeros((size_pp,size_hh))
            for tb_pp_f in prange(self.tb_k.state_pp[ch]):
                index_f = tb_pp_f.index
                for tb_hh_i in self.tb_k.state_hh[ch]:
                    index_i = tb_hh_i.index
                    val = self.Minn_me.cal_V_neu_as(tb_pp_f,tb_hh_i)
                    v_pphh_t[index_f][index_i] = val
            V_pphh_list.append(v_pphh_t)

            # pppp
            #size_ch_pp = len(self.tb_k.state_pp[ch])
            v_pppp_t = np.zeros((size_pp,size_pp))
            for tb_pp_f in prange(self.tb_k.state_pp[ch]):
                index_f = tb_pp_f.index
                for tb_pp_i in self.tb_k.state_pp[ch]:
                    index_i = tb_pp_i.index
                    val = self.Minn_me.cal_V_neu_as(tb_pp_f,tb_pp_i)
                    v_pppp_t[index_f][index_i] = val
            V_pppp_list.append(v_pppp_t)
            # =========++++++++========== #

        self.V_hhhh = np.array(V_hhhh_list)
        self.V_pphh = np.array(V_pphh_list)
        self.V_pppp = np.array(V_pppp_list)
        # print("V_hhhh len : ", len(self.V_hhhh))
        # for ch in self.ch_list:
        #     if( self.V_pphh[ch].shape[0]>0):
        #         print("ch : ",ch,"V_pphh shape : ", self.V_pphh[ch].shape)
        # v_t = self.V_pphh[0]
        # print(v_t[v_t.nonzero()])
        # quit()

        # for ch in ch_list:
        #     if( self.V_pphh[ch].shape[0]>0):
        #         print("ch : ",ch,"V_pphh shape : ", self.V_pphh[ch].shape)
        # for ch in ch_list:
        #     if( self.V_pppp[ch].shape[0]>0):
        #         print("ch : ",ch,"V_pppp shape : ", self.V_pppp[ch].shape)

    def build_V_2(self):
        # V_me_pphh V_me_pppp V_me_hhhh
        # [channel][tb_i][tb_j]

        print("self.tb_k.state_hh : ",len(self.tb_k.state_hh))

        V_hhhh_list = []
        V_pphh_list = []
        V_pppp_list = []


        for ch in self.ch_list:
            # hhhh
            #print("ch : ",ch)
            size_hh = len(self.tb_k.state_hh[ch])
            size_hh_2 = size_hh/2 +1
            size_pp = len(self.tb_k.state_pp[ch])
            size_pp_2 = size_pp/2 +1
            # print("size_hh : ", size_hh)
            # print("size_pp : ", size_pp)

            #self.V_hhhh[ch]=self.V_hhhh[ch].reshape(size_hh)

            v_hhhh_t = np.zeros((size_hh,size_hh))
            nonzero_cont = 0
            cal_flag = []
            for tb_hh_f in self.tb_k.state_hh[ch]:#[0:4]:
                index_f = tb_hh_f.index
                for tb_hh_i in self.tb_k.state_hh[ch]:#[0:4]:
                    index_i = tb_hh_i.index
                    flag_t1 = (index_f,index_i)
                    flag_t2 = (index_i,index_f)
                    if flag_t1 in cal_flag:
                        continue
                    cal_flag.append(flag_t1)
                    cal_flag.append(flag_t2)

                    val = self.Minn_me.cal_V_neu_as(tb_hh_f,tb_hh_i)

                    v_hhhh_t[index_f][index_i] = val
                    v_hhhh_t[index_i][index_f] = val

            V_hhhh_list.append(v_hhhh_t)
            #print("++++++ nonzero : ",nonzero_cont)

            #sys.exit()
            # pphh

            v_pphh_t = np.zeros((size_pp,size_hh))
            for tb_pp_f in self.tb_k.state_pp[ch]:
                index_f = tb_pp_f.index
                for tb_hh_i in self.tb_k.state_hh[ch]:
                    index_i = tb_hh_i.index
                    val = self.Minn_me.cal_V_neu_as(tb_pp_f,tb_hh_i)
                    v_pphh_t[index_f][index_i] = val
            V_pphh_list.append(v_pphh_t)

            # pppp
            #size_ch_pp = len(self.tb_k.state_pp[ch])
            cal_flag = []
            v_pppp_t = np.zeros((size_pp,size_pp))
            for tb_pp_f in self.tb_k.state_pp[ch]:
                index_f = tb_pp_f.index
                for tb_pp_i in self.tb_k.state_pp[ch]:
                    index_i = tb_pp_i.index
                    # flag_t1 = (index_f,index_i)
                    # flag_t2 = (index_i,index_f)
                    # if flag_t1 in cal_flag:
                    #     continue
                    # cal_flag.append(flag_t1)
                    # cal_flag.append(flag_t2)
                    val = self.Minn_me.cal_V_neu_as(tb_pp_f,tb_pp_i)
                    v_pppp_t[index_f][index_i] = val
                    v_pppp_t[index_i][index_f] = val
            V_pppp_list.append(v_pppp_t)
            # =========++++++++========== #

        self.V_hhhh = np.array(V_hhhh_list)
        self.V_pphh = np.array(V_pphh_list)
        self.V_pppp = np.array(V_pppp_list)

    def build_T0(self):

        #print(f_dom_v)

        f_dom_l = []
        for ch in self.ch_list:
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])

            f_dom_ch = np.zeros((size_pp,size_hh))
            for tb_pp in self.tb_k.state_pp[ch]:
                index_pp = tb_pp.index
                a = tb_pp.index_1
                b = tb_pp.index_2
                a_i = a - self.sps_k.fermi_index - 1
                b_i = b - self.sps_k.fermi_index - 1

                e_a = self.f_me_pp[a_i][a_i]
                e_b = self.f_me_pp[b_i][b_i]
                # e_a = self.sps_k.state[a].E
                # e_b = self.sps_k.state[b].E

                e_ab = e_a + e_b
                for tb_hh in self.tb_k.state_hh[ch]:
                    index_hh = tb_hh.index
                    i = tb_hh.index_1
                    j = tb_hh.index_2
                    e_i = self.f_me_hh[i][i]
                    e_j = self.f_me_hh[j][j]
                    # e_i = self.sps_k.state[i].E
                    # e_j = self.sps_k.state[j].E
                    e_ij = e_i + e_j
                    val =  e_ij - e_ab
                    f_dom_ch[index_pp][index_hh] = val
                    if(val == 0.0):
                        print("Wrong happened at build_T0 f_dom_ch = 0")
                        print(a,b,i,j)
                        sys.exit()
            f_dom_l.append(f_dom_ch)
        f_dom_v = np.array(f_dom_l)
        # for ch in self.ch_list:
        #     if( f_dom_v[ch].shape[0]>0):
        #         print("ch : ",ch,"f_dom_v shape : ", f_dom_v[ch].shape )
        #         f_t = f_dom_v[ch]
        #         if ch == 0:
        #             print(f_t[f_t.nonzero()])

        #f_dom_flag = f_dom_v.nonzero()
        #f_dom_num = f_dom_v.shape[0]*f_dom_v.shape[1]*f_dom_v.shape[2]*f_dom_v.shape[3]

        self.f_dom_abij=f_dom_v

        T_pphh_ch = []
        for ch in self.ch_list:
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])

            T_pphh_t = np.zeros((size_pp,size_hh))
            for tb_pp in self.tb_k.state_pp[ch]:
                index_pp = tb_pp.index
                a = tb_pp.index_1 - self.sps_k.fermi_index-1
                b = tb_pp.index_2 - self.sps_k.fermi_index-1

                for tb_hh in self.tb_k.state_hh[ch]:
                    index_hh = tb_hh.index
                    i = tb_hh.index_1
                    j = tb_hh.index_2
                    v_t = self.V_pphh[ch][index_pp][index_hh]
                    f_dom = self.f_dom_abij[ch][index_pp][index_hh]
                    if(f_dom == 0):
                        print("Wrong happened at build_T0 ? f_dom = 0")
                        sys.exit()
                    T_pphh_t[index_pp][index_hh] = v_t / (f_dom*1.0)
            T_pphh_ch.append(T_pphh_t)
        # --------------- for ch ----------------- #
        self.T_pphh = np.array(T_pphh_ch)
        # for ch in ch_list:
        #     if( self.T_pphh[ch].shape[0]>0):
        #         print("ch : ",ch,"T_pppp shape : ", self.T_pphh[ch].shape )


    def build_Hbar(self):

        #h_1 = self.V_pphh

        h_2_l = []
        h_3_l = []
        h_4_l = []

        for ch in self.ch_list:
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])

            h_2_ch = np.zeros((size_pp,size_hh))
            h_3_ch = np.zeros((size_pp,size_hh))
            h_4_ch = np.zeros((size_pp,size_hh))


            if(size_hh == 0 or size_pp == 0):
                h_2_l.append(h_2_ch)
                h_3_l.append(h_3_ch)
                h_4_l.append(h_4_ch)
                continue
            else:
                h_2_ch = -1*self.f_dom_abij[ch] * self.T_pphh[ch]
                h_3_ch = 0.5*np.dot(self.V_pppp[ch],self.T_pphh[ch])
                h_4_ch = 0.5*np.dot(self.T_pphh[ch],self.V_hhhh[ch])
                h_2_l.append(h_2_ch)
                h_3_l.append(h_3_ch)
                h_4_l.append(h_4_ch)
        #-----
        h_2_vec=np.array(h_2_l)
        h_3_vec=np.array(h_3_l)
        h_4_vec=np.array(h_4_l)
        H_pphh_l = []
        for ch in self.ch_list:
            #H_pphh_ch = np.zeros((size_pp,size_hh))
            H_pphh_ch = self.V_pphh[ch]+h_2_vec[ch]+h_3_vec[ch]+h_4_vec[ch]
            H_pphh_l.append(H_pphh_ch)
        H_pphh_vec = np.array(H_pphh_l)

        return H_pphh_vec

    def cal_New_T2_pphh(self,H_bar):
        # t_abij_new = t_abij + (H_bar[ab][ij])/(f+f-f-f)
        #T_pphh_New = np.zeros(self.T_pphh.shape)
        T_l = []
        for ch in self.ch_list:
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])
            T_ch = np.zeros((size_pp,size_hh))
            T_ch_new= self.T_pphh[ch] + np.divide(H_bar[ch],self.f_dom_abij[ch])
            T_ch = 0.2*T_ch_new + 0.8*self.T_pphh[ch]
            T_l.append(T_ch)
        T_vec = np.array(T_l)
        self.T_pphh = T_vec
        #self.T_pphh = (0.2*T_pphh_New + 0.8*self.T_pphh)
        #print("  ====== T_pphh : ======= ")
        #print(self.T_pphh[self.T_pphh.nonzero()])

    def iter(self):
        self.iter_flag = 1
        H_bar = self.build_Hbar()
        #A=np.reshape(H_bar,(H_bar.size,))
        H_sum = 0.0
        for a in H_bar:
            for b in a:
                for i in b:
                    H_sum += i*i
        print(H_sum)
        eps = 0.00000001
        time_i = 0
        #E_tot = self.E_tot_cal()

        while (time_i < 1000):
            self.cal_New_T2_pphh(H_bar)
            H_bar = self.build_Hbar()
            H_sum = 0.0
            for a in H_bar:
                for b in a:
                    for i in b:
                        H_sum += i*i

            print("H_sum : ",H_sum, "Ec_cal : ", self.Ec_cal())
            if(H_sum<eps or H_sum == np.nan):
                 break
            # E_tot_new = self.E_tot_cal()
            # E_tot_diff = abs(E_tot - E_tot_new)
            # E_tot = E_tot_new
            # print("E_tot_diff",E_tot_diff)
            # if(E_tot_diff<eps):
            #     break
            time_i+=1

            if(H_sum > 10000):
                self.iter_flag = -1
                break

    def build_Hbar_ch(self,T_ch,ch):

        size_hh = len(self.tb_k.state_hh[ch])
        size_pp = len(self.tb_k.state_pp[ch])

        if(size_hh == 0 or size_pp == 0):
            H_ch = np.zeros((size_pp,size_hh))
            return H_ch

        h_2_ch = -1*self.f_dom_abij[ch] * T_ch
        h_3_ch = 0.5*np.dot(self.V_pppp[ch],T_ch)
        h_4_ch = 0.5*np.dot(T_ch,self.V_hhhh[ch])

        H_pphh_ch = self.V_pphh[ch]+h_2_ch+h_3_ch+h_4_ch

        return H_pphh_ch

    def cal_New_T2_pphh_ch(self,H_ch,T_ch,ch):
        size_hh = len(self.tb_k.state_hh[ch])
        size_pp = len(self.tb_k.state_pp[ch])
        #T_ch = np.zeros((size_pp,size_hh))
        T_ch_new= T_ch + np.divide(H_ch,self.f_dom_abij[ch])
        T_ch_new_mix = 0.2*T_ch_new + 0.8*T_ch
        return T_ch_new_mix

    def iter_ch(self):
        self.iter_flag = 1

        eps = 0.00001
        T_pphh_l = []
        for ch in self.ch_list:
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])

            if(size_hh == 0 or size_pp == 0):
                T_ch = np.zeros((size_hh,size_pp))
                T_pphh_l.append(T_ch)
                continue
            time_i = -1

            while (time_i < 1000):
                time_i += 1
                if time_i == 0:
                    H_bar_ch = self.build_Hbar_ch(self.T_pphh[ch],ch)
                    T_ch = self.cal_New_T2_pphh_ch(H_bar_ch,self.T_pphh[ch],ch)
                else:
                    H_bar_ch = self.build_Hbar_ch(T_ch,ch)
                    T_ch_new = self.cal_New_T2_pphh_ch(H_bar_ch,T_ch,ch)
                    T_ch = T_ch_new
                H_sum = np.einsum("ij,ij->",H_bar_ch,H_bar_ch)
                #print("ch : ",ch," H_sum : ",H_sum)
                if(H_sum<eps or H_sum == np.nan):
                    T_pphh_l.append(T_ch)
                    break
            #quit()
        self.T_pphh = np.array(T_pphh_l)

        print("H_sum : ",H_sum, "Ec_cal : ", self.Ec_cal())



        #A=np.reshape(H_bar,(H_bar.size,))

    def Ec_cal(self):
        if(self.iter_flag < 0):
            return np.nan
        Ec_1 = 0
        Ec_2 = 0.0
        Ec_3 = 0.0
        for ch in self.ch_list:
            size_hh = len(self.tb_k.state_hh[ch])
            size_pp = len(self.tb_k.state_pp[ch])
            for ab in np.arange(0,size_pp):
                for ij in np.arange(0,size_hh):
                    sum_ch = 1.0/4*self.T_pphh[ch][ab][ij]*self.V_pphh[ch][ab][ij]
                    Ec_2 += sum_ch
            #sum_ch = 1.0/4*np.einsum('ab,ab->',self.T_pphh[ch],self.V_pphh[ch])

        Ec = Ec_1 + Ec_2 + Ec_3
        return Ec

    def Etot_cal(self):
        e_ii = 0
        for i in self.tb_k.A_list:
            e_ii += self.sps_k.state[i].E
        e_v = 0

        #e_v =
        for ch in self.ch_list:
            val = np.trace(self.V_hhhh[ch])
            #print("### \t --- ch : ",ch,"t\t +e_v = ",val)
            e_v += val

        E_ref = e_ii + 0.5*e_v
        print("------ e_ii/14 = ",e_ii/14," ---- 0.5*e_v = ",0.5*e_v)
        # E_ref = 0# = self.f_me_hh.trace()
        # for i in self.tb_k.A_list:
        #     E_ref += self.f_me_hh[i][i]

        E_corr = self.Ec_cal()

        E_tot = E_ref + E_corr
        self.E_ref = E_ref
        self.E_corr = E_corr
        self.E_tot = E_tot
        return E_tot

    def E_HF_test(self):
        e_ii = 0
        for i in self.tb_k.A_list:
            e_ii += self.sps_k.state[i].E
        e_v = 0
        for i in self.sp_A_list:
            a= self.sps_k.state[i]
            for j in self.sp_A_list:
                b= self.sps_k.state[j]
                if(i == j):
                    continue
                a_xyz = a.k_xyz
                b_xyz = b.k_xyz

                ab_xyz = a_xyz + b_xyz

                s2_1 = a.s2_z
                s2_2 = b.s2_z
                k12 = a_xyz - b_xyz
                S12_2 = a.s2_z + b.s2_z
                tb_t = TB_k(a.index,b.index,ab_xyz,k12,s2_1,s2_2,S12_2)
                v_t = self.Minn_me.cal_V_neu_as(tb_t,tb_t)
                e_v+=v_t
                #print("i : ",i,a_xyz,s2_1,b_xyz,s2_2,v_t)
                #quit()

        E_ref = e_ii + 0.5*e_v
        print("------ e_ii/14 = ",e_ii/14," ---- 0.5*e_v = ",0.5*e_v)
        print(" +++ E_ref/14 : ",E_ref/14)
















