#!/usr/bin/python
import numpy as np
import itertools
import sys
#import matplotlib.pyplot as plt
from basis import *
from slater import *
from pairing_me import *
#from hamiltonian import *

class TwoBody:
    #c
    #'Two Body state'
    def __init__(self,a,b,index):
        self.index_a = a
        self.index_b = b
        self.index = index


class state_TwoBody:
    # two body state
    def __init__(self,sp,pair_flag):
        self.sp = sp
        self.sp_size = len(sp.state)
        self.TwP_state = []
        self.pair_flag = pair_flag
    def build(self):
        index = 0
        for sp_t_a in self.sp.state:
            p_a = sp_t_a.p
            index_a = sp_t_a.index
            for sp_t_b in self.sp.state:
                p_b = sp_t_b.p
                index_b = sp_t_b.index
                if(index_a == index_b):
                    continue
                if(self.pair_flag == 1 and p_a != p_b):
                    continue
                TwP = TwoBody(index_a,index_b,index)
                self.TwP_state.append(TwP)
                index += 1
        self.fermi_twb = self.find_fermi_twb()
        print("self.fermi_twb : ",self.fermi_twb)
        print("self.sp.fermi_index : ",self.sp.fermi_index)

    def find_fermi_twb(self):
        # find pair ij <= fermi_twb and ab > fermi_twb
        index = 0
        fermi_twb = 0
        fermi_sp = self.sp.fermi_index
        for twb in self.TwP_state:
            fermi_twb = index
            index +=1
            index_a = twb.index_a
            index_b = twb.index_b
            #print(index_a,index_b)
            if(index_a> fermi_sp or index_b > fermi_sp):
                break

        return fermi_twb-1


    def print_TwoB(self):
        for twp in self.TwP_state:
            print(twp.index, "\t (",twp.index_a,", ",twp.index_b,")" )



class CCSD_ME:
    def __init__(self,sp,V_me):
        self.sp = sp
        #self.twb = twb
        self.V_me = V_me
        self.sp_size = len(sp.state)
        self.nu_size = len(sp.state)-sp.fermi_index-1
        self.A_size = sp.fermi_index+1

        #self.size_twb = len(twb.TwP_state)
        #self.size_ij = twb.fermi_twb
        #self.size_ab = len(twb.TwP_state) - self.size_ij

        self.T_pphh = np.zeros((self.nu_size,self.nu_size,self.A_size,self.A_size))
        self.V_pphh = np.zeros((self.nu_size,self.nu_size,self.A_size,self.A_size))
        self.V_hhhh = np.zeros((self.A_size,self.A_size,self.A_size,self.A_size))
        self.V_pppp = np.zeros((self.nu_size,self.nu_size,self.nu_size,self.nu_size))
        self.f_me_hh = np.zeros((self.A_size,self.A_size))
        self.f_me_pp = np.zeros((self.nu_size,self.nu_size))
        self.f_dom_abij = np.zeros((self.nu_size,self.nu_size,self.nu_size,self.nu_size))
        self.H_bar = np.zeros((self.nu_size,self.nu_size,self.A_size,self.A_size))

        self.A_list = np.arange(0,self.A_size)
        self.nu_list = np.arange(0,self.nu_size)
        print("A_list : ",self.A_list)
        print("nu_list : ",self.nu_list)


        self.build_f_pq()
        self.build_V()
        self.build_T0()


    def cal_f_pq(self,index_1,index_2):
        val = 0
        e_pq = 0.0
        if(index_1 == index_2):
            e_pq = self.sp.state[index_1].energy
        #i_list = np.arange(0,self.sp.fermi_index+1)
            #print("index_1",index_1,"E : ",self.sp.state[index_1].energy)
        v_pq = 0.0
        for i in self.A_list:
            conf = str(index_1) +"-"+ str(i) +"-"+str(index_2)+"-" +str(i)
            v_pq_t = 0
            if conf in self.V_me.me_dic:
                v_pq_t = self.V_me.me_dic[conf]
            #print("index_1: ",index_1,"\t i ",i ,"\t index_2: ",index_2,"\t val",v_pq_t)
            v_pq += v_pq_t

        val = e_pq + v_pq
        # if(abs(val)>0.1):
        #
        #     print("val : ",val,"\t e_pq : ",e_pq,"\t v_pq : ",v_pq)

        return val

    def build_f_pq(self):
        print("self.A_size : ",self.A_size,"self.nu_size : ",self.nu_size)
        print("self.sp.fermi_index : ",self.sp.fermi_index)

        print(self.A_list)
        for i in self.nu_list:
            for j in self.nu_list:
                ii = i + self.sp.fermi_index+1
                jj = j + self.sp.fermi_index+1
                val=self.cal_f_pq(ii,jj)
                if(abs(val) > 0):
                    print(ii,jj,val)
                self.f_me_pp[i][j]=val

        print(self.nu_list)
        for i in self.A_list:
            for j in self.A_list:
                #print("ii , jj ",ii,jj)
                val=self.cal_f_pq(i,j)
                #print(i,j,val)
                self.f_me_hh[i][j]=val
        #print(" f_me_hh : ")
        #print(self.f_me_hh)
        #print(" f_me_pp : ")
        #print(self.f_me_pp)



    def build_V(self):
        # V_me_pphh V_me_pppp V_me_hhhh
        for a in self.nu_list:
            aa = a + self.sp.fermi_index + 1
            for b in self.nu_list:
                bb = b + self.sp.fermi_index + 1
                for i in self.A_list:
                    for j in self.A_list:
                        conf = str(aa)+"-" + str(bb)+"-" +str(i)+"-" +str(j)
                        v_t = 0.0
                        if conf in self.V_me.me_dic:
                            v_t = self.V_me.me_dic[conf]
                        self.V_pphh[a][b][i][j]=v_t

        for a in self.nu_list:
            aa = a + self.sp.fermi_index + 1
            for b in self.nu_list:
                bb = b + self.sp.fermi_index + 1
                for c in self.nu_list:
                    cc = c + self.sp.fermi_index + 1
                    for d in self.nu_list:
                        dd = d + self.sp.fermi_index + 1
                        conf = str(aa)+"-" + str(bb)+"-" +str(cc)+"-" +str(dd)
                        v_t = 0.0
                        if conf in self.V_me.me_dic:
                            v_t = self.V_me.me_dic[conf]
                        self.V_pppp[a][b][c][d]=v_t

        for i in self.A_list:
            for j in self.A_list:
                for k in self.A_list:
                    for l in self.A_list:
                        conf = str(i)+"-" + str(j)+"-" +str(k)+"-" +str(l)
                        v_t = 0.0
                        if conf in self.V_me.me_dic:
                            v_t = self.V_me.me_dic[conf]
                        self.V_hhhh[i][j][k][l]=v_t

        #print("self.V_pphh")
        #print(self.V_pphh)
        # print("self.V_pppp")
        # print(self.V_pppp)
        # print("self.V_hhhh")
        # print(self.V_hhhh)

    def build_T0(self):

        f_dom_v = np.zeros(self.T_pphh.shape)
        for a in self.nu_list:
            for b in self.nu_list:
                for i in self.A_list:
                    for j in self.A_list:
                        val= self.f_me_hh[i][i]+self.f_me_hh[j][j]-self.f_me_pp[a][a]-self.f_me_pp[b][b]
                        f_dom_v[a][b][i][j] =val
                        if(abs(val) <= 0.0 ):
                            print("Wrong happened at f_dom_num val : ",val,'\t a:',a,'\t b:',b,'\t i:',i,"\t j:",j)
                            sys.exit()
        f_dom_flag = f_dom_v.nonzero()
        f_dom_num = f_dom_v.shape[0]*f_dom_v.shape[1]*f_dom_v.shape[2]*f_dom_v.shape[3]

        self.f_dom_abij=f_dom_v
        #print(f_dom_v)

        for a in self.nu_list:
            aa = a + self.sp.fermi_index + 1
            for b in self.nu_list:
                bb = b + self.sp.fermi_index + 1
                for i in self.A_list:
                    for j in self.A_list:
                        conf = str(aa)+"-" + str(bb)+"-" +str(i)+"-" +str(j)
                        v_t = 0
                        if conf in self.V_me.me_dic:
                            v_t = self.V_me.me_dic[conf]
                        #f_dom = self.f_me_pp[a][a]+self.f_me_pp[b][b]-self.f_me_hh[i][i]-self.f_me_hh[j][j]
                        f_dom = self.f_dom_abij[a][b][i][j]
                        if(f_dom == 0):
                            print("Wrong happened at build_T0 f_dom = 0")
                            sys.exit()
                        self.T_pphh[a][b][i][j] = v_t / (f_dom*1.0)
        #print(self.T_pphh)
        #print("  ====== T_pphh : ======= ")
        #print(self.T_pphh[self.T_pphh.nonzero()])


    def build_Hbar(self):

        h_1 = self.V_pphh

        h_2_1 = np.einsum('acij,bc->abij',self.T_pphh,self.f_me_pp)
        h_2_2 = -1*h_2_1.transpose((1,0,2,3))
        h_2 = (h_2_1 + h_2_2)

        h_3_1 = np.einsum('abik,kj->abij',self.T_pphh,self.f_me_hh)
        h_3_2 = -1*h_3_1.transpose((0,1,3,2))
        h_3 = -1*(h_3_1 + h_3_2)

        # h_2 = -1*self.T_pphh*self.f_dom_abij
        # h_3 = np.zeros(h_1.shape)

        h_4 = 0.5*np.einsum('abcd,cdij->abij',self.V_pppp,self.T_pphh)

        h_5 = 0.5*np.einsum('abkl,klij->abij',self.T_pphh,self.V_hhhh)

        h_6 = np.zeros(h_1.shape)

        # h_6_1 = np.einsum('acik,kbcj',self.T_pphh,self.V_hhhh)
        # h_6_2 = -1*h_6_1.transpose((1,0,2,3))
        # h_6_3 = -1*h_6_1.transpose((0,1,3,2))
        # h_6_4 = h_6_1.transpose((0,1,3,2))
        # h_6 = h_6_1 + h_6_2 + h_6_3 + h_6_4

        h_7_A = np.einsum('acik,cdkl->adil',self.T_pphh,self.V_pphh)
        h_7_1 = np.einsum('adil,dblj->abij',h_7_A,self.T_pphh)
        h_7_2 = -1*h_7_1.transpose((1,0,2,3))
        h_7_3 = -1*h_7_1.transpose((0,1,3,2))
        h_7_4 = h_7_1.transpose((1,0,3,2))
        h_7 = 0.5*(h_7_1 +h_7_2 +h_7_3+h_7_4)

        # h_7_A = np.einsum('acik,cdkl->adil',self.T_pphh,self.V_pphh)
        # h_7_1 = np.einsum('adil,bdjl->abij',h_7_A,self.T_pphh)
        # h_7_2 = -1*h_7_1.transpose((0,1,3,2))
        # h_7 = (h_7_1 + h_7_2)

        h_8_A = np.einsum('cdik,cdkl->il',self.T_pphh,self.V_pphh)
        h_8_1 = np.einsum('il,ablj->abij',h_8_A,self.T_pphh )
        h_8_2 = -1*h_8_1.transpose((0,1,3,2))
        h_8 = 0.5*(h_8_1+h_8_2)

        h_9_A = np.einsum('ackl,cdkl->ad',self.T_pphh,self.V_pphh)
        h_9_1 = np.einsum('ad,dbij->abij',h_9_A,self.T_pphh)
        h_9_2 = -1*h_9_1.transpose((1,0,2,3))
        h_9 = 0.5*(h_9_1 + h_9_2)

        h_10_A = np.einsum('abkl,cdkl->abcd',self.T_pphh,self.V_pphh)
        h_10 = 1.0/4 * np.einsum('abcd,cdij->abij',h_10_A,self.T_pphh)

        H_bar = (h_1+h_2+h_3+h_4+h_5+h_6+h_7+h_8+h_9+h_10)
        #self.H_bar = H_bar

        # print("  ====== self.self.T_pphh : ======= ")
        # print(self.T_pphh[self.T_pphh.nonzero()])
        # print("  ====== self.f_me_pp : ======= ")
        # print(self.f_me_pp[self.f_me_pp.nonzero()])
        #
        # print("  ====== h_1 : ======= ")
        # print(h_1[h_1.nonzero()])
        # print("  ====== h_2_1 : ======= ")
        # print(h_2_1[h_2_1.nonzero()])
        # print("  ====== h_2 : ======= ")
        # print(h_2[h_2.nonzero()])
        # print("  ====== h_3_1 : ======= ")
        # print(h_3_1[h_3_1.nonzero()])
        # print("  ====== h_3_2 : ======= ")
        # print(h_3_2[h_3_2.nonzero()])
        # print("  ====== h_3 : ======= ")
        # print(h_3[h_3.nonzero()])
        # print("  ====== h_4 : ======= ")
        # print(h_4[h_4.nonzero()])
        # print("  ====== h_5 : ======= ")
        # print(h_5[h_5.nonzero()])
        # print("  ====== h_6 : ======= ")
        # print(h_6[h_6.nonzero()])
        # print("  ====== h_7 : ======= ")
        # print(h_7[h_7.nonzero()])
        # print("  ====== h_8 : ======= ")
        # print(h_8[h_8.nonzero()])
        # print("  ====== h_8_1 : ======= ")
        # print(h_8_1[h_8_1.nonzero()])
        # print("  ====== h_8_2 : ======= ")
        # print(h_8_2[h_8_2.nonzero()])
        # print("  ====== h_9 : ======= ")
        # print(h_9[h_9.nonzero()])
        # print("  ====== h_9_1 : ======= ")
        # print(0.5*h_9_1[h_9_1.nonzero()])
        # print("  ====== h_10 : ======= ")
        # print(h_10[h_10.nonzero()])
        # print("  ====== H_bar : ======= ")
        # print(H_bar[H_bar.nonzero()])


        #
        # h_2_1_check = self.ein_sum_A2_B4_23(self.f_me_pp,self.T_pphh)
        # h_2_2_check = -1*self.ein_sum_trans_12(h_2_1)
        #
        # h_3_1_check = self.ein_sum_A2_B4_14(self.f_me_pp,self.T_pphh)
        # h_3_2_check = -1*self.ein_sum_trans_34(h_3_1)
        #
        # print((h_2_1 == h_2_1_check).all())
        # print((h_2_2 == h_2_2_check).all())
        # print((h_3_1 == h_3_1_check).all())
        # print((h_3_2 == h_3_2_check).all())
        # print("h_1")
        # print(h_1)
        # print("h_2_1")
        # print(h_2_1)
        # print("h_2_2")
        # print(h_2_2)
        return H_bar

    def build_Hbar2(self):

        h_1 = self.V_pphh

        h_2 = -1*self.T_pphh*self.f_dom_abij

        h_3 = 0.5*np.einsum('abcd,cdij->abij',self.V_pppp,self.T_pphh)

        h_4 = 0.5*np.einsum('abkl,klij->abij',self.T_pphh,self.V_hhhh)

        h_5 = np.zeros(h_1.shape)

        h_6_A = np.einsum('abkl,cdkl->abcd',self.T_pphh,self.V_pphh)
        h_6 = 1.0/4 * np.einsum('abcd,cdij->abij',h_6_A,self.T_pphh)

        h_7_A = np.einsum('acik,cdkl->adil',self.T_pphh,self.V_pphh)
        h_7_1 = np.einsum('adil,bdjl->abij',h_7_A,self.T_pphh)
        h_7_2 = -1*h_7_1.transpose((0,1,3,2))
        h_7 = (h_7_1 + h_7_2)

        h_8_A = np.einsum('dcik,cdkl->il',self.T_pphh,self.V_pphh)
        h_8_1 = np.einsum('il,ablj->abij',h_8_A,self.T_pphh )
        h_8_2 = -1*h_8_1.transpose((0,1,3,2))
        h_8 = -0.5*(h_8_1 + h_8_2)

        h_9_A = np.einsum('aclk,cdkl->ad',self.T_pphh,self.V_pphh)
        h_9_1 = np.einsum('ad,dbij->abij',h_9_A,self.T_pphh)
        h_9_2 = -1*h_9_1.transpose((1,0,2,3))
        h_9 = -0.5*(h_9_1 + h_9_2)

        H_bar = (h_1+h_2+h_3+h_4+h_5+h_6+h_7+h_8+h_9)

        return H_bar

    def cal_New_T2_pphh(self,H_bar):
        # t_abij_new = t_abij + (H_bar[ab][ij])/(f+f-f-f)
        #T_pphh_New = np.zeros(self.T_pphh.shape)
        T_pphh_New = self.T_pphh + np.divide(H_bar,self.f_dom_abij)
        self.T_pphh = (0.2*T_pphh_New + 0.8*self.T_pphh)
        #print("  ====== T_pphh : ======= ")
        #print(self.T_pphh[self.T_pphh.nonzero()])

    def iter(self):
        self.iter_flag = 1
        H_bar = self.build_Hbar()
        #A=np.reshape(H_bar,(H_bar.size,))
        H_sum = 0.0
        for a in H_bar/self.f_dom_abij:
            for b in a:
                for i in b:
                    for j in i:
                        H_sum += j*j
        #print(H_sum)
        eps = 0.00000001
        time_i = 0
        E_tot = self.E_tot_cal()

        while (time_i < 10000):
            self.cal_New_T2_pphh(H_bar)
            H_bar = self.build_Hbar()
            H_sum = 0.0
            for a in H_bar:
                for b in a:
                    for i in b:
                        for j in i:
                            H_sum += j*j

            print("H_sum : ",H_sum, "Ec_cal : ", self.Ec_cal())
            if(H_sum<eps):
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

    def Ec_cal(self):
        if(self.iter_flag < 0):
            return np.nan
        Ec_1 = 0
        Ec_2 = 1.0/4*np.einsum('abij,abij->',self.T_pphh,self.V_pphh)
        Ec_3 = 0.0
        Ec = Ec_1 + Ec_2 + Ec_3
        return Ec

    def E_tot_cal(self):
        E_ref = self.V_hhhh[0][0][0][0] + 2

        E_tot = E_ref + 1.0/4 * np.einsum('abij,abij',self.T_pphh,self.V_pphh)
        return E_tot

    def ein_sum_A2_B4_23(self,A,B):
        # f:bc  t:acij
        print("A dim : ",A.shape)
        print("B dim : ",B.shape)
        res = np.zeros(B.shape)
        l1_list = np.arange(0,len(res))
        l2_list = np.arange(0,len(res[0][0]))

        for a in l1_list:
            for b in l1_list:
                for i in l2_list:
                    for j in l2_list:
                        val = 0.0
                        for c in l2_list:
                            val += B[a][c][i][j]*A[b][c]
                        res[a][b][i][j] = val
        return res

    def ein_sum_A2_B4_14(self,A,B):
        # f:kj  t:abik
        print("A dim : ",A.shape)
        print("B dim : ",B.shape)
        res = np.zeros(B.shape)
        l1_list = np.arange(0,len(res))
        l2_list = np.arange(0,len(res[0][0]))

        for a in l1_list:
            for b in l1_list:
                for i in l2_list:
                    for j in l2_list:
                        val = 0.0
                        for k in l2_list:
                            val += B[a][b][i][k]*A[k][j]
                        res[a][b][i][j] = val
        return res

    def ein_sum_trans_12(self,A):
        res = np.zeros(A.shape)
        l_list = np.arange(0,len(res))
        for i in l_list:
            for j in l_list:
                res[i][j] = A[j][i]
        return res

    def ein_sum_trans_34(self,A):
        res = np.zeros(A.shape)
        l_list = np.arange(0,len(res))
        for i in l_list:
            for j in l_list:
                for k in l_list:
                    for l in l_list:
                        res[i][j][k][l] = A[i][j][l][k]
        return res











