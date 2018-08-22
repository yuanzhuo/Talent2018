#!/usr/bin/python
import numpy as np
import itertools
import sys
#import matplotlib.pyplot as plt
from basis import *
from pairing_me import *
#from hamiltonian import *


class IMSRG_ME:
    def __init__(self,sp,V_me):
        self.sp = sp
        #self.twb = twb
        self.V_pair = V_me
        self.fermi_index = sp.fermi_index
        self.sp_size = len(sp.state)
        self.nu_size = len(sp.state)-sp.fermi_index-1
        self.A_size = sp.fermi_index+1


        self.V_me = np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        self.f_me = np.zeros((self.sp_size,self.sp_size))

        self.V_pqrs_s = np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        self.f_pq_s = np.zeros((self.sp_size,self.sp_size))

        self.N_pqrs_s = np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        self.N_pq_s = np.zeros((self.sp_size,self.sp_size))

        self.Delta_pq_s = np.zeros((self.sp_size,self.sp_size))
        self.Delta_pqrs_s = np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))



        self.E_s = 0


        self.A_list = np.arange(0,self.A_size)
        self.nu_list = np.arange(0,self.nu_size)
        self.sp_list = np.arange(0,self.sp_size)
        self.build_f_pq()
        self.build_V()
        print("A_list : ",self.A_list)
        print("nu_list : ",self.nu_list)



    def cal_f_pq(self,index_1,index_2):
        val = 0
        e_pq = 0.0
        if(index_1 == index_2):
            e_pq = self.sp.state[index_1].energy
        #i_list = np.arange(0,self.sp.fermi_index+1)
            #print("index_1",index_1,"E : ",self.sp.state[index_1].energy)
        #return e_pq
        v_pq = 0.0
        for i in self.A_list:
            conf = str(index_1) +"-"+ str(i) +"-"+str(index_2)+"-" +str(i)
            v_pq_t = 0
            if conf in self.V_pair.me_dic:
                v_pq_t = self.V_pair.me_dic[conf]
            #print("index_1: ",index_1,"\t i ",i ,"\t index_2: ",index_2,"\t val",v_pq_t)
            v_pq += v_pq_t

        val = e_pq + v_pq
        # if(abs(val)>0.1):
        #
        # print("val : ",val,"\t e_pq : ",e_pq,"\t v_pq : ",v_pq)

        return val

    def build_f_pq(self):
        print("self.A_size : ",self.A_size,"self.nu_size : ",self.nu_size)
        print("self.sp.fermi_index : ",self.sp.fermi_index)

        for i in self.sp_list:
            for j in self.sp_list:
                #print("ii , jj ",ii,jj)
                val=self.cal_f_pq(i,j)
                # if i == j:
                #     print(i,j,val)
                self.f_me[i][j]=val
        # print(" f_me : ")
        # print(self.f_me)




    def build_V(self):
        # V_me_pphh V_me_pppp V_me_hhhh
        for a in self.sp_list:
            for b in self.sp_list:
                for i in self.sp_list:
                    for j in self.sp_list:
                        conf = str(a)+"-" + str(b)+"-" +str(i)+"-" +str(j)
                        v_t = 0.0
                        if conf in self.V_pair.me_dic:
                            v_t = self.V_pair.me_dic[conf]

                        #conf_2 = str(a)+"-" + str(b)+"-" +str(j)+"-" +str(i)
                        # v_t_2 = 0.0
                        # if conf in self.V_pair.me_dic:
                        #     v_t_2 = self.V_pair.me_dic[conf]
                        # if(v_t != 0.0):
                        #     print(a,b,i,j,"\t v_t = ",v_t)
                        self.V_me[a][b][i][j]= v_t

    def cal_E_s(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)

        # dE_1 = np.einsum("ab,ab->",self.N_ph_s[ferm_2:sp_size][0:ferm_1],self.f_me[ferm_2:sp_size][0:ferm_1])
        dE_1 = 0.0
        for i in self.sp_list:
            if(i>self.fermi_index):
                continue
            for a in self.sp_list:
                if(a<=self.fermi_index):
                    continue
                val = self.N_pq_s[i][a]*self.f_pq_s[a][i]
                dE_1 += val
                #print(a,i,"\t val = ",val)

        dE_2 = 0.0
        for i in self.sp_list:
            if(i>self.fermi_index):
                continue
            for j in self.sp_list:
                if(j>self.fermi_index or i==j ):
                    continue
                for a in self.sp_list:
                    if(a<=self.fermi_index):
                        continue
                    for b in self.sp_list:
                        if(b<=self.fermi_index or a==b ):
                            continue
                        val_N_1 = self.N_pqrs_s[i][j][a][b]
                        val_V_1 = self.V_me_s[a][b][i][j]
                        val_1 = val_N_1*val_V_1
                        val = val_1
                        #print(a,b,i,j,"\tval_N : ",val_N_1,"\t val_V : ",val_V_1,"\t val = ",val)
                        dE_2 += val #- val_2
        dE = dE_1 + 0.5*dE_2
        #dE =  0.5*dE_2
        print("dE_one = ",dE_1,"\t dE_two = ",dE_2)

        return dE


    def cal_f_pq_s(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)
        print(np.diag(self.f_me)[0:ferm_1])
        print(np.diag(self.f_me)[ferm_2:sp_size])


        #N_ph = self.N_ph_s[ferm_2:sp_size][0:ferm_1]
        #N_pphh = self.N_pphh_s[ferm_2:sp_size][ferm_2:sp_size][0:ferm_1][0:ferm_1]


        df_s = np.zeros((self.sp_size,self.sp_size))


        #df_t_2 = self.f_me[0:sp_size][ferm_2:sp_size]
        df_1_1 = np.einsum("ia,aj->ij",self.N_pq_s,self.f_me_s)
        df_1_2 = np.einsum("ij->ji",df_1_1)
        df_1 = df_1_1 + df_1_2

        V_biaj_ijab = np.einsum("biaj->ijab",self.V_me_s)
        N_biaj_ijab = np.einsum("biaj->ijab",self.N_pqrs_s)

        df_2_ij_A = np.zeros((self.sp_size,self.sp_size))
        df_2_ij_B = np.zeros((self.sp_size,self.sp_size))
        for i in self.sp_list:
            for j in self.sp_list:
                df_2_ij_A[i][j] = np.einsum("ab,ab->",self.N_pq_s[0:ferm_1][ferm_2:sp_size],V_biaj_ijab[i][j][0:ferm_1][ferm_2:sp_size])
                df_2_ij_B[i][j] = np.einsum("ab,ab->",self.f_me_s[0:ferm_1][ferm_2:sp_size],N_biaj_ijab[i][j][0:ferm_1][ferm_2:sp_size])

        df_2 = df_2_ij_A - df_2_ij_B
        #quit()

        V_abcj_jcab = np.einsum("abcj->jcab",self.V_me_s)
        N_ciab_icab = np.einsum("ciab->icab",self.N_pqrs_s)
        df_3_1=np.zeros((self.sp_size,self.sp_size))
        #df_3_2=np.zeros((self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                sum_1 = 0
                sum_2 = 0
                for c in self.sp_list:
                    if(c<=self.fermi_index):
                        sum_1 += np.einsum("ab,ab->",N_ciab_icab[i][c][ferm_2:][ferm_2:],V_abcj_jcab[j][c][ferm_2:][ferm_2:])
                    else:
                        sum_2 += np.einsum("ab,ab->",N_ciab_icab[i][c][0:ferm_1][0:ferm_1],V_abcj_jcab[j][c][0:ferm_1][0:ferm_1])
                # -----
                df_3_1[i][j] = sum_1+sum_2
                # =====
        df_3_2 = np.einsum("ij->ji",df_3_1)
        df_3 = 0.5*(df_3_1 + df_3_2)


        df_pq_s = df_1 + df_2 + df_3
        return df_pq_s

    def cal_G_pqrs_s(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)

        dG_1_1 = np.einsum("ia,ajkl->ijkl",self.N_pq_s,self.V_me_s)
        dG_1_2 = np.einsum("ia,ajkl->ijkl",self.f_me_s,self.N_pqrs_s)
        dG_1_A = dG_1_1 - dG_1_2
        dG_1_B = np.einsum("ijkl->jikl",dG_1_A)
        dG_1 = dG_1_A - dG_1_B

        dG_2_1 = np.einsum("ak,ijal->ijkl",self.N_pq_s,self.V_me_s)
        dG_2_2 = np.einsum("ak,ijal->ijkl",self.f_me_s,self.N_pqrs_s)
        dG_2_A = dG_2_1 - dG_2_2
        dG_2_B = np.einsum("ijkl->ijlk",dG_2_A)
        dG_2 = -1*(dG_2_A - dG_2_B)

        N_abkl_klab = np.einsum("abkl->klab",self.N_pqrs_s)
        V_abkl_klab = np.einsum("abkl->klab",self.V_me_s)
        dG_3_1=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_2=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_3=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_4=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                for k in self.sp_list:
                    for l in self.sp_list:
                        dG_3_1[i][j][k][l] = np.einsum("ab,ab->",self.N_pqrs_s[i][j][0:ferm_1][0:ferm_1],V_abkl_klab[k][l][0:ferm_1][0:ferm_1])
                        dG_3_2[i][j][k][l] = np.einsum("ab,ab->",self.V_me_s[i][j][0:ferm_1][0:ferm_1],N_abkl_klab[k][l][0:ferm_1][0:ferm_1])
                        dG_3_3[i][j][k][l] = np.einsum("ab,ab->",self.N_pqrs_s[i][j][ferm_2:][ferm_2:],V_abkl_klab[k][l][ferm_2:][ferm_2:])
                        dG_3_4[i][j][k][l] = np.einsum("ab,ab->",self.V_me_s[i][j][ferm_2:][ferm_2:],N_abkl_klab[k][l][ferm_2:][ferm_2:])

        dG_3_A = dG_3_1 - dG_3_2
        dG_3_B = dG_3_3 - dG_3_4

        dG_3 = 0.5*(dG_3_B - dG_3_A)

        N_aibk_ikab = np.einsum("aibk->ikab",self.N_pqrs_s)
        V_bjal_jlab = np.einsum("bjal->jlab",self.V_me_s)

        dG_4_1=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                for k in self.sp_list:
                    for l in self.sp_list:
                        dG_4_1[i][j][k][l] = np.einsum("ab,ab->",N_aibk_ikab[i][k][0:ferm_1][ferm_2:sp_size],V_bjal_jlab[j][l][0:ferm_1][ferm_2:sp_size])

        dG_4_2 = np.einsum("ijkl->ijlk",dG_4_1)

        dG_4 = 2*(dG_4_1 - dG_4_2)

        dG_pqrs_s = dG_1 + dG_2 + dG_3 +dG_4

        return dG_pqrs_s

    def cal_E_s2(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)

        # dE_1 = np.einsum("ab,ab->",self.N_ph_s[ferm_2:sp_size][0:ferm_1],self.f_me[ferm_2:sp_size][0:ferm_1])
        dE_1 = 0.0
        for i in self.sp_list:
            if(i>self.fermi_index):
                continue
            for a in self.sp_list:
                if(a<=self.fermi_index):
                    continue
                val_1 = self.N_pq_s[i][a]*self.f_pq_s[a][i]
                val_2 = self.f_pq_s[i][a]*self.N_pq_s[a][i]

                dE_1 += val_1 - val_2
                #print(a,i,"\t val = ",val)

        dE_2 = 0.0
        for i in self.sp_list:
            if(i>self.fermi_index):
                continue
            for j in self.sp_list:
                if(j>self.fermi_index or i==j ):
                    continue
                for a in self.sp_list:
                    if(a<=self.fermi_index):
                        continue
                    for b in self.sp_list:
                        if(b<=self.fermi_index or a==b ):
                            continue
                        val_N_1 = self.N_pqrs_s[i][j][a][b]
                        val_V_1 = self.V_me_s[a][b][i][j]
                        val_1 = val_N_1*val_V_1
                        val_N_2 = self.N_pqrs_s[a][b][i][j]
                        val_V_2 = self.V_me_s[i][j][a][b]
                        val_2 = val_N_2*val_V_2
                        val = val_1 - val_2
                        #print(a,b,i,j,"\tval_N : ",val_N_1,"\t val_V : ",val_V_1,"\t val = ",val)
                        dE_2 += val #- val_2
        dE = dE_1 + 0.25*dE_2
        #dE =  0.5*dE_2
        print("dE_one = ",dE_1,"\t dE_two = ",dE_2)

        return dE
    def cal_f_pq_s2(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        #print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)
        #print(np.diag(self.f_me)[0:ferm_1])
        #print(np.diag(self.f_me)[ferm_2:sp_size])


        #N_ph = self.N_ph_s[ferm_2:sp_size][0:ferm_1]
        #N_pphh = self.N_pphh_s[ferm_2:sp_size][ferm_2:sp_size][0:ferm_1][0:ferm_1]


        df_s = np.zeros((self.sp_size,self.sp_size))


        #df_t_2 = self.f_me[0:sp_size][ferm_2:sp_size]
        df_1_1 = np.einsum("ia,aj->ij",self.N_pq_s,self.f_me_s)
        df_1_2 = np.einsum("ia,aj->ij",self.f_me_s,self.N_pq_s)
        df_1 = df_1_1 - df_1_2

        V_biaj_ijab = np.einsum("biaj->ijab",self.V_me_s)
        N_biaj_ijab = np.einsum("biaj->ijab",self.N_pqrs_s)

        df_2_ij_A = np.zeros((self.sp_size,self.sp_size))
        df_2_ij_B = np.zeros((self.sp_size,self.sp_size))
        for i in self.sp_list:
            for j in self.sp_list:
                res_A_1 = np.einsum("ab,ab->",self.N_pq_s[0:ferm_1,ferm_2:sp_size],V_biaj_ijab[i][j][0:ferm_1,ferm_2:sp_size])
                res_A_2 = np.einsum("ab,ab->",self.f_me_s[0:ferm_1,ferm_2:sp_size],N_biaj_ijab[i][j][0:ferm_1,ferm_2:sp_size])
                df_2_ij_A[i][j] = res_A_1 - res_A_2
                res_B_1 = np.einsum("ab,ab->",self.N_pq_s[ferm_2:sp_size,0:ferm_1],V_biaj_ijab[i][j][ferm_2:sp_size,0:ferm_1])
                res_B_2 = np.einsum("ab,ab->",self.f_me_s[ferm_2:sp_size,0:ferm_1],N_biaj_ijab[i][j][ferm_2:sp_size,0:ferm_1])
                df_2_ij_B[i][j] = res_B_1 - res_B_2

        df_2 = df_2_ij_A + df_2_ij_B
        #print(df_2_ij_A[df_2_ij_A.nonzero()])
        #quit()

        V_abcj_jcab = np.einsum("abcj->jcab",self.V_me_s)
        V_ciab_icab = np.einsum("ciab->icab",self.V_me_s)
        N_abcj_jcab = np.einsum("abcj->jcab",self.N_pqrs_s)
        N_ciab_icab = np.einsum("ciab->icab",self.N_pqrs_s)
        df_3_1=np.zeros((self.sp_size,self.sp_size))
        #df_3_2=np.zeros((self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                sum_1 = 0
                sum_2 = 0
                for c in self.sp_list:
                    if(c<=self.fermi_index):
                        res_1 = np.einsum("ab,ab",N_ciab_icab[i][c][ferm_2:,ferm_2:],V_abcj_jcab[j][c][ferm_2:,ferm_2:])
                        res_2 = np.einsum("ab,ab",V_ciab_icab[i][c][ferm_2:,ferm_2:],N_abcj_jcab[j][c][ferm_2:,ferm_2:])
                        val= res_1 - res_2
                        sum_1 += val

                    else:
                        res_3 = np.einsum("ab,ab->",N_ciab_icab[i][c][0:ferm_1,0:ferm_1],V_abcj_jcab[j][c][0:ferm_1,0:ferm_1])
                        res_4 = np.einsum("ab,ab->",V_ciab_icab[i][c][0:ferm_1,0:ferm_1],N_abcj_jcab[j][c][0:ferm_1,0:ferm_1])
                        sum_2 += res_3-res_4
                # -----
                df_3_1[i][j] = sum_1+sum_2
                # =====
        #df_3_2 = np.einsum("ij->ji",df_3_1)
        df_3 = 0.5*df_3_1


        df_pq_s = df_1 + df_2 + df_3
        #print(df_3[df_3.nonzero()])
        #print(df_3)
        #sys.exit()
        return df_pq_s

    def cal_G_pqrs_s2(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)

        dG_1_1 = np.einsum("ia,ajkl->ijkl",self.N_pq_s,self.V_me_s)
        dG_1_2 = np.einsum("ia,ajkl->ijkl",self.f_me_s,self.N_pqrs_s)
        dG_1_A = dG_1_1 - dG_1_2
        dG_1_B = np.einsum("ijkl->jikl",dG_1_A)
        dG_1 = dG_1_A - dG_1_B

        dG_2_1 = np.einsum("ak,ijal->ijkl",self.N_pq_s,self.V_me_s)
        dG_2_2 = np.einsum("ak,ijal->ijkl",self.f_me_s,self.N_pqrs_s)
        dG_2_A = dG_2_1 - dG_2_2
        dG_2_B = np.einsum("ijkl->ijlk",dG_2_A)
        dG_2 = -1*(dG_2_A - dG_2_B)

        N_abkl_klab = np.einsum("abkl->klab",self.N_pqrs_s)
        V_abkl_klab = np.einsum("abkl->klab",self.V_me_s)
        dG_3_1=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_2=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_3=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_4=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                for k in self.sp_list:
                    for l in self.sp_list:
                        dG_3_1[i][j][k][l] = np.einsum("ab,ab->",self.N_pqrs_s[i][j][0:ferm_1,0:ferm_1],V_abkl_klab[k][l][0:ferm_1,0:ferm_1])
                        dG_3_2[i][j][k][l] = np.einsum("ab,ab->",self.V_me_s[i][j][0:ferm_1,0:ferm_1],N_abkl_klab[k][l][0:ferm_1,0:ferm_1])
                        dG_3_3[i][j][k][l] = np.einsum("ab,ab->",self.N_pqrs_s[i][j][ferm_2:,ferm_2:],V_abkl_klab[k][l][ferm_2:,ferm_2:])
                        dG_3_4[i][j][k][l] = np.einsum("ab,ab->",self.V_me_s[i][j][ferm_2:,ferm_2:],N_abkl_klab[k][l][ferm_2:,ferm_2:])

        dG_3_A = dG_3_1 - dG_3_2
        dG_3_B = dG_3_3 - dG_3_4

        dG_3 = 0.5*(dG_3_B - dG_3_A)

        N_aibk_ikab = np.einsum("aibk->ikab",self.N_pqrs_s)
        V_bjal_jlab = np.einsum("bjal->jlab",self.V_me_s)

        dG_4_1=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                for k in self.sp_list:
                    for l in self.sp_list:
                        res_1 = np.einsum("ab,ab->",N_aibk_ikab[i][k][0:ferm_1,ferm_2:sp_size],V_bjal_jlab[j][l][0:ferm_1,ferm_2:sp_size])
                        res_2 = np.einsum("ab,ab->",N_aibk_ikab[i][k][ferm_2:sp_size,0:ferm_1],V_bjal_jlab[j][l][ferm_2:sp_size,0:ferm_1])
                        dG_4_1[i][j][k][l] = res_1 - res_2

        dG_4_2 = np.einsum("ijkl->ijlk",dG_4_1)

        dG_4 = 2*(dG_4_1 - dG_4_2)

        dG_pqrs_s = dG_1 + dG_2 + dG_3 #+dG_4

        return dG_pqrs_s

    def cal_f_pq_s3(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size

        df_pq_s = np.zeros((self.sp_size,self.sp_size))
        df_1 = np.zeros((self.sp_size,self.sp_size))
        df_2 = np.zeros((self.sp_size,self.sp_size))
        df_3 = np.zeros((self.sp_size,self.sp_size))


        for i in self.sp_list:
            for j in self.sp_list:
                df_1_t = 0.0
                for a in self.sp_list:
                    df_1_A =  self.N_pq_s[i][a]*self.f_me_s[a][i]
                    df_1_B =  self.f_me_s[i][a]*self.N_pq_s[a][i]
                    df_1_t += df_1_A - df_1_B
                df_1[i][j] = df_1_t

                df_2_t_A = 0.0
                df_2_t_B = 0.0
                for a in self.sp_list:
                    fermi_flag = 0
                    if(a>self.fermi_index):
                        fermi_flag = 1
                    for b in self.sp_list:
                        if(fermi_flag == 0 and b>self.fermi_index):
                            A_1= self.N_pq_s[a][b]*self.V_me_s[b][i][a][j]
                            A_2= self.f_me_s[a][b]*self.N_pqrs_s[b][i][a][j]
                            df_2_t_A += A_1 - A_2
                        if(fermi_flag == 1 and b<=self.fermi_index):
                            B_1= self.N_pq_s[a][b]*self.V_me_s[b][i][a][j]
                            B_2= self.f_me_s[a][b]*self.N_pqrs_s[b][i][a][j]
                            df_2_t_B += B_1 - B_2
                df_2[i][j] = df_2_t_A-df_2_t_B

                df_3_t_A = 0.0
                df_3_t_B = 0.0
                for a in self.sp_list:
                    fermi_flag = 0
                    if(a>self.fermi_index):
                        fermi_flag = 1
                    for b in self.sp_list:
                        if(fermi_flag == 0 and b > self.fermi_index):
                            continue
                        if(fermi_flag == 1 and b <= self.fermi_index):
                            continue
                        for c in self.sp_list:
                            if(fermi_flag == 0 and c > self.fermi_index):
                                A_1 = self.N_pqrs_s[c][i][a][b]*self.V_me_s[a][b][c][j]
                                A_2 = self.V_me_s[c][i][a][b]*self.N_pqrs_s[a][b][c][j]
                                df_3_t_A += A_1 - A_2

                            if(fermi_flag == 1 and c <= self.fermi_index):
                                B_1 = self.N_pqrs_s[c][i][a][b]*self.V_me_s[a][b][c][j]
                                B_2 = self.V_me_s[c][i][a][b]*self.N_pqrs_s[a][b][c][j]
                                df_3_t_B += B_1 - B_2

                df_3[i][j] = 0.5*(df_3_t_A+df_3_t_B)

        df_pq_s = df_1 + df_2 + df_3
        print(df_3[df_3.nonzero()])
        print(df_3)
        sys.exit()
        return df_pq_s


    def cal_G_pqrs_s3(self):
        ferm_1 = self.fermi_index+1
        ferm_2 = self.fermi_index+1
        sp_size = self.sp_size
        print("ferm_1 : ",ferm_1,"ferm_2 : ",ferm_2)

        dG_1_1 = np.einsum("ia,ajkl->ijkl",self.N_pq_s,self.V_me_s)
        dG_1_2 = np.einsum("ia,ajkl->ijkl",self.f_me_s,self.N_pqrs_s)
        dG_1_A = dG_1_1 - dG_1_2
        dG_1_B = np.einsum("ijkl->jikl",dG_1_A)
        dG_1 = dG_1_A - dG_1_B

        dG_2_1 = np.einsum("ak,ijal->ijkl",self.N_pq_s,self.V_me_s)
        dG_2_2 = np.einsum("ak,ijal->ijkl",self.f_me_s,self.N_pqrs_s)
        dG_2_A = dG_2_1 - dG_2_2
        dG_2_B = np.einsum("ijkl->ijlk",dG_2_A)
        dG_2 = -1*(dG_2_A - dG_2_B)

        N_abkl_klab = np.einsum("abkl->klab",self.N_pqrs_s)
        V_abkl_klab = np.einsum("abkl->klab",self.V_me_s)
        dG_3_1=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_2=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_3=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))
        dG_3_4=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                for k in self.sp_list:
                    for l in self.sp_list:
                        dG_3_1[i][j][k][l] = np.einsum("ab,ab->",self.N_pqrs_s[i][j][0:ferm_1][0:ferm_1],V_abkl_klab[k][l][0:ferm_1][0:ferm_1])
                        dG_3_2[i][j][k][l] = np.einsum("ab,ab->",self.V_me_s[i][j][0:ferm_1][0:ferm_1],N_abkl_klab[k][l][0:ferm_1][0:ferm_1])
                        dG_3_3[i][j][k][l] = np.einsum("ab,ab->",self.N_pqrs_s[i][j][ferm_2:][ferm_2:],V_abkl_klab[k][l][ferm_2:][ferm_2:])
                        dG_3_4[i][j][k][l] = np.einsum("ab,ab->",self.V_me_s[i][j][ferm_2:][ferm_2:],N_abkl_klab[k][l][ferm_2:][ferm_2:])

        dG_3_A = dG_3_1 - dG_3_2
        dG_3_B = dG_3_3 - dG_3_4

        dG_3 = 0.5*(dG_3_B - dG_3_A)

        N_aibk_ikab = np.einsum("aibk->ikab",self.N_pqrs_s)
        V_bjal_jlab = np.einsum("bjal->jlab",self.V_me_s)

        dG_4_1=np.zeros((self.sp_size,self.sp_size,self.sp_size,self.sp_size))

        for i in self.sp_list:
            for j in self.sp_list:
                for k in self.sp_list:
                    for l in self.sp_list:
                        res_1 = np.einsum("ab,ab->",N_aibk_ikab[i][k][0:ferm_1][ferm_2:sp_size],V_bjal_jlab[j][l][0:ferm_1][ferm_2:sp_size])
                        res_2 = np.einsum("ab,ab->",N_aibk_ikab[i][k][ferm_2:sp_size][0:ferm_1],V_bjal_jlab[j][l][ferm_2:sp_size][0:ferm_1])
                        dG_4_1[i][j][k][l] = res_1 - res_2

        dG_4_2 = np.einsum("ijkl->ijlk",dG_4_1)

        dG_4 = 2*(dG_4_1 - dG_4_2)

        dG_pqrs_s = dG_1 + dG_2 + dG_3 #+dG_4

        return dG_pqrs_s

    def cal_N_s(self):
        N_pq_s = np.zeros(self.N_pq_s.shape)
        for p in self.sp_list:
            fermi_flag = 0
            if(p>self.fermi_index):
                fermi_flag = 1
            for q in self.sp_list:
                if(fermi_flag==0 and q<=self.fermi_index):
                    continue
                if(fermi_flag==1 and q>self.fermi_index):
                    continue
                delta = self.Delta_pq_s[p][q]
                if(delta == 0):
                    print("Wrong happened here at cal_N_s 1")
                    sys.exit()
                val = self.f_me_s[p][q]/delta
                N_pq_s[p][q] = val
                #print("a : ",a,"\t i : ",i,val)
        N_qp_s = np.einsum("pq->qp",N_pq_s)
        self.N_pq_s = N_pq_s #- N_qp_s
        N_pqrs_s = np.zeros(self.N_pqrs_s.shape)
        for p in self.sp_list:
            fermi_flag = 0
            if(p>self.fermi_index):
                fermi_flag = 1
            for q in self.sp_list:
                if(fermi_flag==1 and q<=self.fermi_index):
                    continue
                if(fermi_flag==0 and q>self.fermi_index):
                    continue
                if(p==q):
                    continue
                for r in self.sp_list:
                    if(fermi_flag==0 and r<=self.fermi_index):
                        continue
                    if(fermi_flag==1 and r>self.fermi_index):
                        continue
                    for s in self.sp_list:
                        if(fermi_flag==0 and s<=self.fermi_index):
                            continue
                        if(fermi_flag==1 and s>self.fermi_index):
                            continue
                        if(r==s):
                            continue
                        delta = self.Delta_pqrs_s[p][q][r][s]
                        if(delta == 0):
                            print("Wrong happened here at cal_N_s 2")
                            sys.exit()
                        val = self.V_me_s[p][q][r][s]/delta
                        N_pqrs_s[p][q][r][s] = val
        N_rspq_s = np.einsum("pqrs->rspq",N_pqrs_s)
        self.N_pqrs_s = N_pqrs_s #- N_rspq_s



    def cal_Delta(self):
        #  self.Delta_pq[a][i]
        Delta_pq_s = np.zeros(self.Delta_pq_s.shape)
        for p in self.sp_list:
            fermi_flag = 0
            if(p>self.fermi_index):
                fermi_flag = 1
            for q in self.sp_list:
                if(fermi_flag==0 and q<=self.fermi_index):
                    continue
                if(fermi_flag==1 and q>self.fermi_index):
                    continue
                val = self.f_me_s[p][p] - self.f_me_s[q][q] #+ self.V_me_s[a][i][a][i]
                Delta_pq_s[p][q] = val
        self.Delta_pq_s=Delta_pq_s
        #  self.Delta_pqrs[a][b][i][j]
        Delta_pqrs_s = np.zeros(self.Delta_pqrs_s.shape)
        for p in self.sp_list:
            fermi_flag = 0
            if(p>self.fermi_index):
                fermi_flag = 1
            for q in self.sp_list:
                if(fermi_flag==1 and q<=self.fermi_index):
                    continue
                if(fermi_flag==0 and q>self.fermi_index):
                    continue
                for r in self.sp_list:
                    if(fermi_flag==0 and r<=self.fermi_index):
                        continue
                    if(fermi_flag==1 and r>self.fermi_index):
                        continue
                    for s in self.sp_list:
                        if(fermi_flag==0 and s<=self.fermi_index):
                            continue
                        if(fermi_flag==1 and s>self.fermi_index):
                            continue
                        val = self.f_me_s[p][p] + self.f_me_s[q][q] - self.f_me_s[r][r] - self.f_me_s[s][s]
                        # val += self.V_me_s[i][j][i][j] + self.V_me_s[a][b][a][b]
                        # val += -1*self.V_me_s[a][i][a][i]-1*self.V_me_s[b][j][b][j]
                        # val += -1*self.V_me_s[a][j][a][j]-1*self.V_me_s[b][i][b][i]
                        Delta_pqrs_s[p][q][r][s] = val
        self.Delta_pqrs_s=Delta_pqrs_s

    def cal_MBPT_2(self):
        res = 0.0
        for i in self.sp_list:
            if(i>self.fermi_index):
                continue
            for j in self.sp_list:
                if(j>self.fermi_index):
                    continue
                for a in self.sp_list:
                    if(a<=self.fermi_index):
                        continue
                    for b in self.sp_list:
                        if(b<=self.fermi_index):
                            continue
                        v2_ijab = self.V_me_s[i][j][a][b]*self.V_me_s[a][b][i][j]
                        delta = self.f_me_s[i][i] + self.f_me_s[j][j] - self.f_me_s[a][a] - self.f_me_s[b][b]
                        if(delta == 0):
                            print("Wrong happened here at cal_MBPT_2 ")
                            sys.exit()
                        res += v2_ijab/delta
        res = res/4.0
        return res


    def build(self):
        ds = 0.1

        self.V_me_s = self.V_me
        self.f_me_s = self.f_me
        self.cal_Delta()
        self.cal_N_s()
        sp_size = self.sp_size
        sp_size_2 = self.sp_size*self.sp_size
        self.srg_conv = np.zeros([100,sp_size_2,sp_size_2])
        self.srg_conv_f = np.zeros([100,sp_size_2,sp_size_2])
        self.srg_conv_G = np.zeros([100,sp_size_2,sp_size_2])
        #self.srg_f_s = np.zeros([100,sp_size,sp_size])
        #self.srg_V_s = np.zeros([100,sp_size,sp_size,sp_size,sp_size])
        self.srg_f_s = []
        self.srg_V_s = []


        dE_d = self.cal_E_s2()
        mbpt_2 = self.cal_MBPT_2()
        print("\t mbpt_2 :",mbpt_2,"\t dE :",dE_d)
        dE_1 = dE_d * 0.1
        #quit()
        time_i = -1
        eps = 0.001
        while (time_i < 100):
            time_i += 1
            index_1 = -1
            h_sum = 0
            self.srg_f_s.append(self.f_me_s)
            self.srg_V_s.append(self.V_me_s)

            for a in self.sp_list:
                for b in self.sp_list:
                    if (a==b):
                        continue
                    index_1 += 1
                    index_2 = -1
                    for c in self.sp_list:
                        for d in self.sp_list:
                            if (c==d):
                                continue

                            index_2 += 1


                            val_f = self.f_me_s[a][c] + self.f_me_s[b][d]
                            val_v = self.V_me_s[a][b][c][d]
                            val =  val_v + val_f
                            #print(time_i,index_1,index_2)
                            self.srg_conv[time_i][index_1][index_2] = abs(val)
                            self.srg_conv_f[time_i][index_1][index_2] = abs(val_f)
                            self.srg_conv_G[time_i][index_1][index_2] = abs(val_v)

                            if(a!=c and b!=d):
                                h_sum += abs(val_v)

            #self.srg_conv[time_i] = (self.f_me_s + self.V_me_s) * (self.f_me_s + self.V_me_s)
            df_s = self.cal_f_pq_s2()
            dG_s = self.cal_G_pqrs_s2()

            self.f_me_s = self.f_me_s + df_s*ds
            self.V_me_s = self.V_me_s + dG_s*ds
            V_diff_2 = self.V_me_s*self.V_me_s
            V_me_s_diff = np.einsum("abcd->",V_diff_2) - np.einsum("ijij->",V_diff_2)

            self.cal_Delta()
            self.cal_N_s()
            dE_d = self.cal_E_s2()
            dE_1 += dE_d*ds
            mbpt_2 = self.cal_MBPT_2()
            print(time_i,"\t mbpt_2 :",mbpt_2,"\t dE_d : ",dE_d,"\t dE_1  :",dE_1 ,h_sum)
            if(abs(dE_d) < eps or V_me_s_diff < eps):
                self.dE = dE_1
                break






