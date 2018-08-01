#!/usr/bin/python
from basis_k import *

class Minn_ME:
    # 'Minnesota potential\begin{align}
    # \langle k_p k_q \vert V_\alpha\vert k_r k_s \rangle = {V_\alpha\over L^3} \left({\pi\over\kappa_\alpha}\right)^{3/2}
    # e^{- {q^2 \over 4\kappa_\alpha}} \delta_{k_p+k_q}^{k_r+k_s} .
    # \label{_auto40}
    # \end{align}'
    # \begin{align}
    # \langle	k_p s_p k_q s_q\vert V_\alpha\vert k_r s_r k_s s_s\rangle &= \langle	k_p k_q \vert V_\alpha\vert k_r k_s \rangle
    # \left(\delta_{s_p}^{s_r}\delta_{s_q}^{s_s} - \delta_{s_p}^{s_s}\delta_{s_q}^{s_r}\right) ,
    # \label{_auto41}
    # \end{align}
    def __init__(self,sps_k,tb_k):
        self.sps_k = sps_k
        self.tb_k = tb_k
        self.k_step = sps_k.k_step
        self.L_k3 = sps_k.L_k*sps_k.L_k*sps_k.L_k

        self.hbarc =  197.3269678792965;    # Mev.fm
        self.M_n = 938.918725 # Mev/c2
        self.hc2_2M = self.hbarc**2/(2.0*self.M_n)

        self.V_R = 200 #MeV
        self.V_S = -91.85 #MeV
        self.V_T = -178 #MeV
        self.V_RST = np.array([self.V_R,self.V_S,self.V_T])

        self.Kp_R = 1.487 #fm-2
        self.Kp_S = 0.465 #fm-2
        self.Kp_T = 0.639 #fm-2
        self.Kp_RST = np.array([self.Kp_R,self.Kp_S,self.Kp_T])

        self.pre_factor_1 = self.V_RST/self.L_k3
        self.pre_factor_2_2 = np.pi/self.Kp_RST
        self.pre_factor_2 = np.power(self.pre_factor_2_2,1.5)
        self.pre_factor = self.pre_factor_1*self.pre_factor_2

    def cal_V_sub(self,alpha,kpq,krs):
        # kpq,krs are two body state
        # kpq.k12 k1 - k2
        res = 0.0
        if(kpq.K_tot[0] != krs.K_tot[0]):
            return res
        if(kpq.K_tot[1] != krs.K_tot[1]):
            return res
        if(kpq.K_tot[2] != krs.K_tot[2]):
            return res

        q_trans = 0.5*self.k_step*(kpq.k12 - krs.k12)
        #print(q_trans)
        #q_trans_2 = q_trans[0]**2 + q_trans[1]**2 + q_trans[2]**2
        q_trans_2 = np.dot(q_trans,q_trans)
        #print(q_trans_2)

        #print("q_trans_2 : ",q_trans_2)
        x_alpha = -1*q_trans_2/(4.0*self.Kp_RST[alpha])
        #print("x_alpha : ",x_alpha)

        res = self.pre_factor[alpha]*np.exp(x_alpha)
        return res

    def cal_V_neu(self,kpq,krs):
        # kpq,krs are two body state
        res = 0.0
        delta_1 = 0
        delta_2 = 0
        # if(kpq.K_tot[] != krs.K_tot):
        #     return res
        #print("\tS12p:",kpq.S12_2,"\t S12:",krs.S12_2)
        if(kpq.S12_2 != krs.S12_2):
            return res
        # ==== !!! becare only neutron matter !!!====
        if(kpq.S12_2 != 0 ):
            return res
        #print(" ---- S12p:",kpq.S12_2,"\t S12:",krs.S12_2)
        # ==== !!! becare only neutron matter !!!====
        if((kpq.s2_1 == krs.s2_1) and (kpq.s2_2 == krs.s2_2)):
            delta_1 = 1
        if((kpq.s2_1 == krs.s2_2) and (kpq.s2_2 == krs.s2_1)):
            delta_2 = 1
        delta = delta_1 - delta_2
        # if(delta != 10):
        #     print(" ---- s_p:",kpq.s2_1,"\t s_q:",kpq.s2_2,"\t s_r:",krs.s2_1,"\t s_s:",krs.s2_2)
        #     print(" ---- delta_1 = ", delta_1, "\t delta_2 = ",delta_2,"\t delta : ",delta)

        if(delta == 0.0):
            return res
        res_0 = self.cal_V_sub(0,kpq,krs)
        res_1 = self.cal_V_sub(1,kpq,krs)
        res = 0.5*delta*(res_0 + res_1)
        #print("res_0 : ", res_0,"res_2")
        return res

    def cal_V_neu_as(self,kpq,krs):
        #ksr = krs
        index_1 = krs.index_2
        index_2 = krs.index_1
        s2_1 = krs.s2_2
        s2_2 = krs.s2_1
        K_tot = krs.K_tot
        k12 = -1*krs.k12
        S12_2 = krs.S12_2
        #self.index = -1
        V_pqrs = self.cal_V_neu(kpq,krs)
        ksr = TB_k(index_1,index_2,K_tot,k12,s2_1,s2_2,S12_2)
        #print("V_pqrs : ",V_pqrs)
        V_pqsr = self.cal_V_neu(kpq,ksr)

        V = V_pqrs - V_pqsr
        #print(V_pqrs,V_pqsr)
        return V



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

        k2 = kp_xyz[0]**2 + kp_xyz[1]**2 + kp_xyz[2]**2
        res = self.hc2_2M*k2
        return res











