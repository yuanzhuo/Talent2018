#!/usr/bin/python
from basis_k import *

class minn_me:
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
        self.R_const = 200 #MeV
        self.S_const = -91.85 #MeV
        self.T_const = -178 #MeV
