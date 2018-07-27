#!/usr/bin/python
import numpy as np
import itertools
import sys
import matplotlib.pyplot as plt
from basis_k import *
#from ccsd_me import *



particle_num = 14

N_max = 1
sps=Basis_SP_K(N_max,particle_num)
sps.build()
#sps.print_state()
#sps.print_state()

tb_k = Basis_TB_k(sps)
tb_k.build()

# ==============++++++++++==============#


# ==============++++++++++==============#





