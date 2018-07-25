from sympy import *
from pylab import *
import matplotlib.pyplot as pltm
import numpy as np
import copy
from bit_configs import *
g = Symbol('g')
f = Symbol('f')
V = Symbol('V')
t = Symbol('t')

def nucleus(part_no, part_states, d, truncate):
    class single_particle_states:
        def __init__(self, total_number, total_levels, energy_spacing_d):
            self.total_number=total_number
            self.total_levels=total_levels
            self.occupations=2
            self.levels=[]
            part_no=-1
            fill=1

            fermions=self.total_number
            this_occupation=self.occupations
            for level in np.arange(0, self.total_levels, 1):
                for occupation in np.arange(0, self.total_levels, 1):
                    if this_occupation > 0:
                        if fermions > 0:
                            fermions-=1
                            part_no+=1
                        else:
                            part_no+=1
                            fill=0
                        if part_no==-1:
                            self.levels.append([fill, part_no, level, 0, level*energy_spacing_d])
                        elif this_occupation%2:
                            self.levels.append([fill, part_no, level, -1, level*energy_spacing_d])
                        else:
                            self.levels.append([fill, part_no, level, 1, level*energy_spacing_d])
                        this_occupation -= 1
                    else:
                        break

                this_occupation=self.occupations

    ground_string='1'*part_no+'0'*(part_states-part_no)
    single_part_configs=list(perm_unique(ground_string))
    single_part_configs=[''.join(config) for config in single_part_configs]
    no_single_configs=len(single_part_configs)
    print(no_single_configs,"single_part_configs")
    pair_configs=[config for config in single_part_configs if config.count('11')==part_no/2 and config.count('00')==(part_states-part_no)/2]
    no_pairs_configs=len(pair_configs)
    print(no_pairs_configs, "pair_configs")
    sp_0=single_particle_states(part_no, part_states/2, d)
    hole_states=[level[1] for level in sp_0.levels if level[0]] #ijk..
    particle_states=[level[1] for level in sp_0.levels if not level[0]] #abc...

    def count_excitation(config):
        r_h=0
        for b_i, b1 in enumerate(config[:part_no]):
            if b1 != ground_string[b_i]:
                r_h+=1
        return r_h

    # def P(a,b,i,j):
    #
    #     return 0

    def make_slater_dets(configs, part_no, part_states, d):
        #returns dictionary of slater determinants for the basis
        config_slater_det = {}
        orders=[]
        config_total_energies=[]
        for slater_no, config in enumerate(configs):
            sp_c = copy.deepcopy(sp_0)
            config_energies=[]
            for part_i, fill in enumerate(config):
                sp_c.levels[part_i] = [int(fill), sp_c.levels[part_i][1], sp_c.levels[part_i][2], sp_c.levels[part_i][3], sp_c.levels[part_i][4]]
                if int(fill): config_energies.append(sp_c.levels[part_i][4])
            order=count_excitation(config)
            orders.append(order)

            config_total_energy = sum(config_energies)
            config_total_energies.append((config_total_energy, config))
            config_slater_det[slater_no] = sp_c

        pair_configs=np.array(sorted(config_total_energies))[:,1]

        return config_slater_det, pair_configs

    Ham=Matrix(np.ones((no_pairs_configs, no_pairs_configs))*g)
    pairs_slater_dets, pair_configs  = make_slater_dets(pair_configs, part_no, part_states, d)
    no_pairs_configs=len(pair_configs)

    def basic_ham(c1, c2): # not using state information, just counts pairs..
        match, h=0, 0
        for lev_i, level1 in enumerate(c1):
            if level1 == c2[lev_i] and int(level1) == 1:
                match += 1
        if match >= part_no:
            config_energies=[]
            for part_i, fill in enumerate(c1):
                if int(fill):
                    config_energies.append(sp_0.levels[part_i][4])
            h+=sum(config_energies)

        h-=2.0*g*float(match/4.0)

        return h


    def V_int(a,b,i,j):
        #lookup table
        return V

    def P(p,q):
        if p in part_states:
            parity_states=part_states #p, q = a, b
        if p in hole_states:
            parity_states=hole_states #p, q = i, j
        if parity_states.index(p) < parity_states:
            return 1
        if parity_states.index(p) > parity_states:
            return -1

    #make fock_matrix
    fock_matrix=Matrix(np.ones((part_states, part_states))*V)
    for p in np.arange(0, part_states, 1):
        for q in np.arange(0, part_states, 1):
            val=0
            if p==q:
                val+=sp_0.levels[p][4]
            for i in hole_states:
                val+=V_int(p,i,q,i)
            fock_matrix[p,q] = val


    def bar_H_ij_ab(i,j,a,b):
        h=0
        h+=V(a,b,i,j)
        h+=P(a,b)*1

        #fock_matrix[p,q]

        return h

    for row in range(no_pairs_configs):
        for col in range(no_pairs_configs):
            c1=pair_configs[row]
            c2=pair_configs[col]
            c1_ex=count_excitation(c1)
            c2_ex=count_excitation(c2)

            Ham[row,col]=basic_ham(c1, c2) #HF diagonalization hamiltonian
            # fock_matrix[row,col]=make_fock_matrix(c1, c2)

    # fock_matrix()
    # print(np.matrix(Ham))
    print("fock_matrix", np.matrix(fock_matrix))

    for i in hole_states: #"in"
        for a in particle_states: #"out"
            print(i, a)
            # bar_H_ij_ab()

    return Ham, fock_matrix

nucleus(4,8,1,-1)

# Ham, fock_matrix=nucleus(4,8,1,-1)

# print(fock_matrix)

# fock_matrix=np.diag(np.diag(np.matrix(Ham)))

# print(np.diagnp.matrix(Ham))

# print(fock_matrix)

# # for no_particles, no_states, truncate in [[4, 10, 0], [4, 10, 0], [6, 10, 0], [8, 10, 0]]:
# # for no_particles, no_states, truncate in [[6, 4*2, -1], [6, 8*2, -1], [6, 12*2, -1]]:
# # for no_particles, no_states, truncate in [[4, 4*2, -1], [4, 8*2, -1], [4, 12*2, -1]]:
# for no_particles, no_states, truncate in [[8, 8*2, -1], [8, 10*2, -1]]:#, [8, 12*2, -1]
#
#     Ham = nucleus(no_particles, no_states, 1, truncate)
#
#     print(np.matrix(Ham))
#
#     e1, e2 = [],[]
#     gs = linspace(-1, 1, 20)
#
#     for g_val in gs:
#         print("g", g_val)
#         Ham_time=np.matrix(matrix(Ham.subs({g: g_val*1.0})), dtype=float)
#
#         if g_val == gs[0] and 0:
#             plt.matshow(np.matrix(Ham_time))
#             plt.show()
#
#         u1, v1 = linalg.eig(Ham_time)
#         e1.append(min(u1))
#         Ham_time_0=np.diag(np.diag(Ham_time))
#         u1, v1 = linalg.eig(Ham_time_0)
#         e2.append(min(u1))
#
#     corr=np.array(e1)-np.array(e2)
#
#     exact = plt.plot(gs, corr,'-o',linewidth = 2.0, label = str(no_particles) + 'particles'+ str(no_states) +'states'+str(truncate)+"truncate")
#
# # plt.axis([-1.1, 1.1, -7, 1])
# # plt.axis([-1.1, 1.1, -11, 2])
# plt.legend()
# plt.xlabel("g")
# plt.ylabel('Correlation energy, E')
# plt.show()
#
# #CC theory lecture notes in phys ch. 8
# #could add 2nd order PT by looping with app. int. sum
