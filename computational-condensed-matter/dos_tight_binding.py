import numpy as np
import math
from itertools import permutations
from matplotlib import pyplot as plt
from dos import *

# GLOBAL VALS
DIM = 1
ERROR = 0.1
NUM_K = 1000
NUM_E = 1000
HOP = 1
E_0 = 2 * DIM
SPACING = 0.01

# NONSENSE COMMENT
# # Dispersion E(k) fxn
# def dispersion(k_perms, dimensions, num_k, t, spacing, e_0):
#     # Starting E
#     val = e_0

#     for i in range(dimensions):
#         val -= 2 * t * math.cos(k_perms[i] * spacing)

#     return val



# # Function for k matrix -> k permutation matrix
# def k_perms(k_mat, dimensions, num_k):
#     k_vals = []
#     for i in range(num_k):
#         k_vals.append(k_mat[i][0])

#     # Permutations (as list) of list of k_vals, choose dimensions
#     k_perms = list(permutations(k_vals, dimensions))
    
#     # We must add list of k_vals which have all same k vals
#     k_append = k_perms.append
#     for i in range(num_k):
#         k_append(tuple(k_mat[i])) 

#     #print(k_perms)
#     return k_perms


# def gen_E(start, final, n):
#     step = (final - start) / n

#     Es = []
#     for i in range(n + 1):
#         Es.append(start + (i * step))

#     return Es


def main():
    k_mat = np.zeros((NUM_K, DIM))
    stepsize = ((2 * math.pi) / SPACING) / NUM_K
    for i in range(NUM_K):
        k_mat[i] = (-math.pi / SPACING) + (i * stepsize)

    k_perm_mat = k_perms(k_mat, DIM, NUM_K)

    # Generate list of dispersion, E(k), vals
    Dispersions = []
    for i in range(len(k_perm_mat)):
        Dispersions.append(dispersion(k_perm_mat[i], DIM, NUM_K, HOP, SPACING, E_0))
    Dispersions = np.asarray(Dispersions)


    Es = gen_E(0, 4, NUM_E)

    E_max = max(Dispersions)

    D = np.zeros((NUM_E + 1))
    for i in range(len(Es)):
        states_at_E = 0
        for disp in Dispersions:
            if abs(Es[i] - disp) <= ERROR:
                states_at_E += 1

        D[i] = (states_at_E * E_max)
   
    E_norm = np.zeros((NUM_E + 1))
    for i in range(len(Es)):
        Es[i] = Es[i] / E_max
        E_norm[i] = Es[i]
        
    plt.plot(E_norm, D)
    plt.show()

if __name__ == '__main__':
    main()
