from itertools import permutations
from math import cos
import numpy as np
cimport numpy as np



cpdef double dispersion(np.ndarray[double, ndim=1] k_perms, int dimensions, int num_k, double hop, double spacing, double e_0):
    cdef double val = e_0

    for i in range(dimensions):
        val -= 2 * hop * cos(k_perms[i] * spacing)
    
    return val

cpdef np.ndarray[double, ndim=2] k_perms(np.ndarray[double, ndim=2] k_mat, int dimensions, int num_k):
    k_vals = np.zeros((num_k))

    for i in range(num_k):
        k_vals[i] = k_mat[i][0]

    k_permutations = list(permutations(k_vals, dimensions))
    for i in range(num_k):
        k_permutations.append(tuple(k_mat[i]))

    k_permutations = np.asarray(k_permutations)

    return k_permutations

cpdef np.ndarray[double, ndim=1] gen_E(double start, double final, int n):
    cdef double step = (final - start) / n

    Es = np.zeros((n + 1))
    for i in range(n + 1):
        Es[i] = start + (i * step)
    
    return Es
