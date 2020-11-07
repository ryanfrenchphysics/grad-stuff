import numpy as np
from scipy.optimize import fsolve
from numba import jit#, types
from numba.extending import overload#, register_jitable
# from numba.errors import TypingError


def func(x):
    return x + 2*np.cos(x)


'''
Used as:
np.absolute(2. * np.pi * fsolve(FullG, [0.9], args=T))[0]

FullG: Float64
[0.9]: ndarray
args: Tuple (T: Float64)
'''
@overload(fsolve)
def jit_fsolve(fxn, x0, args=None):
    @jit(nopython=True, fastmath=True)
    def use(f, x, _a):
        if args is None:
            return fsolve(fxn, x0)
        else:
            return fsolve(fxn, x0, args)

def main():
    # These take the same time
    a = fsolve(func, 0.2)
    #b = jit_fsolve(func, 0.2)
