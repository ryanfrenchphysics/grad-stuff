import math
import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as inte

def dispersion(e_0, k, a, t):
	#e_0 : ground state energy?
	#k : wave number
	#a : lattice spacing 
	#t : nearest neighbor hopping integral value

	#dispersion relation for 1D tight binding

	ret = e_0 + 2.*t*math.cos(k*a)
	return ret

def gaussian(e, e_k, f):
	#e : energy
	#e_k : dispersion relation
	#f : full width half max of gaussian

	#gaussian used to approximate delta function

	ret = (2./f)*math.sqrt(math.log(2)/math.pi)*math.exp(-(4.*math.log(2.)*(e-e_k)**2.)/f**2.)
	return ret

def main(a, e_0, t, f, num_E, num_k):
	#a : lattice spacing
	#num_k : number of k_steps to go between -pi/a to pi/a
	#num_E : number of energy steps	

	E_k = []
	k_list = []
	k = -math.pi/a
	dk = 2.*math.pi/(a*num_k)
	for i in xrange(0,num_k + 1):
		k_list.append(k)
		E_k.append(dispersion(e_0, k, a, t))
		k = k + dk
	
	E_max = max(E_k)
	dE = E_max/num_E
	E = 0.
	E_list = []
	D_list = []
	for i in xrange(0,num_E):
		func = lambda k : gaussian(E,dispersion(e_0,k,a,t),f)
		D_list.append(E_max*inte.quad(func,-math.pi/a,math.pi/a)[0])
		E_list.append(E/E_max)
		E = E + dE
	print D_list[0]
	plt.plot(E_list,D_list)
	plt.show()

main(1,0.2,0.1,.01,1000,1000) #e_0 should be twice t

