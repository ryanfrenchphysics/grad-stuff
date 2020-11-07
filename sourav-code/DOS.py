import math
import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as inte

hbar = 1.
m_e = 1.

def dispersion(mu,k_x,k_y,a,t1,t2,t3):
	#mu : chemical potential
	#k_x : x - momentum (p = hbar*k)
	#k_y : y - momentum
	#a : lattice parameter
	#t1 : nearest neighbor hopping integral
	#t2 : next nearest neighbor hopping
	#t3 : next next nearest neighbor hopping

	ret = mu + 0.5*t1*(np.cos(k_x*a)+np.cos(k_y*a)) + t2*np.cos(k_x*a)*np.cos(k_y*a) + 0.5*t3*(np.cos(2.*k_x*a)+np.cos(2.*k_y*a))
	#ret=  mu-t1*(np.cos(k_x*a)+np.cos(k_y*a))-t2*(np.cos(k_x*a)*np.cos(k_y*a))
	return np.absolute(ret)

def SDW_dispersion(M,Q,mu,k_x,k_y,a,t1,t2,t3):
    #M : Spin Density Wave order parameter
    #Q : nesting vector


    xi_plus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) + dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))

    xi_minus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) - dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))
    ret_a = xi_plus + np.sqrt(M**2. + xi_minus**2.)
    ret_b = xi_plus - np.sqrt(M**2. + xi_minus**2.)
    return np.abs(ret_a), np.abs(ret_b)

def integrator_delta_herm(integrate,dispersion,dk,E,n,f):
	#integrate : 2D array of the function we wish to integrate over k-space
	#dispersion : 2D array with all dispersion relation values over k-space
	#dk : distance between k values in our grid
	#E : energy that we're currently integrating at
	#n : number of terms in hermite polynomial expansion

	coeff = np.zeros(n+1)
	for i in range(0,len(coeff)):
	    if i % 2 == 0:
	        j=float(i/2)
	        coeff[i] = ((-1)**float(j))/(math.factorial(float(j))*4**float(i)*np.sqrt(np.pi))
	    else:
	        coeff[i] = 0.
	dispersion_array = np.array(dispersion)
	integrate_array = np.array(integrate)

	##the variuos choices of the gaussian##
	exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(2.*f**2.))))
	#exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(4.*f))))
	hermite = np.array(np.polynomial.hermite.hermval(E-dispersion_array,coeff))

	##Multiplying the gaussian with the hermite polynomial expansion##
	func = (1./np.sqrt(2*np.pi*f))*exp*hermite
	#func = exp*hermite*(1./(2.*np.sqrt(np.pi*f)))

	ret = func*integrate_array
	#print integrate_array, func, ret
	#print np.sum(ret)
	return np.sum(ret)*dk*dk

def DOS(M, Q, mu,a,t1,t2,t3,num_E,E_min,E_max,num_k,n,f):
	#see above for already defined parameters
	#num_E : number of energy steps
	#E_min : lower bound on energy
	#E_max : upper bound on energy

	dE = (E_max - E_min)/num_E
	E = E_min
	E_list = []
	D_list = []
	D_lista = []
	D_listb = []

	k_x = np.linspace(-np.pi/a,np.pi/a,num_k)
	k_y = np.linspace(-np.pi/a,np.pi/a,num_k)
	X,Y = np.meshgrid(k_x,k_y)
	function = dispersion(mu,X,Y,a,t1,t2,t3)
	functiona, functionb = SDW_dispersion(M, Q, mu,X,Y,a,t1,t2,t3)
	dk = 2.*np.pi/(a*num_k)


	for i in range(0,num_E):
		#func = lambda k_x, k_y : gaussian(E,dispersion(mu,k_x,k_y,a,t1,t2,t3),f)
		#D_list.append(((a/np.pi)**2.)*inte.dblquad(func,-np.pi/a,np.pi/a,lambda k_x: -np.pi/a,lambda k_x: np.pi/a)[0])
		D_list.append(((a/np.pi)**2.)*integrator_delta_herm(1.,function,dk,E,n,f))
		D_lista.append(((a/np.pi)**2.)*integrator_delta_herm(1.,functiona,dk,E,n,f))
		D_listb.append(((a/np.pi)**2.)*integrator_delta_herm(1.,functionb,dk,E,n,f))
		E_list.append(E)
		E = E + dE

	plt.plot(E_list,D_list, color='r')
	plt.plot(E_list,D_lista, color='b')
	plt.plot(E_list,D_listb, color='g')
	plt.show()

#DOS(mu,a,   t1,t2,t3,num_E, E_min,E_max,num_k,n,f))
DOS(0.1, [np.pi, np.pi],0.29, 1.0,-1.0,0.4,-0.1,500,-0.5,2.,300,10,.1)
