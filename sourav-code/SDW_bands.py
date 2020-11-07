import math
import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
import scipy.integrate as inte
from mpl_toolkits.mplot3d import Axes3D

def dispersion(mu,k_x,k_y,a,t1,t2,t3):
	#mu : chemical potential
	#k_x : x - momentum (p = hbar*k)
	#k_y : y - momentum
	#a : lattice parameter
	#t1 : nearest neighbor hopping integral
	#t2 : next nearest neighbor hopping
	#t3 : next next nearest neighbor hopping

	ret = mu - t1*(np.cos(k_x*a)+np.cos(k_y*a)) - t2*np.cos(k_x*a)*np.cos(k_y*a) - t3*(np.cos(2.*k_x*a)+np.cos(2.*k_y*a))
	return ret

def SC_dispersion(delta,e_k):
	#delta : super conducting order parameter
	#e_k = : dispersion relation in normal state

	ret_plus = np.sqrt(delta**2 + e_k**2)
	ret_minus = -np.sqrt(delta**2 + e_k**2)
	return ret_plus, ret_minus


def SC_dispersion_D(delta,k_x, k_y,e_k):
	#delta : super conducting order parameter
	#e_k = : dispersion relation in normal state

	delta_k = delta * (np.cos(k_x) - np.cos(k_y)) / 2
	ret_plus = np.sqrt(delta_k**2 + e_k**2)
	ret_minus = -np.sqrt(delta_k**2 + e_k**2)
	return ret_plus, ret_minus


def SDW_dispersion_comp_OAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3):
	if (k_x >= 0. and k_y >= 0.) or (k_x <= 0. and k_y <= 0.):
		m = M
	else:
		m = 0.
	E = dispersion(mu,k_x,k_y,a,t1,t2,t3)
	E_q = dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3)
	H = np.array([[E,m],[m,E_q]])
	evals, evecs = la.eig(H)
	return np.sort(evals)

def SDW_dispersion(M,Q,mu,k_x,k_y,a,t1,t2,t3):
	#M : Spin Density Wave order parameter
	#Q : nesting vector

	X = np.array(k_x)
	Y = np.array(k_y)

	ix = int(np.sqrt(X.size))
	iy = int(np.sqrt(Y.size))

	ret = np.zeros((ix,iy))

	for i in range(0,ix):
		for j in range(0,iy):

			if (X[i][j] > 0. and Y[i][j] > 0.) or (X[i][j] < 0. and Y[i][j] < 0):
				xi_plus = 0.5*(dispersion(mu,X[i][j],Y[i][j],a,t1,t2,t3) + dispersion(mu,X[i][j]+Q[0],Y[i][j]+Q[1],a,t1,t2,t3))

				xi_minus = 0.5*(dispersion(mu,X[i][j],Y[i][j],a,t1,t2,t3) - dispersion(mu,X[i][j]+Q[0],Y[i][j]+Q[1],a,t1,t2,t3))
				ret[i][j] = xi_plus + np.sign(xi_minus)*np.sqrt(M**2. + xi_minus**2.)
	#ret_minus = xi_plus - np.sqrt(M**2. + xi_minus**2.)
			else:
				ret[i][j] = dispersion(mu,X[i][j],Y[i][j],a,t1,t2,t3)
			#print X[i][j], Y[i][j], ret[i][j]
	return ret

def SDW_dispersionOAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3):
	#M : Spin Density Wave order parameter
	#Q : nesting vector

	"""
	Since I couldn't have if-statements if k_x & k_y were numpy arrays, I had this function
	just make its own k-grid and calculate the dispersion relation array, but this program
	requires the dispersion One element At A Time (OOAT)
	"""

	"""
	X = np.array(k_x)
	Y = np.array(k_y)

	ix = int(np.sqrt(X.size))
	iy = int(np.sqrt(Y.size))

	ret = np.zeros((ix,iy))

	for i in range(0,ix):
		for j in range(0,iy):
	"""
	if (k_x > 0. and k_y > 0.) or (k_x < 0. and k_y < 0):
		xi_plus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) + dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))

		xi_minus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) - dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))
		ret = xi_plus + np.sign(xi_minus)*np.sqrt(M**2. + xi_minus**2.)
	#ret_minus = xi_plus - np.sqrt(M**2. + xi_minus**2.)
	else:
		ret = dispersion(mu,k_x,k_y,a,t1,t2,t3)
			#print X[i][j], Y[i][j], ret[i][j]
	return ret

def SDW_dispersion_RBZ(M,Q,mu,k_x,k_y,a,t1,t2,t3):
	#M : Spin Density Wave order parameter
	#Q : nesting vector

	xi_plus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) + dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))

	xi_minus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) - dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))
	ret_plus = xi_plus + np.sqrt(M**2. + xi_minus**2.)
	ret_minus = xi_plus - np.sqrt(M**2. + xi_minus**2.)
	return ret_plus, ret_minus

def fermi_surface(mu,a,t1,t2,t3,num_k):
	#see above for what these parameters mean
	#num_k : number of points in 1d in k-space

	k_x = np.linspace(-np.pi/a,np.pi/a,num_k)
	k_y = np.linspace(-np.pi/a,np.pi/a,num_k)
	X,Y = np.meshgrid(k_x,k_y)
	function = dispersion(mu,X,Y,a,t1,t2,t3) - mu
	X_new = X/np.pi
	Y_new = Y/np.pi
	plt.contour(X_new,Y_new,function,[0])
	plt.show()


def band_structure3D(mu,a,t1,t2,t3,delta,num_k):
	x = np.linspace(0.,np.pi/a,num_k)
	y = np.linspace(0.,np.pi/a,num_k)
	X,Y = np.meshgrid(x,y)

	# S-wave:
	e_k_plus, e_k_minus = SC_dispersion(delta,dispersion(mu,X,Y,a,t1,t2,t3))

	# D-wave:
	#e_k_plus, e_k_minus = SC_dispersion_D(delta,X, Y,dispersion(mu,X,Y,a,t1,t2,t3))

	e_k = dispersion(mu,X,Y,a,t1,t2,t3)


	e_k_SDW = SDW_dispersion(delta,[np.pi,np.pi],mu,X,Y,a,t1,t2,t3)
	e_k_SDW_2 = SDW_dispersion(delta,[np.pi,0.],mu,X,Y,a,t1,t2,t3)

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1,projection='3d')
	#surf = ax.plot_surface(X,Y,e_k)
	surf2 = ax.plot_surface(X,Y,e_k_plus)
	surf3 = ax.plot_surface(X,Y,e_k_minus)
	ax.contour(X, Y, e_k, 0, colors = ['black'])
	#surf4 = ax.scatter(X,Y,e_k_SDW, color = 'b')
	#surf5 = ax.scatter(X,Y,e_k_SDW_2,color = 'r')
	plt.show()

band_structure3D(0.,1.,1.,0.6,0.,0.3,150)
