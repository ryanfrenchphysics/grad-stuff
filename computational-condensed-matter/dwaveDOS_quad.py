import math
import scipy.integrate as inte
import scipy.special as sci
from matplotlib import pyplot as plt

def DOS(stepsize, start, finish, epsilon):
	n = int((finish-start)/stepsize)
	X = []
	N = []
	x = start
	for i in xrange(0,n):
		if abs(x) < 1. and x != 0.:
				X.append(x)
				func = lambda phi: x/math.sqrt(x**2. - math.cos(2.*phi)**2.)
				N.append(inte.quad(func,0.5*math.acos(x)+epsilon,0.5*(math.pi-math.acos(x))-epsilon)[0])
		elif abs(x) > 1:
			X.append(x)
			#N.append(inte.quad(func,0,2.*math.pi)[0])
			N.append((2./math.pi)*sci.ellipk(1./abs(x)))
		x = x + stepsize
	plt.plot(X,N)
	plt.show()
	
DOS(.01,0.,2.,0.001)
