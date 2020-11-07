import scipy.special as sci
import scipy.integrate as inte
from matplotlib import pyplot as plt
import math

def DOS(stepsize, start, finish):
	n = int((finish-start)/stepsize)
	X = []
	N = []
	x = start
	for i in xrange(0,n):
		X.append(x)
		if abs(x) < 1:
			N.append((2./math.pi)*abs(x)*sci.ellipk(abs(x)))
		elif abs(x) > 1:
			N.append((2./math.pi)*sci.ellipk(1./abs(x)))
		x = x + stepsize
	plt.plot(X,N)
	for i in xrange(0,len(X)):
		print X[i], N[i]
	plt.show()

def DOS2(stepsize, start, finish, epsilon):
	n = int((finish-start)/stepsize)
	X = []
	N = []
	x = start
	for i in xrange(0,n):
		X.append(x)
		if abs(x) < 1 and x != 0:
			print x
			func = lambda phi: x/math.sqrt(x**2. - math.cos(2.*phi)**2.)
			N.append(inte.quad(func,0.5*(math.pi-math.acos(x))-epsilon,0.5*math.acos(x)+epsilon)[0]/(2.*math.pi))
		elif abs(x) > 1:
			N.append((2./math.pi)*sci.ellipk(1./abs(x)))
		x = x + stepsize
	plt.plot(X,N)
	plt.show()

DOS(.01,0.,2.)	
