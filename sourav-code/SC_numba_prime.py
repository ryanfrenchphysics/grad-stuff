#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')


# In[2]:


import math
import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
import scipy.integrate as inte
from scipy.optimize import fsolve
from numba import jit

#hbar = 1.05e-34
#m_e = 9.11e-31
#k_b = 1.38e-23

hbar = 1.
m_e = 1.
k_b = 1.

#Temperature = 0.1


# In[3]:


@jit
def dispersion(mu,k_x,k_y,a,t1,t2,t3):
    #This is the dispersion relation for the band.
    #(contains both tight binding as well as quadratic dispersion relations)
    #mu : chemical potential
    #k_x : x - momentum (p = hbar*k)
    #k_y : y - momentum
    #a : lattice parameter
    #t1 : nearest neighbor hopping integral
    #t2 : next nearest neighbor hopping
    #t3 : next next nearest neighbor hopping
    ret = mu + 0.5*t1*(np.cos(k_x*a)+np.cos(k_y*a)) + t2*np.cos(k_x*a)*np.cos(k_y*a) + 0.5*t3*(np.cos(2.*k_x*a)+np.cos(2.*k_y*a))
    #ret = mu - t1*(np.cos(k_x*a)+np.cos(k_y*a)) - t2*np.cos(k_x*a)*np.cos(k_y*a)
    #ret = 0.5*(np.power(k_x,2.) + np.power(k_y,2.))
    #return np.absolute(ret)
    return ret


# In[4]:


@jit
def SC_dispersion(delta,e_k):
    #This is the analytic expression for the superconductor dispersion relation
    #delta : super conducting order parameter
    #e_k = : dispersion relation in normal state
    ret_plus = np.sqrt(delta**2. + e_k**2.)
    ret_minus = -np.sqrt(delta**2. + e_k**2.)
    return ret_plus, ret_minus

@jit
def SC_dispersion_plus(delta, e_k):
    #only returns the positive root of the superconducting dispersion relation
    ret = np.sqrt(np.power(delta,2.) + np.power(e_k,2.))
    return ret


# In[5]:


@jit
def SC_dispersion_comp_OAAT(delta,e_k):
    #This computes the eigenvalues and eigenvectors of the Supecconductor Hamiltonian matrix H
    #The eigenvalues determines the SC dispersion numerically only one k point at a time
    H = np.array([[e_k,-delta],[-delta,-e_k]])
    evals,evecs = la.eig(H)
    return evals,evecs
@jit
def SC_dispersion_comp(delta,e_k):
    #This determines the SC dispersion numerically at each k point on the grid.
    #It only returns the positive branch
    #by making SC_dispersion_comp_OAAT step through each k point on the grid.
    #the eigen values are sorted using np.max() as positive first and negative second.
    #This is done because Pythpn flips the order of the eigenvalues to (-,+) when the sign of xi flips from + to -
    #therefore to make all the positive/(negative) values to go together
    #the eigen values are sorted using np.max()
    num_k = len(e_k[0])
    ret = np.zeros((num_k,num_k))
    for i in range(0,num_k):
        for j in range(0,num_k):
            evals, evecs = SC_dispersion_comp_OAAT(delta,e_k[i][j])
            Evals = evals.astype(float)
            ret[i][j] = np.max(Evals)
    return ret


# In[6]:


@jit
def SC_coherence_comp_OAAT(delta,k_x,k_y,k_xp,k_yp,mu,a,t1,t2,t3):
    #This determines the coherence matrix from the eigenvectors of H only one k point at a time
    #The eigen vectors need to be sorted based on the sign of E
    #The Bogoliubov matrix (evecs_k) is contructed with the eigenvectors belonging to +E as the first coloumn
    #and those belonging to _E as the second coloumn.
    #(evecs_k) is the matrix of eigenvectors that python autamtically returns.
    #Because of the above mentioned flipping the  vectors need to be sorted too.
    #So if python flips the order of the eigen values to (-,+) we flip back the order of the eigenvectors
    #by putting the first column for the +E and the second coloumn for -E
    #the coherence matrix C is constructed by forming the product transpose(evecs_kp).S.evecs_k


    evals_k, evecs_k = SC_dispersion_comp_OAAT(delta,dispersion(mu,k_x,k_y,a,t1,t2,t3))
    if evals_k[0] < 0.:
        v1 = evecs_k[:,1]
        v2 = evecs_k[:,0]
        evecs_k[:,0] = v1
        evecs_k[:,1] = v2
    evals_kp, evecs_kp = SC_dispersion_comp_OAAT(delta,dispersion(mu,k_xp,k_yp,a,t1,t2,t3))
    if evals_kp[0] < 0.:
        v1 = evecs_kp[:,1]
        v2 = evecs_kp[:,0]
        evecs_kp[:,0] = v1
        evecs_kp[:,1] = v2
    S = np.array([[1.,0.],[0.,-1.]])
    return (np.dot(np.transpose(evecs_kp),np.dot(S,evecs_k))[0][0])**2.
@jit
def SC_coherence_comp(delta,k_x,k_y,k_xp,k_yp,mu,a,t1,t2,t3):
    #This determines the coherence matrix numerically at each k point on the grid.
    #by making SC_coherence_comp_OAAT step through each k point on the grid.
    num_k = len(k_xp[0])
    ret = np.zeros((num_k,num_k))
    for i in range(0,num_k):
        for j in range(0,num_k):
            ret[i][j] = SC_coherence_comp_OAAT(delta,k_x,k_y,k_xp[i][j],k_yp[i][j],mu,a,t1,t2,t3)
    return ret

def coherence(delta,E,k_x,k_y):
    #gives the analytic expression for the coherence factors
    ret = 1. + (np.sqrt(E**2. - delta**2.)*dispersion(0.,k_x,k_y,1.,1.,0.,0.) - delta**2.)/(E*SC_dispersion_plus(delta,dispersion(0.,k_x,k_y,1.,1.,0.,0.)))
    return 0.5*ret

#print coherence(.3,0.4,1.,1.)
#print SC_coherence_comp_OAAT(.3,0.4,1.,1.)


# In[7]:


@jit
def integrator_delta_herm(integrate,dispersion,dk,E,n,f):
    #integrate : 2D array of the function we wish to integrate over k-space
    #dispersion : 2D array with all dispersion relation values over k-space i.e E(k')
    #dk : distance between k values in our grid
    #E : energy that we're currently integrating at ie E(k)
    #n : number of terms in hermite polynomial expansion
    #To summarize we fix E(k) and integrate over all E(k')

    #creating the coefficients of even hermiite polynomials
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
    #exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(2.*f**2.))))
    exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(4.*f))))
    hermite = np.array(np.polynomial.hermite.hermval(E-dispersion_array,coeff))

    ##Multiplying the gaussian with the hermite polynomial expansion##
    #func = (1./np.sqrt(2*np.pi*f))*exp*hermite
    func = exp*hermite*(1./(2.*np.sqrt(np.pi*f)))

    ret = func*integrate_array
    #print integrate_array, func, ret
    #print np.sum(ret)
    return np.sum(ret)*dk*dk
@jit
def integrator_delta_grad(integrate,dispersion,dk,E):
    dispersion_array = np.array(dispersion)
    integrate_array = np.array(integrate)

    delta_mask = dispersion_array == E
    grad_disp_x = np.gradient(dispersion_array,dk)[1]
    grad_disp_y = np.gradient(dispersion_array,dk)[0]

    grad_disp = np.sqrt(grad_disp_x**2. + grad_disp_y**2.)
    grad_disp_delta = grad_disp[delta_mask]

    if (type(integrate) is float) or (type(integrate) is int):
        integrate_delta = integrate_array
    else:
        integrate_delta = integrate_array[delta_mask]
    #print grad_disp_delta
    #print integrate_delta
    if not integrate_delta.all():
        return 0.
    else:
        ret = integrate_delta/grad_disp_delta
        return np.sum(ret)*dk*dk


@jit
def BZ_integrator(function,dk):
    #takes an array over the BZ (function), and integrates it by summing the points up

    func = np.array(function)
    return np.sum(func)*dk*dk


# In[8]:


# def gaussian(e,e_k,f):
#     #e : energy
#     #e_k : dispersion relation
#     #f : full width half max of gaussian

#     #gaussian used to approximate delta function

#     ret = (2./f)*math.sqrt(math.log(2)/math.pi)*math.exp(-(4.*math.log(2.)*(e-e_k)**2.)/f**2.)
#     return ret

@jit
def scattering_SC(delta, a, k_x, k_y, mu, t1, t2, t3, num_k, n, f):
    #delta : superconducting order parameter
    #a : lattice parameter
    #num_k : number of k opints along a single axis

    #forming the kprime grid with meshgrid

    k_xp = np.linspace(-np.pi/a,np.pi/a,num_k)
    k_yp = np.linspace(-np.pi/a,np.pi/a,num_k)
    X,Y = np.meshgrid(k_xp,k_yp)


    dk = 2.*np.pi/(a*num_k)

    # the energy at which we are integrating E(k)
    evals,evecs = SC_dispersion_comp_OAAT(delta,dispersion(mu,k_x,k_y,a,t1,t2,t3))
    E = np.max(evals)

    #computing the inverse scattering rate by integrating over the kprime grid
    Tau = integrator_delta_herm(SC_coherence_comp(delta,k_x,k_y,X,Y,mu,a,t1,t2,t3),SC_dispersion_comp(delta,dispersion(mu,X,Y,a,t1,t2,t3)),dk,E,n,f)
    #Tau = integrator_delta_grad(SC_coherence_comp(delta,E,X,Y,mu,a,t1,t2,t3),SC_dispersion_comp(delta,dispersion(mu,X,Y,a,t1,t2,t3)),dk,E)

    return np.power(Tau,-1)
@jit
def scattering_n(a,E,mu,t1,t2,t3,num_k,n,f):
    #scattering rate in the normal (non-superconducting) state

    k_x = np.linspace(-np.pi/a,np.pi/a,num_k)
    k_y = np.linspace(-np.pi/a,np.pi/a,num_k)
    X,Y = np.meshgrid(k_x,k_y)
    dk = 2.*np.pi/(a*num_k)

    Tau = integrator_delta_herm(1.,dispersion(mu,X,Y,a,t1,t2,t3),dk,E,n,f)
    #Tau = integrator_delta_grad(1.,dispersion(mu,X,Y,a,t1,t2,t3),dk,E)


    return np.power(Tau,-1)



# In[9]:


@jit
def fermi_deriv(E,T):
    #derivative of the fermi function
    energy = np.array(E)
    ret = np.exp(energy/(k_b*T))/((1.+np.exp(energy/(k_b*T)))**2.)
    return ret


# In[10]:


@jit
def FullG(gap,T):
    #selfconsistent equation that determines the gap Delta
    m = 500
    gaps = gap**2.
    f = np.log(T)
    for n in range(1,m+1):
        nh = n - 0.5
        sq2 = np.sqrt((T*nh)**2. + gaps)
        f = 1./nh - T/sq2 + f
    return f
@jit
def Delta(T):
    #solving the self consistent equation using fsolve
    if T < 1.0:
        return np.absolute(2.*np.pi*fsolve(FullG,[0.9],args = T))[0]
    else:
        return 0.0


# In[11]:

@jit
def thermal_conductivity(a,num_k,n,f,num_T,T_end,mu,t1,t2,t3):

    Temp = []
    k_xx = []
    k_xy = []

    dT = T_end/num_T
    T = dT
    #forming the k grid with meshgrid
    k_x = np.linspace(-np.pi/a,np.pi/a,num_k)
    k_y = np.linspace(-np.pi/a,np.pi/a,num_k)
    X,Y = np.meshgrid(k_x,k_y)
    dk = 2.*np.pi/(a*num_k)

    #dispersion values on the k grid of the normal state
    disp = dispersion(mu,X,Y,a,t1,t2,t3)

    disp_SC_comp = np.zeros((int(num_k),int(num_k)))

    #creating the tau values on the k grid for the normal state
    Tau = np.zeros((len(disp),len(disp[0])))
    for i in range(0,len(disp)):
        for j in range(0,len(disp[0])):
            Tau[i][j] = scattering_n(a,disp[i][j],mu,t1,t2,t3,num_k,n,f)

    #normal state E**2 on the grid
    E_squared = np.array(disp)**2.
    #print Tau

    #Velocities in the normal state on the k grid
    #np.gradient returns 2 arrays (the x-deriv & y_deriv), so I just matched them up here
    v_x = np.gradient(disp,dk)[1]
    v_y = np.gradient(disp,dk)[0]



    #Everything for the Sc should be in the temperature loop
    for t in range(0,num_T):

        #The gap values for all T
        delta = Delta(T)

        #the values of the superconducting dispersion on the grid
        disp_SC = SC_dispersion_comp(delta,disp)
        #disp_SC = SC_dispersion_plus(delta,disp)

        """
        for i in range(0,int(num_k)):
            for j in range(0,int(num_k)):
                disp_SC_comp[i][j], junkData = np.max(SC_dispersion_comp_OAAT(delta,disp[i][j]))
        """
        #print disp_SC_comp
        #print " "

        #creating the tau values on the grid for the superconducting state
        Tau_SC = np.zeros((len(disp_SC),len(disp_SC[0])))
        for i in range(0,len(disp_SC)):
            for j in range(0,len(disp_SC[0])):
                Tau_SC[i][j] = scattering_SC(delta,a,X[i][j],Y[i][j],mu,t1,t2,t3,num_k,n,f)
        #print Tau_SC

        #Superconducting E**2 on the grid
        E_squared_SC = np.power(disp_SC,2.)


        #Velocities in superconducting state on the k grid
        v_x_SC = np.gradient(disp_SC,dk)[1]
        v_y_SC = np.gradient(disp_SC,dk)[0]

        #the fermi function derivatives of the normal and superconducting state on the grid
        fermi = fermi_deriv(disp,T)
        fermi_SC = fermi_deriv(disp_SC,T)

        #storing the values of the temperature
        Temp.append(T)

        #normal state thermal conductivities
        k_xx_n = BZ_integrator(E_squared*Tau*fermi*v_x*v_x,dk)/(T*T)
        k_xy_n = BZ_integrator(E_squared*Tau*fermi*v_x*v_y,dk)/(T*T)

        #superconducting state thermal conductivities
        k_xx_SC = BZ_integrator(E_squared_SC*Tau_SC*fermi_SC*v_x_SC*v_x_SC,dk)/(T*T)
        k_xy_SC = BZ_integrator(E_squared_SC*Tau_SC*fermi_SC*v_x_SC*v_y_SC,dk)/(T*T)

        #storing normalized thermal conductivities
        k_xx.append(k_xx_SC/k_xx_n)
        k_xy.append(k_xy_SC/k_xx_n)

        print("{0:<20} {1:<20} {2:<20} {3:<20} {4:<}".format(T, delta, k_xx_SC/k_xx_n, k_xy_SC/k_xx_n,k_xx_SC, k_xx_n))
        #print (T, delta, k_xx_SC/k_xx_n, k_xy_SC/k_xx_n,k_xx_SC, k_xx_n)
        T = T + dT
        #Temperature = T

#    get_ipython().run_line_magic('matplotlib', 'qt')
    # plt.plot(Temp,k_xx)
    # plt.plot(Temp,k_xy)
    # plt.show()

#thermal_conductivity(a,num_k,n,f,num_T,T_end,mu,t1,t2,t3)
thermal_conductivity(1.,10.,3,.01,30,1.1,0.29,-1.,0.4,0.)
#scattering_SC(0.3,1.,2.,200,3,.1)
