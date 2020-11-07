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
from numba import jit, prange
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import multiprocessing as mp

#%matplotlib qt
#hbar = 1.05e-34
#m_e = 9.11e-31
#k_b = 1.38e-23

hbar = 1.
m_e = 1.
k_b = 1.
NUMK = 10
fps = 10 # frame per sec
frn = 50 # frame number of the animation

#Temperature = 0.1


# In[3]:
@jit
def dispersion(mu,k_x,k_y,a,t1,t2,t3):
    #mu : chemical potential
    #k_x : x - momentum (p = hbar*k)
    #k_y : y - momentum
    #a : lattice parameter
    #t1 : nearest neighbor hopping integral
    #t2 : next nearest neighbor hopping
    #t3 : next next nearest neighbor hopping
    #ret = mu + 0.5*t1*(np.cos(k_x*a)+np.cos(k_y*a)) + t2*np.cos(k_x*a)*np.cos(k_y*a) + 0.5*t3*(np.cos(2.*k_x*a)+np.cos(2.*k_y*a))

    ret = mu - t1*(np.cos(k_x*a)+np.cos(k_y*a)) - t2*np.cos(k_x*a)*np.cos(k_y*a)
    return ret
#print(dispersion(0,1.0471975511965974, -3.141592653589793,1,-1,0,0))


# In[4]:

@jit
def SDW_dispersion_comp_OAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3):
    if (k_x >= 0. and k_y >= 0.) or (k_x <= 0. and k_y <= 0.):
        m = M
    else:
        m = 0.
    E = dispersion(mu,k_x,k_y,a,t1,t2,t3)
    E_q = dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3)
    H = np.array([[E,m],[m,E_q]])
    evals, evecs = la.eig(H)
    xi_minus = 0.5*(E - E_q)
    if xi_minus < 0.:
        e_alpha = evals[1]
        e_beta = evals[0]
        evals[0] = e_alpha
        evals[1] = e_beta
        evec_alpha = evecs[:,1]
        evec_beta = np.array([evecs[0][0],evecs[1][0]]) #I wrote evec_beta this way because it was just setting the pointers equal to each other then overwriting evec_beta later since we changed evecs[:,0] to evec_alpha
        evecs[:,0] = evec_alpha
        evecs[:,1] = evec_beta
    return evals,evecs

@jit(parallel=True)
def SDW_dispersion_comp(M,Q,mu,k_x,k_y,a,t1,t2,t3,num_k):
    disp_a = np.zeros((int(num_k),int(num_k)))
    disp_b = np.zeros((int(num_k),int(num_k)))
    disp_ooat = SDW_dispersion_comp_OAAT
    for i in prange(0,int(num_k)):
        for j in prange(0,int(num_k)):
            evals,evecs = disp_ooat(M,Q,mu,k_x[i][j],k_y[i][j],a,t1,t2,t3)
            disp_a[i][j] = evals[0]
            disp_b[i][j] = evals[1]
    return disp_a, disp_b


@jit
def SDW_xi(Q,mu,k_x,k_y,a,t1,t2,t3):
    #function that just returns xi_k_plus/minus mostly for use in SDW_coherence()
    xi_plus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) + dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))

    xi_minus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) - dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))
    return xi_plus, xi_minus


# In[5]:


# def SDW_dispersion(M,Q,mu,k_x,k_y,a,t1,t2,t3):
#     #M : Spin Density Wave order parameter
#     #Q : nesting vector

#     X = np.array(k_x)
#     Y = np.array(k_y)

#     ix = int(np.sqrt(X.size))
#     iy = int(np.sqrt(Y.size))

#     ret = np.zeros((ix,iy))

#     for i in range(0,ix):
#         for j in range(0,iy):

#             if (X[i][j] > 0. and Y[i][j] > 0.) or (X[i][j] < 0. and Y[i][j] < 0):
#                 xi_plus = 0.5*(dispersion(mu,X[i][j],Y[i][j],a,t1,t2,t3) + dispersion(mu,X[i][j]+Q[0],Y[i][j]+Q[1],a,t1,t2,t3))

#                 xi_minus = 0.5*(dispersion(mu,X[i][j],Y[i][j],a,t1,t2,t3) - dispersion(mu,X[i][j]+Q[0],Y[i][j]+Q[1],a,t1,t2,t3))
#                 ret[i][j] = xi_plus + np.sign(xi_minus)*np.sqrt(M**2. + xi_minus**2.)
#     #ret_minus = xi_plus - np.sqrt(M**2. + xi_minus**2.)
#             else:
#                 ret[i][j] = dispersion(mu,X[i][j],Y[i][j],a,t1,t2,t3)
#             #print X[i][j], Y[i][j], ret[i][j]
#     return ret


# In[6]:

@jit
def SDW_coherence_comp_OAAT(M,Q,k_x,k_y,k_xp,k_yp,mu,a,t1,t2,t3):
    #numerically computes the coherence factors one k point at a time
    evals,B = SDW_dispersion_comp_OAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3)
    evals,Bprime = SDW_dispersion_comp_OAAT(M,Q,mu,k_xp,k_yp,a,t1,t2,t3)
    S = np.array([[1.,1.],[1.,1.]])
    U = np.dot(np.transpose(Bprime),np.dot(S,B))
    Caa = U[0][0]**2.
    Cba = U[0][1]**2
    Cab = U[1][0]**2
    Cbb = U[1][1]**2.
    #print (Cbb)
    return Caa, Cba, Cab, Cbb

@jit(parallel=True)
def SDW_coherence_comp(M,Q,k_x,k_y,k_xp,k_yp,mu,a,t1,t2,t3,num_k):
    #takes the numerical coherence factors and steps through all points on the k prime grid
    Caa = np.zeros((num_k,num_k))
    Cba = np.zeros((num_k,num_k))
    Cab = np.zeros((num_k,num_k))
    Cbb = np.zeros((num_k,num_k))
    coh_oaat = SDW_coherence_comp_OAAT
    for i in prange(0,num_k):
        for j in prange(0,num_k):
            Caa[i][j],Cba[i][j],Cab[i][j],Cbb[i][j] = coh_oaat(M,Q,k_x,k_y,k_xp[i][j],k_yp[i][j],mu,a,t1,t2,t3)

    # Caait = ray.put(Caa)
    # Cabit = ray.put(Cab)
    # Cbait = ray.put(Cba)
    # Cbbit = ray.put(Cbb)

    return Caa,Cba,Cab,Cbb


# In[7]:

@jit
def SDW_coherence_OAAT(M,Q,k_x,k_y,k_xp,k_yp,mu,a,t1,t2,t3):
    #Analytically computes the coherence factors one k point at a time at a time

    if (k_x >= 0. and k_y >= 0.) or (k_x <= 0. and k_y <= 0.):
        m = M
    else:
        m = 0.

    if (k_xp >= 0. and k_yp >= 0.) or (k_xp <= 0. and k_yp <= 0.):
        m_prime = M
    else:
        m_prime = 0.

    E = dispersion(mu,k_x,k_y,a,t1,t2,t3)
    E_q = dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3)
    xi_minus = 0.5*(E - E_q)

    Ep = dispersion(mu,k_xp,k_yp,a,t1,t2,t3)
    E_qp = dispersion(mu,k_xp+Q[0],k_yp+Q[1],a,t1,t2,t3)
    xi_minusp = 0.5*(Ep - E_qp)

    Lambda = np.sqrt((xi_minus)**2. + m**2.)
    Lambda_prime = np.sqrt((xi_minusp)**2. + m_prime**2.)

    Caa = (1. + (m/Lambda) + (m_prime/Lambda_prime) + (m*m_prime)/(Lambda*Lambda_prime))
    Cba = (1. - (m/Lambda) + (m_prime/Lambda_prime) - (m*m_prime)/(Lambda*Lambda_prime))
    Cab = (1. + (m/Lambda) - (m_prime/Lambda_prime) - (m*m_prime)/(Lambda*Lambda_prime))
    Cbb = (1. - (m/Lambda) - (m_prime/Lambda_prime) + (m*m_prime)/(Lambda*Lambda_prime))
    #print(Caa)
    return Caa, Cba, Cab, Cbb

@jit(parallel=True)
def SDW_coherence(M,Q,k_x,k_y,k_xp,k_yp,mu,a,t1,t2,t3,num_k):
    #takes the numerical coherence factors and steps through all points on the k prime grid

    Caa = np.zeros((num_k,num_k))
    Cba = np.zeros((num_k,num_k))
    Cab = np.zeros((num_k,num_k))
    Cbb = np.zeros((num_k,num_k))
    sdw_oaat = SDW_coherence_OAAT
    for i in prange(0,num_k):
        for j in prange(0,num_k):
            Caa[i][j],Cba[i][j],Cab[i][j],Cbb[i][j] = sdw_oaat(M,Q,k_x,k_y,k_xp[i][j],k_yp[i][j],mu,a,t1,t2,t3)

    return Caa,Cba,Cab,Cbb


# In[8]:


# def gaussian(e,e_k,f):
#     #e : energy
#     #e_k : dispersion relation
#     #f : full width half max of gaussian

#     #gaussian used to approximate delta function

#     ret = (2./f)*np.sqrt(np.log(2)/np.pi)*np.exp(-(4.*np.log(2.)*(e-e_k)**2.)/f**2.)
#     return ret

@jit(parallel=True)
def integrator_delta_herm(integrate,dispersion,dk,E,n,f):
    #integrate : 2D array of the function we wish to integrate over k-space
    #dispersion : 2D array with all dispersion relation values over k-space
    #dk : distance between k values in our grid
    #E : energy that we're currently integrating at
    #n : number of terms in hermite polynomial expansion

    coeff = np.zeros(n+1)
    for i in prange(0,len(coeff)):
        if i % 2 == 0:
            coeff[i] = ((-1)**float(i))/(math.factorial(float(i))*4**float(i)*np.sqrt(np.pi))
        else:
            coeff[i] = 0.

    dispersion_array = np.array(dispersion)
    integrate_array = np.array(integrate)

    #exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(2.*f**2.))))
    exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(4.*f))))

    hermite = np.array(np.polynomial.hermite.hermval(E-dispersion_array,coeff))
    #func = (1./np.sqrt(2*np.pi*f))*exp*hermite
    func = (1./(2.*np.sqrt(np.pi*f)))*exp*hermite
    ret = func*integrate_array

    return np.sum(ret)*dk*dk

@jit
def BZ_integrator(function,dk):
    #takes an array over the BZ (function), and integrates it by summing the points up

    func = np.array(function)
    return np.sum(func)*dk*dk


# In[9]:


@jit
def scattering_SDW_OAAT(M,Q,k_x,k_y,mu,a,t1,t2,t3,num_k,n,f):
    #k' grid to be integrated over for scattering
    k_xp = np.linspace(-(np.pi)/a,(np.pi)/a,num_k)
    k_yp = np.linspace(-(np.pi)/a,(np.pi)/a,num_k)
    Xp,Yp = np.meshgrid(k_xp,k_yp)
    dk = 2.*np.pi/(a*num_k)

    Caa,Cba,Cab,Cbb = SDW_coherence_comp(M,Q,k_x,k_y,Xp,Yp,mu,a,t1,t2,t3,num_k)
    #Caa,Cba,Cab,Cbb = SDW_coherence.remote(M,Q,k_x,k_y,Xp,Yp,mu,a,t1,t2,t3,num_k)

    disp_prime_a, disp_prime_b = SDW_dispersion_comp(M,Q,mu,Xp,Yp,a,t1,t2,t3,num_k)
    E,evecs = SDW_dispersion_comp_OAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3)

    E_a = E[0]
    E_b = E[1]
    #these are the inverse scattering rates
    Tau_aa = integrator_delta_herm(Caa,disp_prime_a,dk,E_a,n,f)
    Tau_ba = integrator_delta_herm(Cba,disp_prime_a,dk,E_b,n,f)
    Tau_ab = integrator_delta_herm(Cab,disp_prime_b,dk,E_a,n,f)
    Tau_bb = integrator_delta_herm(Cbb,disp_prime_b,dk,E_b,n,f)



    #print(Tau_bb)

    Tau_a = Tau_aa  #+ Tau_ab #these are the effective inverse scattering rates
    Tau_b = Tau_bb  #+ Tau_ba


    threshold =0 #set this as the minimum value for the inverse scattering rates

    if Tau_a < threshold:
        Tau_a = 0.
    else:
        Tau_a = 1./Tau_a

    if Tau_b < threshold:
        Tau_b = 0.
    else:
        Tau_b = 1./Tau_b

    return Tau_a, Tau_b

@jit(parallel=True)
def scattering_SDW(M,Q,k_x,k_y,mu,a,t1,t2,t3,num_k,n,f):
    """
    #if you want to ONLY call this function replace k_x[i][j] & k_y[i][j] with X[i][j] & Y[i][j] below
    x = np.linspace(-np.pi/a,np.pi/a)
    y = np.linspace(-np.pi/a,np.pi/a)
    X,Y = np.meshgrid(x,y)
    """
    Tau_a = np.zeros([num_k,num_k])
    Tau_b = np.zeros([num_k,num_k])
    scat = scattering_SDW_OAAT
    for i in prange(0,num_k):
        for j in prange(0,num_k):
            Tau_a[i][j],Tau_b[i][j] = scat(M,Q,k_x[i][j],k_y[i][j],mu,a,t1,t2,t3,num_k,n,f)

    return Tau_a, Tau_b

@jit
def scattering_n(a,E,mu,t1,t2,t3,num_k,n,f):
    #scattering rate in the normal (non-superconducting) state

    k_x = np.linspace(-(np.pi)/a,(np.pi)/a,num_k)
    k_y = np.linspace(-(np.pi)/a,(np.pi)/a,num_k)
    X,Y = np.meshgrid(k_x,k_y)
    dk = 2.*np.pi/(a*num_k)

    #Tau = integrator_delta_herm(1.,quad_dispersion(X,Y),dk,E,n,f)
    Tau = integrator_delta_herm(1.,dispersion(mu,X,Y,a,t1,t2,t3),dk,E,n,f)

    return np.power(Tau,-1)


# In[10]:



@jit
def fermi_deriv(E,T):
    #derivative of the fermi function
    energy = np.array(E)
    ret = np.exp(energy/(k_b*T))/((1.+np.exp(energy/(k_b*T)))**2.)
    return ret


# In[11]:


@jit(parallel=True)
def FullG(gap,T):
    #snagged from Sourav's code
    m = 500
    #global Temperature
    gaps = gap**2.
    f = np.log(T)
    for n in prange(1,m+1):
        nh = n - 0.5
        sq2 = np.sqrt((T*nh)**2. + gaps)
        f = 1./nh - T/sq2 + f
    return f

@jit
def Delta(T):
    #global Temperature
    if T < 1.0:
        return np.absolute(2.*np.pi*fsolve(FullG,[0.9],args = T))[0]
    else:
        return 0.0

#def SDW_dispersion_comp(M,Q,mu,k_x,k_y,a,t1,t2,t3,num_k):

'''
TODO: GET THIS TO WORK
'''
def update_plot(frame_num, zarray, plot):
    plot[0].remove()
    plot[0] = ax.plot_surface(X, Y, zarray[:,:,frame_num], cmap=plt.cm.coolwarm, linewidth=0)

def run_animation(X, Y, z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    plot = [ax.plot_surface(X, Y, z[:,:,0], rstride=1, cstride=1)]

    ani = animation.FuncAnimation(fig, update_plot, frn, fargs=(z, plot), interval=1000)
    ani.save('mp4'+'.mp4',writer='ffmpeg',fps=fps)

##@jit
def normv(vxa, vya, vxb, vyb):
    va = np.sqrt(vxa**2 + vya**2)
    vb = np.sqrt(vxb**2 + vyb**2)
    return va, vb


@jit(parallel=True)
def thermal_conductivity(a,num_k,n,f,num_T,T_end,Q,mu,t1,t2,t3):


    #global Temperature

    Temp = []
    k_xx = []
    k_xy = []

    dT = T_end/num_T
    T = .1
    #Temperature = T

    k_x = np.linspace(-(np.pi)/a,(np.pi)/a,num_k)
    k_y = np.linspace(-(np.pi)/a,(np.pi)/a,num_k)
    X,Y = np.meshgrid(k_x,k_y)
    dk = 2.*np.pi/(a*num_k)

    #disp = quad_dispersion(X,Y)

    disp_SDW_a = np.zeros((int(num_k),int(num_k)))
    disp_SDW_b = np.zeros((int(num_k),int(num_k)))


    disp = dispersion(mu,X,Y,a,t1,t2,t3)
    #Tau = np.zeros((len(disp),len(disp[0])))
    # for i in range(0,len(disp)):
    #     for j in range(0,len(disp[0])):
    #         Tau[i][j] = scattering_n(a,disp[i][j],mu,t1,t2,t3,num_k,n,f)
    #E_squared = np.array(disp)**2.

    v_x = np.gradient(disp,dk)[1] #np.gradient returns 2 arrays (the x-deriv & y_deriv), so I just matched them up here
    v_y = np.gradient(disp,dk)[0]
    '''
    disp_SDW_a,disp_SDW_b = SDW_dispersion_comp(1.76,Q,mu,X,Y,a,t1,t2,t3,num_k)
    v_x_SDW_a = np.gradient(disp_SDW_a,dk)[1]
    v_y_SDW_a = np.gradient(disp_SDW_a,dk)[0]
    v_x_SDW_b = np.gradient(disp_SDW_b,dk)[1]
    v_y_SDW_b = np.gradient(disp_SDW_b,dk)[0]

    vnormSDWa, vnormSDWb = normv(v_x_SDW_a, v_y_SDW_a, v_x_SDW_b, v_y_SDW_b)
    vnormSDWa.reshape((NUMK,NUMK))
    vnormSDWb.reshape((NUMK,NUMK))
    # vnorm, nonsense = normv(v_x, v_y, v_x, v_y)
    # vnorm = vnorm.reshape((NUMK, NUMK))
    # vnorminv = 1/vnorm
    vnormSDWainv = 1/vnormSDWa
    vnormSDWbinv = 1/vnormSDWb

    for i in range(NUMK):
        for j in range(NUMK):
            print("kx: " + str(X[i,j]) + "\tky: " + str(Y[i,j]) + "\t1/normSDW: " + str(vnormSDWbinv[i,j]))
            print("\n")
    '''

    disp_a, disp_b = SDW_dispersion_comp(0.,Q,mu,X,Y,a,t1,t2,t3,num_k)
    Tau_a, Tau_b = scattering_SDW(0.,Q,X,Y,mu,a,t1,t2,t3,num_k,n,f)
    E_squared_a = np.power(disp_a,2.)
    E_squared_b = np.power(disp_b,2.)
    v_x_a = np.gradient(disp_a,dk)[1]
    v_y_a = np.gradient(disp_a,dk)[0]
    v_x_b = np.gradient(disp_b,dk)[1]
    v_y_b = np.gradient(disp_b,dk)[0]

    z = np.zeros((NUMK, NUMK, frn))


    for t in prange(0,num_T):
        M = Delta(T)
        disp_SDW_a,disp_SDW_b = SDW_dispersion_comp(M,Q,mu,X,Y,a,t1,t2,t3,num_k)

        Tau_SDW_a = np.zeros((num_k,num_k))
        Tau_SDW_b = np.zeros((num_k,num_k))
        Tau_SDW_a, Tau_SDW_b = scattering_SDW(M,Q,X,Y,mu,a,t1,t2,t3,num_k,n,f)
        z[:,:,t] = Tau_SDW_b

        E_squared_SDW_a = np.power(disp_SDW_a,2.)
        E_squared_SDW_b = np.power(disp_SDW_b,2.)

        v_x_SDW_a = np.gradient(disp_SDW_a,dk)[1]
        v_y_SDW_a = np.gradient(disp_SDW_a,dk)[0]
        v_x_SDW_b = np.gradient(disp_SDW_b,dk)[1]
        v_y_SDW_b = np.gradient(disp_SDW_b,dk)[0]

        fermi_SDW_a = fermi_deriv(disp_SDW_a,T)
        fermi_SDW_b = fermi_deriv(disp_SDW_b,T)

        """
        fermi = fermi_deriv(disp,T)
        """


        fermi_a = fermi_deriv(disp_a,T)
        fermi_b = fermi_deriv(disp_b,T)

        Temp.append(T)
        """
        k_xx_nwm = BZ_integrator(E_squared*Tau*fermi*v_x*v_x,dk)/(T*T)
        k_xy_nwm = BZ_integrator(E_squared*Tau*fermi*v_x*v_y,dk)/(T*T)
        """
        k_xx_a = BZ_integrator(E_squared_a*Tau_a*fermi_a*v_x_a*v_x_a,dk)/(T*T)
        k_xx_b = BZ_integrator(E_squared_b*Tau_b*fermi_b*v_x_b*v_x_b,dk)/(T*T)

        k_xy_a = BZ_integrator(E_squared_a*Tau_a*fermi_a*v_x_a*v_y_a,dk)/(T*T)
        k_xy_b = BZ_integrator(E_squared_b*Tau_b*fermi_b*v_x_b*v_y_b,dk)/(T*T)

        #k_xx_SDW_a = BZ_integrator(E_squared_SDW_a*(Tau_SDW_a+Tau_SDW_b)*fermi_SDW_a*v_x_SDW_a*v_x_SDW_a,dk)/(T*T)
        #k_xy_SDW_a = BZ_integrator(E_squared_SDW_a*(Tau_SDW_a+Tau_SDW_b)*fermi_SDW_a*v_x_SDW_a*v_y_SDW_a,dk)/(T*T)

        k_xx_SDW_a = BZ_integrator(E_squared_SDW_a*(Tau_SDW_a)*fermi_SDW_a*v_x_SDW_a*v_x_SDW_a,dk)/(T*T)
        k_xy_SDW_a = BZ_integrator(E_squared_SDW_a*(Tau_SDW_a)*fermi_SDW_a*v_x_SDW_a*v_y_SDW_a,dk)/(T*T)


        #k_xx_SDW_b = BZ_integrator(E_squared_SDW_b*(Tau_SDW_b+Tau_SDW_a)*fermi_SDW_b*v_x_SDW_b*v_x_SDW_b,dk)/(T*T)
        #k_xy_SDW_b = BZ_integrator(E_squared_SDW_b*(Tau_SDW_b+Tau_SDW_a)*fermi_SDW_b*v_x_SDW_b*v_y_SDW_b,dk)/(T*T)

        k_xx_SDW_b = BZ_integrator(E_squared_SDW_b*(Tau_SDW_b)*fermi_SDW_b*v_x_SDW_b*v_x_SDW_b,dk)/(T*T)
        k_xy_SDW_b = BZ_integrator(E_squared_SDW_b*(Tau_SDW_b)*fermi_SDW_b*v_x_SDW_b*v_y_SDW_b,dk)/(T*T)



        k_xx_SDW = k_xx_SDW_a + k_xx_SDW_b
        k_xy_SDW = k_xy_SDW_a + k_xy_SDW_b

        k_xx_n = k_xx_a + k_xx_b
        k_xy_n = k_xy_a+k_xy_b

#         k_xx.append(k_xx_SDW_b/k_xx_n)
#         k_xy.append(k_xy_SDW_b/k_xx_n)

        k_xx.append(k_xx_SDW_b)
        k_xy.append(k_xy_SDW_b)


        # print (T, M, k_xx_SDW/k_xx_n, k_xy_SDW/k_xx_n, k_xx_SDW_a, k_xx_SDW_b, k_xx_n)
        # print(Tau_a)
        #print(Tau_b)
        #print(v_x_SDW_a)
        #print(v_y_SDW_a)
        #print(v_x_SDW_b)
        #print(v_y_SDW_b)
        #print(v_x_a)
        #print(v_y_a)
        #print(v_x_b)
        #print(v_y_b)
        #print(fermi_a)
        #print(fermi_b)
        #print("{0:<20} {1:<20} {2:<20} {3:<20}".format(k_xx_n, k_xx_nwm,k_xy_n, k_xy_nwm))

        #print("{0:<20} {1:<20} {2:<20} {3:<20} {4:<}".format(T, M, k_xx_SDW/k_xx_n, k_xy_SDW/k_xx_n, k_xx_n))
        T = T + dT
        #Temperature = T

        va, vb = normv(v_x, v_y, v_x, v_y)



        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax.plot_surface(X, Y, Tau_SDW_b, cmap=plt.cm.coolwarm, linewidth=0)
        ax.plot_surface(X, Y, (Tau_SDW_b), cmap=plt.cm.coolwarm, linewidth=0)

        plt.title("Temp:" + str(T))
        plt.show()


if __name__ == '__main__':
    #thermal_conductivity(a,num_k,n,f,num_T,T_end,Q,mu,t1,t2,t3)
    thermal_conductivity(1.,NUMK,3,.01,1,1.02,[np.pi,np.pi],0.29,-1.,0.4,0.)
