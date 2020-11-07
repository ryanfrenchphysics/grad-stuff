#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
import scipy.integrate as inte
from mpl_toolkits.mplot3d import Axes3D
import ray
#%matplotlib


# In[2]:

@jit
def dispersion(mu,k_x,k_y,a,t1,t2,t3):
    #mu : chemical potential
    #k_x : x - momentum (p = hbar*k)
    #k_y : y - momentum
    #a : lattice parameter
    #t1 : nearest neighbor hopping integral
    #t2 : next nearest neighbor hopping
    #t3 : next next nearest neighbor hopping
    ret=  mu-t1*(np.cos(k_x*a)+np.cos(k_y*a))-t2*(np.cos(k_x*a)*np.cos(k_y*a))
    #ret = mu + 0.5*t1*(np.cos(k_x*a)+np.cos(k_y*a)) + t2*np.cos(k_x*a)*np.cos(k_y*a) + 0.5*t3*(np.cos(2.*k_x*a)+np.cos(2.*k_y*a))
    return ret


# In[3]:


# def SDW_dispersion_comp_OAAT(M,Q,mu,k_x,k_y,a,t1,t2,t3):
#     if (k_x >= 0. and k_y >= 0.) or (k_x <= 0. and k_y <= 0.):
#         m = M
#     else:
#         m = 0.
#     E = dispersion(mu,k_x,k_y,a,t1,t2,t3)
#     E_q = dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3)
#     H = np.array([[E,m],[m,E_q]])
#     evals, evecs = la.eig(H)
#     return np.sort(evals)


# In[4]:

@ray.remote()
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
    return evals


# In[5]:


# def SDW_dispersion_RBZ(M,Q,mu,k_x,k_y,a,t1,t2,t3):
#     #M : Spin Density Wave order parameter
#     #Q : nesting vector

#     xi_plus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) + dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))

#     xi_minus = 0.5*(dispersion(mu,k_x,k_y,a,t1,t2,t3) - dispersion(mu,k_x+Q[0],k_y+Q[1],a,t1,t2,t3))
#     ret_plus = xi_plus + np.sqrt(M**2. + xi_minus**2.)
#     ret_minus = xi_plus - np.sqrt(M**2. + xi_minus**2.)
#     return ret_plus, ret_minus


# In[6]:

@ray.remote()
@jit
def fermi_surface(mu,a,t1,t2,t3,num_k):
    #see above for what these parameters mean
    #num_k : number of points in 1d in k-space

    k_x = np.linspace(-np.pi/a,np.pi/a,num_k)
    k_y = np.linspace(-np.pi/a,np.pi/a,num_k)
    X,Y = np.meshgrid(k_x,k_y)
    function = dispersion(mu,X,Y,a,t1,t2,t3) - mu
    X_new = X/np.pi
    Y_new = Y/np.pi
#     %matplotlib inline
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cpf = ax.contourf(X_new,Y_new,function,10)
    #colours = ['w' if level<0 else 'k' for level in cpf.levels]
    # cp = ax.contour(X, Y, function, 20)
    # ax.clabel(cp, fontsize=12)
    plt.show()
#fermi_surface(mu,a,t1,t2,t3,num_k)
#fermi_surface(0.29,1.0,-1.0,0.4,-0.1,100)



@jit
def band_structure(mu,a,t1,t2,t3,M,Q,num_steps):
    #most parameters are already defined in above functions
    #num_steps : number of steps along each path

    """
        |-------------------|
        |        /        |
        |       /        |
        |      /        |
        |      ----------|
        |            |
        |            |
        |            |
        |-------------------|

    """
    #Path 1 is the diagonal
    #Path 2 is down
    #Path 3 is to the left

    #steps_1 = int(np.sqrt(2)*num_steps)
    n = 0

    N = []
    BS = []
    SDW_plus = []
    SDW_minus = []
    SDW_plus_comp = []
    SDW_minus_comp = []


    X = []
    Y = []

    for i in range(0,num_steps):
    #this is the for loop for path 1
        k_x = (((np.pi/a)/float(num_steps-1))*float(i))
        k_y = (k_x) #we're going along a 45 degree angle

        X.append(k_x)
        Y.append(k_y)

#         plus, minus = SDW_dispersion_RBZ(M,Q,mu,k_x,k_y,a,t1,t2,t3)
#         SDW_plus.append(plus)
#         SDW_minus.append(minus)

        minus_comp, plus_comp = SDW_dispersion_comp_OAAT.remote(M,Q,mu,k_x,k_y,a,t1,t2,t3)
        SDW_plus_comp.append(plus_comp)
        SDW_minus_comp.append(minus_comp)



        N.append(n)
        BS.append(dispersion(mu,k_x,k_y,a,t1,t2,t3))

        n = n + 1
    pi_pi = n

    for i in range(0,num_steps):
    #this is the loop for path 2
        k_x = (np.pi)
        k_y = (np.pi/a - float(i)*((np.pi/a)/float(num_steps-1)))

        X.append(k_x)
        Y.append(k_y)

#         plus, minus = SDW_dispersion_RBZ(M,Q,mu,k_x,k_y,a,t1,t2,t3)
#         SDW_plus.append(plus)
#         SDW_minus.append(minus)

        minus_comp, plus_comp = SDW_dispersion_comp_OAAT.remote(M,Q,mu,k_x,k_y,a,t1,t2,t3)
        SDW_plus_comp.append(plus_comp)
        SDW_minus_comp.append(minus_comp)


        N.append(n)
        BS.append(dispersion(mu,k_x,k_y,a,t1,t2,t3))

        n = n + 1
    pi_0 = n

    for i in range(0,num_steps):
    #this is the loop for path 3
        k_x = (np.pi/a - float(i)*((np.pi/a)/float(num_steps-1)))
        k_y = 0.

        X.append(k_x)
        Y.append(k_y)

#         plus, minus = SDW_dispersion_RBZ(M,Q,mu,k_x,k_y,a,t1,t2,t3)
#         SDW_plus.append(plus)
#         SDW_minus.append(minus)

        minus_comp, plus_comp = SDW_dispersion_comp_OAAT.remote(M,Q,mu,k_x,k_y,a,t1,t2,t3)
        SDW_plus_comp.append(plus_comp)
        SDW_minus_comp.append(minus_comp)


        N.append(n)
        BS.append(dispersion(mu,k_x,k_y,a,t1,t2,t3))

        n = n + 1
    plt.figure(2)
#     %matplotlib qt
    plt.plot(N,BS,label='Tight Binding')
#     plt.plot(N,SDW_plus,label='Tight Binding with SDW (plus band)')
#     plt.plot(N,SDW_minus,label = 'Tight Binding with SDW (minus band)')
    plt.plot(N,SDW_plus_comp,label = 'Tight Binding with SDW (plus_comp)')
    plt.plot(N,SDW_minus_comp,label = 'Tight Binding with SDW (minus_comp)')

    plt.axhline(y=0)
    x_axis = [0,pi_pi,pi_0,n] #the only tic marks we want on the x-axis are when (k_x,k_y) = [(0,0) , (pi,pi), (pi, 0), (0,0)]
    x_labels = ['(0,0)','($\pi$/a,$\pi$/a)','($\pi$/a,0)','(0,0)'] # these are the labels at the proper points
    plt.xticks(x_axis,x_labels)
    plt.axis([0,n,-2.,2.])
    #plt.legend()
    #plt.scatter(X,Y)
    plt.show()



# In[8]:


#fermi_surface(mu,a,t1,t2,t3,num_k)
#fermi_surface(0.29,1.0,-1.0,0.4,-0.1,100)


# In[9]:

if __name__ == '__main__':

    #band_structure(mu,a,t1,t2,t3,M,Q,num_steps)
    band_structure(0.29,1.0,-1.0,0.4,0.0,0.1,[np.pi,np.pi],100)


# In[ ]:
