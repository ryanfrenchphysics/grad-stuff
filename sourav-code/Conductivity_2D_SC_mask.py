import math
import numpy as np
from numpy import linalg as la
from matplotlib import pyplot as plt
import scipy.integrate as inte
from scipy.optimize import fsolve
from numba import jit
from numba.extending import overload
import multiprocessing as mp
import time
import warnings
warnings.filterwarnings('ignore')

# m_e = 9.11e-31
# k_b = 1.38e-23


@overload(fsolve)
def jit_fsolve(fxn, x0, args=None):
    @jit(nopython=True)
    def use(f, x, _a):
        if args is None:
            return fsolve(fxn, x0)
        else:
            return fsolve(fxn, x0, args)


hbar = 1.
m_e = 1.
k_b = 1.

# Temperature = 0.1
@jit
def dispersion(mu, k_x, k_y, a, t1, t2, t3):
    # mu : chemical potential
    # k_x : x - momentum (p = hbar*k)
    # k_y : y - momentum
    # a : lattice parameter
    # t1 : nearest neighbor hopping integral
    # t2 : next nearest neighbor hopping
    # t3 : next next nearest neighbor hopping

    ret = mu + t1 * (np.cos(k_x * a) + np.cos(k_y * a)) + t2 * np.cos(
        k_x * a) * np.cos(k_y * a) + 0.5 * t3 * (np.cos(2. * k_x * a) + np.cos(2. * k_y * a))
    return np.absolute(ret)


@jit
def SC_dispersion(delta, e_k):
    # delta : super conducting order parameter
    # e_k = : dispersion relation in normal state

    ret_plus = np.sqrt(delta**2. + e_k**2.)
    ret_minus = -np.sqrt(delta**2. + e_k**2.)
    return ret_plus, ret_minus


def quad_dispersion(k_x, k_y):
    ret = 0.5 * (np.power(k_x, 2.) + np.power(k_y, 2.))
    return ret


@jit
def SC_dispersion_plus(delta, e_k):
    # only returns the positive root of the superconducting dispersion relation
    ret = np.sqrt(np.power(delta, 2.) + np.power(e_k, 2.))
    return ret


@jit
def SC_dispersion_comp_OAAT(delta, e_k):
    M = np.array([[e_k, -delta], [-delta, -e_k]])
    evals, evecs = la.eig(M)
    return evals, evecs


@jit
def SC_dispersion_comp(delta, e_k):
    num_k = len(e_k[0])
    ret = np.zeros((num_k, num_k))
    for i in range(0, num_k):
        for j in range(0, num_k):
            evals, evecs = SC_dispersion_comp_OAAT(delta, e_k[i][j])
            Evals = evals.astype(float)
            ret[i][j] = np.max(Evals)
    return ret


@jit
def SC_dispersion_comp_new(delta, e_k):
    num_k = len(e_k)
    ret = np.zeros(num_k)
    for i in range(0, num_k):
        evals, evecs = SC_dispersion_comp_OAAT(delta, e_k[i])
        Evals = evals.astype(float)
        ret[i] = np.max(Evals)
    return ret


def gaussian(e, e_k, f):
    # e : energy
    # e_k : dispersion relation
    # f : full width half max of gaussian

    # gaussian used to approximate delta function

    ret = (2. / f) * math.sqrt(math.log(2) / math.pi) * \
        math.exp(-(4. * math.log(2.) * (e - e_k)**2.) / f**2.)
    return ret


@jit
def SC_coherence_comp_OAAT(delta, E, k_x, k_y, mu, a, t1, t2, t3):
    if E < delta:
        xi_k = 0
    else:
        xi_k = np.sqrt(E**2. - delta**2.)
    evals_k, evecs_k = SC_dispersion_comp_OAAT(delta, xi_k)
    if evals_k[0] < 0.:
        v1 = evecs_k[:, 1]
        v2 = evecs_k[:, 0]
        evecs_k[:, 0] = v1
        evecs_k[:, 1] = v2
    evals_kp, evecs_kp = SC_dispersion_comp_OAAT(
        delta, dispersion(mu, k_x, k_y, a, t1, t2, t3))
    if evals_kp[0] < 0.:
        v1 = evecs_kp[:, 1]
        v2 = evecs_kp[:, 0]
        evecs_kp[:, 0] = v1
        evecs_kp[:, 1] = v2
    S = np.array([[1., 0.], [0., -1.]])
    return (np.dot(np.transpose(evecs_kp), np.dot(S, evecs_k))[0][0])**2.


@jit
def SC_coherence_comp(delta, E, k_x, k_y, mu, a, t1, t2, t3):
    num_k = len(k_x)
    ret = np.zeros(num_k)
    for i in range(0, num_k):
        ret[i] = SC_coherence_comp_OAAT(
                delta, E, k_x[i], k_y[i], mu, a, t1, t2, t3)

    return ret


def coherence(delta, E, k_x, k_y):

    ret = 1. + (np.sqrt(E**2. - delta**2.) * dispersion(0., k_x, k_y, 1., 1., 0., 0.) -
                delta**2.) / (E * SC_dispersion_plus(delta, dispersion(0., k_x, k_y, 1., 1., 0., 0.)))

    return 0.5 * ret

# # print coherence(.3,0.4,1.,1.)
# # print SC_coherence_comp_OAAT(.3,0.4,1.,1.)


@jit
def integrator_delta_herm(integrate, dispersion, dk, E, n, f):
    # integrate : 2D array of the function we wish to integrate over k-space
    # dispersion : 2D array with all dispersion relation values over k-space
    # dk : distance between k values in our grid
    # E : energy that we're currently integrating at
    # n : number of terms in hermite polynomial expansion
    # print("At integrator")
    coeff = np.zeros(int(n + 1))
    for j in range(0, len(coeff)):
        if j % 2 == 0:
            coeff[int(j)] = ((-1)**float(j)) / \
                (math.factorial(float(j)) * 4**float(j) * np.sqrt(np.pi))
        else:
            coeff[int(j)] = 0.

    # print("after integrator loop")
    dispersion_array = np.array(dispersion)
    integrate_array = np.array(integrate)


    # print("after array cast")
    # exp = np.array(np.exp(-(((E-dispersion_array)**2.)/(2.*f**2.))))
    exp = np.array(np.exp(-(((E - dispersion_array)**2.) / (4. * f))))
    # print("after exp")
    hermite = np.array(np.polynomial.hermite.hermval(
        E - dispersion_array, coeff))
    # print("after hermite")
    # func = (1./np.sqrt(2*np.pi*f))*exp*hermite
    func = exp * hermite * (1. / (2. * np.sqrt(np.pi * f)))
    # print("after func")
    ret = func * integrate_array
    # print("after ret")
    # # print integrate_array, func, ret
    # # print np.sum(ret)
    return np.sum(ret) * dk * dk


@jit
def integrator_delta_grad(integrate, dispersion, dk, E):
    dispersion_array = np.array(dispersion)
    integrate_array = np.array(integrate)

    delta_mask = dispersion_array == E
    grad_disp_x = np.gradient(dispersion_array, dk)[1]
    grad_disp_y = np.gradient(dispersion_array, dk)[0]

    grad_disp = np.sqrt(grad_disp_x**2. + grad_disp_y**2.)
    grad_disp_delta = grad_disp[delta_mask]

    if (type(integrate) is float) or (type(integrate) is int):
        integrate_delta = integrate_array
    else:
        integrate_delta = integrate_array[delta_mask]
    # # print grad_disp_delta
    # # print integrate_delta
    if not integrate_delta.all():
        return 0.
    else:
        ret = integrate_delta / grad_disp_delta
        return np.sum(ret) * dk * dk


@jit
def BZ_integrator(function, dk):
    # takes an array over the BZ (function), and integrates it by summing the points up

    func = np.array(function)
    return np.sum(func) * dk * dk


@jit
def scattering_SC(delta, a, E, mu, t1, t2, t3, num_k, kx, ky, n, f):
    # delta : superconducting order parameter
    # a : lattice parameter
    # num_k : number of k opints along a single axis

    dk = 2. * np.pi / (a * num_k)


    Tau = integrator_delta_herm(SC_coherence_comp(delta, E, kx, ky, mu, a, t1, t2, t3), SC_dispersion_comp_new(
        delta, dispersion(mu, kx, ky, a, t1, t2, t3)), dk, E, n, f)
    # Tau = integrator_delta_grad(SC_coherence_comp(delta,E,X,Y,mu,a,t1,t2,t3),SC_dispersion_comp(delta,dispersion(mu,X,Y,a,t1,t2,t3)),dk,E)
    return np.power(Tau, -1)


@jit
def scattering_n(a, E, mu, t1, t2, t3, num_k, kx, ky, n, f):
    # scattering rate in the normal (non-superconducting) state
    # kx, ky: 1D arrays
    # print("len disp: ", len(dispersion(mu, kx, ky, a, t1, t2, t3)))
    dk = 2. * np.pi / (a * num_k)
    # print("scattering_n before Tau")
    Tau = integrator_delta_herm(1., dispersion(
        mu, kx, ky, a, t1, t2, t3), dk, E, n, f)
    # Tau = integrator_delta_grad(1.,dispersion(mu,X,Y,a,t1,t2,t3),dk,E)
    # print(np.power(Tau, -1))
    return np.power(Tau, -1)


@jit
def fermi_deriv(E, T):
    # derivative of the fermi function
    energy = np.array(E)
    ret = np.exp(energy / (k_b * T)) / ((1. + np.exp(energy / (k_b * T)))**2.)
    return ret


@jit
def FullG(gap, Temperature):
    # snagged from Sourav's code
    m = 500
    # global Temperature
    gaps = gap**2.
    f = np.log(Temperature)
    for n in range(1, m + 1):
        nh = n - 0.5
        sq2 = np.sqrt((Temperature * nh)**2. + gaps)
        f = 1. / nh - Temperature / sq2 + f
    return f


def Delta(T):
    # global Temperature
    if T < 1.0:
        return np.absolute(2. * np.pi * fsolve(FullG, [0.9], args=T))[0]
    else:
        return 0.0


@jit
def thermal_conductivity(a, num_k, n, f, num_T, T_end, mu, t1, t2, t3):

    # global Temperature

    Temp = []
    k_xx = []
    k_xy = []

    dT = T_end / num_T
    T = dT
    # Temperature = T

    k_x = np.linspace(-np.pi / a, np.pi / a, num_k)
    k_y = np.linspace(-np.pi / a, np.pi / a, num_k)
    X, Y = np.meshgrid(k_x, k_y)
    dk = 2. * np.pi / (a * num_k)

    # disp = quad_dispersion(X,Y)
    disp = dispersion(mu, X, Y, a, t1, t2, t3)

    Tau = np.zeros((len(disp), len(disp[0])))
    for i in range(0, len(disp)):
        for j in range(0, len(disp[0])):
            Tau[i][j] = scattering_n(
                a, disp[i][j], mu, t1, t2, t3, num_k, n, f)
    E_squared = np.array(disp)**2.
    # # print Tau

    # np.gradient returns 2 arrays (the x-deriv & y_deriv), so I just matched them up here
    v_x = np.gradient(disp, dk)[1]
    v_y = np.gradient(disp, dk)[0]

    for t in range(0, num_T):
        delta = Delta(T)
        disp_SC = SC_dispersion_comp(delta, disp)
        # disp_SC = SC_dispersion_plus(delta,disp)
        """
		for i in range(0,int(num_k)):
			for j in range(0,int(num_k)):
				disp_SC_comp[i][j], junkData = np.max(
				    SC_dispersion_comp_OAAT(delta,disp[i][j]))
		"""
        # # print disp_SC_comp
        # # print " "

        Tau_SC = np.zeros((len(disp_SC), len(disp_SC[0])))
        for i in range(0, len(disp_SC)):
            for j in range(0, len(disp_SC[0])):
                Tau_SC[i][j] = scattering_SC(
                    delta, a, disp_SC[i][j], mu, t1, t2, t3, num_k, n, f)
        # # print Tau_SC

        E_squared_SC = np.power(disp_SC, 2.)

        v_x_SC = np.gradient(disp_SC, dk)[1]
        v_y_SC = np.gradient(disp_SC, dk)[0]

        fermi = fermi_deriv(disp, T)
        fermi_SC = fermi_deriv(disp_SC, T)

        Temp.append(T)

        k_xx_n = BZ_integrator(
            E_squared * Tau * fermi * v_x * v_x, dk) / (T * T)
        k_xy_n = BZ_integrator(
            E_squared * Tau * fermi * v_x * v_y, dk) / (T * T)
        k_xx_SC = BZ_integrator(
            E_squared_SC * Tau_SC * fermi_SC * v_x_SC * v_x_SC, dk) / (T * T)
        k_xy_SC = BZ_integrator(
            E_squared_SC * Tau_SC * fermi_SC * v_x_SC * v_y_SC, dk) / (T * T)

        k_xx.append(k_xx_SC / k_xx_n)
        # k_xx.append(k_xx_n)
        k_xy.append(k_xy_SC / k_xx_n)
        # # print T, delta, k_xx_SC/k_xx_n, k_xy_SC/k_xx_n,k_xx_SC, k_xx_n
        T = T + dT
        # Temperature = T

    # plt.plot(Temp,k_xx)
    # plt.plot(Temp,k_xy)
    # plt.show()
    return Temp, k_xx, k_xy


@jit
def thermal_conductivity_OAAT(args):
    a, num_k, n, f, T, mu, t1, t2, t3, Tc, Ec = args

    k_x = np.linspace(-np.pi / a, np.pi / a, num_k)
    k_y = np.linspace(-np.pi / a, np.pi / a, num_k)
    X, Y = np.meshgrid(k_x, k_y)
    dk = 2. * np.pi / (a * num_k)
    # print("1")
    disp = dispersion(mu, X, Y, a, t1, t2, t3)
    E_mask = disp <= Ec
    disp_new = disp[E_mask]
    x_new = X[E_mask]
    y_new = Y[E_mask]
    # print("2")

    # # print(x_new)
    # # print("\n")
    # # print(y_new)

    Tau = np.zeros(len(disp_new))
    for i in range(0, len(disp_new)):
        Tau[i]=scattering_n(a, disp_new[i], mu, t1, t2, t3,
                            num_k, x_new, y_new, n, f)
    E_squared=np.array(disp_new)**2.

    # print("3")
    # np.gradient returns 2 arrays (the x-deriv & y_deriv), so I just matched them up here
    v_x=np.gradient(disp, dk)[1]
    v_y=np.gradient(disp, dk)[0]
    v_x_new = v_x[E_mask]
    v_y_new = v_y[E_mask]
    # print("4")
    # This is typically in the T-loop, but now we just want one T-step
    delta=Delta(T) * Tc
    disp_SC=SC_dispersion_comp_new(delta, disp_new)
    disp_SC_old = SC_dispersion_comp(delta, disp)
    # print("5")
    Tau_SC=np.zeros(len(disp_SC))
    for i in range(0, len(disp_SC)):
        Tau_SC[i]=scattering_SC(delta, a, disp_SC[i],
                                mu, t1, t2, t3, num_k, x_new, y_new, n, f)
    # print("6")
    E_squared_SC=np.power(disp_SC, 2.)

    v_x_SC=np.gradient(disp_SC_old, dk)[1]
    v_y_SC=np.gradient(disp_SC_old, dk)[0]
    v_x_SC_new = v_x_SC[E_mask]
    v_y_SC_new = v_y_SC[E_mask]
    # print("7")
    fermi=fermi_deriv(disp_new, T)
    fermi_SC=fermi_deriv(disp_SC, T)
    # print("8")
    k_xx_n=BZ_integrator(E_squared * Tau * fermi * v_x_new * v_x_new, dk) / (T * T)
    k_xy_n=BZ_integrator(E_squared * Tau * fermi * v_x_new * v_y_new, dk) / (T * T)
    k_xx_SC=BZ_integrator(E_squared_SC * Tau_SC *
                            fermi_SC * v_x_SC_new * v_x_SC_new, dk) / (T * T)
    k_xy_SC=BZ_integrator(E_squared_SC * Tau_SC *
                            fermi_SC * v_x_SC_new * v_y_SC_new, dk) / (T * T)
    # print("9")
    # print("finished fxn")
    return T, k_xx_SC / k_xx_n, k_xy_SC / k_xx_n


def main():
    # t0 = time.time()
    # T, k_xx, k_xy = thermal_conductivity(
    #     1., 10, 3, .1, 4, 1.2, 0.29, -1., 0.4, 0.)
    # t1 = time.time()
    # ## print("Time for one process: ", t1 - t0, " seconds")
    #
    T=[]
    k_xx=[]
    k_xy=[]

    pool=mp.Pool(processes=8)
    args=[(1., 16, 3., .05, .06 * i, 0., 100., 0., 0., 1., 30.)
            for i in range(1, 20)]
    # print(args)
    t0=time.time()
    results=pool.map(thermal_conductivity_OAAT, args)
    t1=time.time()
    print results
    # # print("Time for eight processes: ", t1 - t0, " seconds")
    for i in range(0, len(results)):
        T.append(results[i][0])
        k_xx.append(results[i][1])
        k_xy.append(results[i][2])
    plt.plot(T, k_xx)
    plt.plot(T, k_xy)
    plt.show()


main()
# scattering_SC(0.3,1.,2.,200,3,.1)
