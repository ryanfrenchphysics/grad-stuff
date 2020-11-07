##############################################################################
# PHSX 567 homework assignment 5: Chaos in Driven Pendulum
#
#
# Developed and confirmed execution on:
# OS:                 Manjaro Linux 18.1.5
# Kernel:             Linux 4..19.96-1-MANJARO
# Architecture:       x86-64
# Python Relesase:    3.7.4
#
#
#
# Created by: Ryan French
# Date: 04-03-2020 #############################################################################

# from math import *
import warnings
import numpy as np
from numpy.testing import suppress_warnings
from scipy.integrate import solve_ivp
from numba import jit, NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import nolds as nd
import seaborn as sns

# If we want to ignore all warnings:
warnings.filterwarnings("ignore")
with suppress_warnings() as sup:
    sup.filter(np.ComplexWarning)

warnings.simplefilter(
    'ignore', category=NumbaDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaWarning)


# constants
M = 0.0005  # kg
G = 9.81  # m s^-2
μ = 0.1  # s^-s
L = 0.1  # m
Q = 10**(-8)  # C
E0 = 10**6  # V m^-1
ω = 11  # rad s^-1

TIME = 30
NUM_STEPS = 2048
NUM_ω = 30
NUM_ES = NUM_ω

# Plot?
PLOT_FLAG = True

ω_E = np.array(
    [np.linspace(0.0, 30.0, NUM_ω),
    np.linspace(0.0, 1e6, NUM_ES)]
)

def print_info():
    print(f"""\
    Time length: {TIME} s
    Numer of steps: {NUM_STEPS}
    Number of ω: {NUM_ω}
    Number of E: {NUM_ES}
    """
    )

@jit
def progress_bar(current, total, end_str=""):
    """Simple Percentage Progress bar

        Args:
            current (number): The current value
            total (number): The total number of values
            end_str (str): String to follow the progress bar

        Returns:
            None
    """

    # Ignore all of this, it's a straight copy of my code that I wrote for my former company

    percentage = int(float(current / total) * 100)
    if percentage > 100:
        percentage = 100

    progress_string = '['
    hash_count = int(percentage / 2)

    for _ in range(hash_count):
        progress_string += "#"

    space_count = 50 - hash_count

    for _ in range(space_count):
        progress_string += " "

    progress_string += "] " + str(percentage) + "% "
    if end_str is not "":
        print(str(progress_string) + " " + (end_str), end='\r')
    else:
        print(progress_string, end='\r')

    if percentage > 99:
        print("")



@jit
def init_conds(ϕ0, Ω0):
    '''\
    Return init conds in a numpy array
    '''
    return np.array([ϕ0, Ω0])


@jit(nopython=True, fastmath=True)
def derivs(t, fxns, omega=ω, efield=E0):
    '''\
    Returns derivatives of the decoupled ODE to pass to ODEsolver
    '''
    # fxns[0] -> phi'(t)
    # fxns[1] -> f'(t)
    return np.array([
        fxns[1],
        ((-G / L) * np.sin(fxns[0]) -
         μ * fxns[1] + ((Q * efield) / (M * L)) * np.cos(fxns[0]) * np.cos(omega * t))
    ])


@jit(fastmath=True)
def solver(inits, time_arr, teval, omega=ω, energy=E0):
    '''\
    Solve ODE using solve_ivp method from scipy
    '''
    return solve_ivp(
        fun=derivs, t_span=tuple(time_arr), t_eval=teval, y0=inits, method='RK45', dense_output=True, args=(omega, energy)
    )


# @jit()
def lyapunov(Ω):
    '''\
    Use nolds lyap_r to find lyapunov exponent
    '''
    try:
        return nd.lyap_r(Ω)
    except:
        return 0.0


@jit
def phi_plotter(t, ϕ_stable, ϕ_unstable):
    '''\
    Plots a stable and unstable ϕ(t)
    '''
    # f is array of solutions
    if PLOT_FLAG == False:
        return
    else:
        plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(t, ϕ_stable, 'r', label=r"$\phi=0,ω=0$")
        ax.plot(t, ϕ_unstable, 'b', label=r"$\phi=10^{-4},ω=0$")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("ϕ (rad)")
        ax.legend(loc="best")
        ax.set_title(r"Driven Pendulum Amplitude, ϕ")
        plt.plot()
        plt.show()


def plot_str(omega, efield):
    return f"Ω = {omega:.2e}, E = {efield:.2e}"

# @jit(parallel=True)
def phase_plotter(num_bins=64):
    '''\
    Plot unstable phases of ϕ / ω
    '''
    if PLOT_FLAG == False:
        return
    else:
        # Clear all plots
        plt.close("all")

        t = np.linspace(0, TIME, NUM_STEPS)

        inits1 = np.array([0.0, 0.0])
        inits2 = np.array([0.0001, 0.0])

        # Array of ω - E values
        phase_vals = np.array(
        [[x,y] for x in ω_E[0] for y in ω_E[1]]
        )

        # Generate ϕ and ω from IVP
        print("Generating IVP values...\n\n")
        ϕ1 = np.array(
            [solver(inits1, (0, TIME), t, a, b).y[0] for a in ω_E[0] for b in ω_E[1]]
        )
        ϕ2 = np.array(
            [solver(inits2, (0, TIME), t, a, b).y[0] for a in ω_E[0] for b in ω_E[1]]
        )
        #
        # ϕprime1 = np.array(
        #     [solver(inits1, (0, TIME), t, a, b).y[1] for a in ω_E[0] for b in ω_E[1]]
        # )
        # ϕprime2 = np.array(
        #     [solver(inits2, (0, TIME), t, a, b).y[1] for a in ω_E[0] for b in ω_E[1]]
        # )

        # Find chaos for each ϕ(t)
        print("Determining chaos...\n\n")
        chaos1 = [lyapunov(np.array(a)) for a in list(ϕ1[:])]
        chaos2 = [lyapunov(np.array(a)) for a in list(ϕ2[:])]



        # Containers for unstable ω,E vals
        stable_vals_phi = []
        stable_vals_prime = []
        unstable_vals_phi = []
        unstable_vals_prime = []

        xmin = ω_E[0][0]
        xmax = ω_E[0][-1]
        ymin = ω_E[1][0]
        ymax = ω_E[1][-1]

        # Fill container for lyapunov > 0.0
        # Append Ω, E, if given ϕ(t)'s Lyapunov's
        # Exponent > 0
        omegas = []
        Es = []
        for i, v in enumerate(chaos1):
            if v > 0.0:
                omegas.append(phase_vals[i][0])
                Es.append(phase_vals[i][1])


        # Use PyPlot's hist2d funtion to create heatmap
        # plt.scatter(omegas, Es)#, bins=num_bins)
        # h, x, y = np.histogram2d(omegas, Es, bins=num_bins)

        # ax = sns.heatmap([omegas, Es])
        # extent= [x[0], x[-1], y[0], y[-1]]
        # plt.clf()
        plt.title("Parameter Space Plot, \"Stable\" Solution")
        plt.xlabel("Ω (rad/s)")
        plt.ylabel("E (V/m)")
        plt.imshow(np.array(chaos1).reshape((NUM_ω, NUM_ω)), interpolation='quadric', extent=[xmin, xmax, ymin, ymax], cmap=cm.hot, aspect='auto')
        plt.savefig("parameter_space1.png")


        # Fill container for lyapunov > 0.0
        omegas = []
        Es = []
        for i, v in enumerate(chaos2):
            if v > 0.0:
                omegas.append(phase_vals[i][0])
                Es.append(phase_vals[i][1])


        plt.title("Parameter Space Plot, \"Unstable\" Solution")
        plt.xlabel("Ω (rad/s)")
        plt.ylabel("E (V/m)")
        what = plt.imshow(np.array(chaos2).reshape((NUM_ω, NUM_ω)), interpolation='quadric', extent=[xmin, xmax, ymin, ymax], cmap=cm.hot, aspect='auto')
        plt.savefig("parameter_space2.png")


def main():

    print_info()

    t = np.linspace(0, TIME, NUM_STEPS)
    inits = init_conds(0, 0.0)
    inits2 = init_conds(0.0001, 0.0)

    # Solve IVP
    soln = solver(inits, (0, TIME), t)
    soln2 = solver(inits2, (0, TIME), t)

    # Get ϕ,ω values from IVP
    phi_0 = np.array(soln.y[0])
    phi_1 = np.array(soln.y[1])
    phi2_0 = np.array(soln2.y[0])
    phi2_1 = np.array(soln2.y[1])


    print("\nPlotting ϕ(t)...")
    phi_plotter(t, phi_0, phi2_0)

    print("\nGenerating phase plot...")
    phase_plotter()


if __name__ == "__main__":
    main()
