#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
# PHSX 567 homework assignment 2: Integration
#
#
# Developed and confirmed execution on:
# OS:                 Manjaro Linux 18.1.5
# Kernel:             Linux 4..19.96-1-MANJARO
# Architecture:       x86-64
# Python Relesase:    3.7.4
#
#
# SciPy Quad Integrator Documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html#scipy.integrate.quad
#
#
# Created by: Ryan French
# Date: 01/31/2020
##############################################################################

import math
import warnings
import numpy as np
from scipy import integrate
import matplotlib as mpl
import matplotlib.pyplot as plt

# QUAD integrator throws warning for low number of intervals.
# Turn them off:
warnings.filterwarnings("ignore")

N = 10000  # Number of data points for plot


def integrate_one(num_intervals: int = 50, suppress_print=False):
    """Use scipy's QUAD integrator scheme, combined with a given number of points, to integrate our integrand from 0 to 1

    Args:
        num_intervals (Optional[int]): Limit on number of intervals
        suppress_print (Optional[bool]): Whether or not to print result

    Returns:
        None
    """

    # Use QUAD integration on function, from 0 to 1, with given number of intervals.
    result = integrate.quad(
        lambda x: (1 - x)**(1. / 3) / math.log(x), 0.0, 1.0, limit=num_intervals
    )

    # The QUAD integration scheme returns a tuple: (result, error)
    val = result[0]
    err = result[1]

    result_rounded = round(val, 10)

    if suppress_print == False:
        # We want to avoid printing the result when we evaluate this integral for the error plot
        print('*'*28," PROBLEM 2 ", '*'*29)
        print("-" * 70)
        print()
        print(f"Result of integral (N=100):\t{result_rounded}")
        print(f"Result from Mathematica:\t-2.5562482320")
        print("")
        print("-" * 70)
        print("\n")

    # Return tuple, just like QUAD
    return val, err


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

    progress_string += f"]{percentage}% "
    if end_str is not "":
        print(f"{progress_string}({end_str})", end='\r')
    else:
        print(progress_string, end='\r')

    if percentage > 99:
        print("")


def error_plot() -> None:
    """Integrate the given integrand with interval values ranging from 1 to N and plot them against the error

    Args:
         None

    Returns:
        None
    """

    x = range(1, N + 1)
    y = []                  # Error value container

    print()
    print('*'*28," PROBLEM 3 ", '*'*29)
    print("-" * 70)
    print()

    print("Creating data points for error plot")
    print("Progress:")
    for i in x:
        # Generate progress bar (because there are so many points)
        # and append error values to the error container
        progress_bar(i, N, f"{i}/{N}")
        ans, err = integrate_one(i, True)
        y.append(err)

    print()
    print("You can see a steep drop in error at about N = 10. This is because the QUAD integrator is based on Fortran's QUADPACK integration scheme, which converges extremely quickly.")
    print()
    print("End of script. I hope you've enjoyed it!")
    print("")
    print("-" * 70)
    print("\n")


    # Create log-log plot
    mpl.style.use("bmh")
    plt.plot(x, y, 'r-', label='interpolated reflectivity')
    plt.title(f"Error Plot, N = {N}")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("log(N)")
    plt.ylabel("log(Error)")
    plt.show()


if __name__ == "__main__":
    # Allow this file to be run as a script
    integrate_one(100)
    error_plot()
