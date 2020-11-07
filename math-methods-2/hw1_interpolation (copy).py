#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
# PHSX 567 homework assignment 1: Cubic spline interpolation of CXRO reflectivity data
#
# Script usage:
# ------------------------------------------------------
# ~~~Command Line:
# python hw1_interpolation.py datafile.dat.txt spacing
#
# You may need to substitute python -> python3 if you're running othe  r versions of Linux,
# such as Ubuntu. The "spacing" argument is a float < 0.1 and is optional.
#
# inputdata.dat.txt is the text file containing the CXRO data
#
# ~~~Module:
#
# import hw1_interpolation
# csi(filename, spacing)
#
# Again, spacing is optional and filename is a string.
# ------------------------------------------------------
#
#
# Developed and confirmed execution on:
# OS:                 Manjaro Linux 18.1.5
# Kernel:             Linux 4..19.96-1-MANJARO
# Architecture:       x86-64
# Python Relesase:    3.7.4
#
#
# SciPy Cubic Spline documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html
#
#
# Created by: Ryan French
# Date: 01/21/2020
##############################################################################


import sys  # For command-line use
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib as mpl
import matplotlib.pyplot as plt


def csi(filename: str, spacing: float = 0.01) -> None:
    """Perform cubic spline interpolation on tabulated CXRO data

    Args:
        filename (str): The name of the data file containing reflectivity data
        spacing (Optional[float]): Spacing between points when plotting the interpolated function (this controls how smooth the interpolated function is)

    Returns:
        None
    """

    # Create empty containers for the data
    wavelengths = []
    reflectivities = []

    # The CXRO files aren't in any standard format, so just read line-by-line
    with open(filename) as f:
        for line in f:
            data = line.split()
            try:
                # Append data row-by-row, unless the data isn't a float
                wavelengths.append(float(data[0]))
                reflectivities.append(float(data[1]))
            except ValueError:
                pass

    if not wavelengths or not reflectivities:
        # Exit if data containers are empty, exit
        print(f"Fatal error: No data found in file {filename}! Exiting.")
        sys.exit()

    if spacing > wavelengths[1] - wavelengths[0]:
        print(
            f"Error: spacing ({spacing}) is larger than initial spacing between data ({wavelengths[1] - wavelengths[0]})! Exiting.")
        sys.exit()

    # Create a cubic spline object from data points
    cs = CubicSpline(wavelengths, reflectivities)

    # Generate a more finely-spaced array of wavelength values (allowing us to "smooth" our interpolated function)
    xs = np.arange(wavelengths[0], wavelengths[-1] + spacing, spacing)

    # Create plot
    mpl.style.use("bmh")
    plt.scatter(wavelengths, reflectivities, 75, c='k',
                marker='+', label='reflectivity data')
    plt.plot(xs, cs(xs), 'r-', label='interpolated reflectivity')
    plt.title(f"Multilayer Reflectivity Data, CXRO ({filename})")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Reflectivity")
    plt.legend(loc='upper left')
    plt.show()


if __name__ == "__main__":
    # Executes the following when run as a script rather than a module
    if len(sys.argv) == 2:
        csi(sys.argv[1])
    else:
        csi(sys.argv[1], float(sys.argv[2]))
