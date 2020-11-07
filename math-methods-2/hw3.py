#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
# PHSX 567 homework assignment 3: Integration of Mandelbrot Set
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
# Date: ##############################################################################

import math
import warnings
import numpy as np
from scipy import integrate
from numba import jit, NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt

# QUAD integrator throws warning for low number of intervals.
# Turn them off:
# warnings.filterwarnings("ignore")

warnings.simplefilter(
    'ignore', category=NumbaDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaWarning)

MAX_IT = 1000       # Max iterations

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
        print(str(progress_string) + " " +(end_str), end='\r')
    else:
        print(progress_string, end='\r')

    if percentage > 99:
        print("")


@jit
def mandelbrot(c,maxiter):
    z = c
    for n in range(maxiter):
        if abs(z) > 2:
            return n
        z = z*z + c
    return 0



@jit
def mandelbrot_set(xmin,xmax,ymin,ymax,width,height,maxiter):
    r1 = np.linspace(xmin, xmax, width)
    r2 = np.linspace(ymin, ymax, height)
    n3 = np.empty((width,height))
    max_pts = width * height
    print("Generating grid of " + str(max_pts) + " points...")
    pt = 0
    for i in range(width):
        for j in range(height):
            pt += 1
            n3[i,j] = mandelbrot(r1[i] + 1j*r2[j],maxiter)
            progress_bar(pt, max_pts, "grid points")
    return (r1,r2,n3)


def mandelbrot_image(xmin,xmax,ymin,ymax,width=10,height=10,\
                     maxiter=80,cmap='jet',gamma=0.3):
    dpi = 72
    img_width = dpi * width
    img_height = dpi * height
    x,y,z = mandelbrot_set(xmin,xmax,ymin,ymax,img_width,img_height,maxiter)

    fig, ax = plt.subplots(figsize=(width, height),dpi=72)
    ticks = np.arange(0,img_width,3*dpi)
    x_ticks = xmin + (xmax-xmin)*ticks/img_width
    plt.xticks(ticks, x_ticks)
    y_ticks = ymin + (ymax-ymin)*ticks/img_width
    plt.yticks(ticks, y_ticks)
    ax.set_title(cmap)

    norm = colors.PowerNorm(gamma)
    ax.imshow(z.T,cmap=cmap,origin='lower',norm=norm)
    plt.show()


if __name__ == "__main__":
    # Allow this file to be run as a script
    mandelbrot_image(-2.0,2.0,-2.0,2.0, maxiter = 1000, cmap='hot')
