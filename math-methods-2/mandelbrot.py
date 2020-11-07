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
# Created by:   Ryan French
# Date:         02-15-2020 ##############################################################################

import math
from statistics import stdev
import warnings
import random
import numpy as np
from scipy import integrate
from numba import jit, NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt

warnings.simplefilter(
    'ignore', category=NumbaDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaWarning)

# Default Monte Carlo reps
REP = 1000


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

#
# Functions to check if point is in mandelbrot set
#
@jit
def in_mandelbrot(c, lem):
    z = c
    for n in range(lem):
        if abs(z) > 2:
            return 0
        z = z * z + c
    return 1

@jit
def mandelbrot(c, lem):
    z = c
    for n in range(lem):
        if abs(z) > 2:
            return n
        z = z*z + c
    return 0

@jit
def smooth_mandelbrot(c, lem, horizon, log_horizon):
    z = c
    for n in range(lem):
        az = abs(z)
        if az > horizon:
            return n - np.log(np.log(az))/np.log(2) + log_horizon
        z = z*z + c
    return 0


#Calculate area of mandelbrot set, up to a certain lemniscate

@jit
def mand_area(lem, xmin=-2.0, xmax=2.0, ymin=-2.0, ymax=2.0, rep=REP):
    rand = random.uniform
    num_in_lem = 0
    for i in range(rep):
        y = rand(xmin, xmax) + rand(ymin, ymax)*1j
        if in_mandelbrot(y, lem) == 1:
            num_in_lem += 1

    return num_in_lem * (xmax - xmin) * (ymax - ymin) / (rep)


# Create a grid and overlay numbers corresponding to values in a lemniscate, for plotting purposes

@jit
def mandelbrot_set(xmin=-2.0,xmax=2.0,ymin=-2.0,ymax=2.0,pts_x=100,pts_y=100,maxlem=100):
    r1 = np.linspace(xmin, xmax, pts_x)
    r2 = np.linspace(ymin, ymax, pts_y)
    n3 = np.empty((pts_x,pts_y))
    max_pts = pts_x * pts_y
    print("Generating grid of " + str(max_pts) + " points...")
    pt = 0
    for i in range(pts_x):
        for j in range(pts_y):
            pt += 1
            n3[i,j] = mandelbrot(r1[i] + 1j*r2[j], maxlem)
            progress_bar(pt, max_pts, "grid points")
    return (r1,r2,n3)

@jit
def smooth_mandelbrot_set(xmin=-2.0,xmax=2.0,ymin=-2.0,ymax=2.0,pts_x=100,pts_y=100,maxlem=100):
    horizon = 2.0 ** 40
    log_horizon = np.log(np.log(horizon))/np.log(2)
    r1 = np.linspace(xmin, xmax, pts_x)
    r2 = np.linspace(ymin, ymax, pts_y)
    n3 = np.empty((pts_x,pts_y))
    max_pts = pts_x * pts_y
    print("Generating grid of " + str(max_pts) + " points...")
    pt = 0
    for i in range(pts_x):
        for j in range(pts_y):
            pt += 1
            n3[i,j] = smooth_mandelbrot(r1[i] + 1j*r2[j],maxlem,horizon,log_horizon)
            progress_bar(pt, max_pts, "grid points")
    return (r1,r2,n3)



# Generate heatmap

def mandelbrot_image(xmin=-2.0,xmax=2.0,ymin=-2.0,ymax=2.0,pts_x=10,pts_y=10,maxlem=100,dpi=72):
    img_width = dpi * pts_x
    img_height = dpi * pts_y
    x,y,z = smooth_mandelbrot_set(xmin,xmax,ymin,ymax,img_width,img_height,maxlem)

    fig,ax = plt.subplots(figsize=(pts_x, pts_y), dpi=dpi)

    ticks = np.arange(0,img_width+dpi,dpi)
    x_ticks = xmin + (xmax-xmin)*ticks/img_width
    plt.xticks(ticks, x_ticks)
    y_ticks = ymin + (ymax-ymin)*ticks/img_width
    plt.yticks(ticks, y_ticks)
    plt.title(f"Mandelbrot Set, maximum lemniscas: {maxlem}, grid size: ({pts_x}, {pts_y}), dpi: {dpi}")

    light = colors.LightSource(azdeg=315, altdeg=10)
    norm = colors.PowerNorm(0.3)
    heatmap = ax.imshow(z.T, cmap="hot", norm=norm, origin="lower")
    bnds = np.arange(0, maxlem)
    plt.colorbar(heatmap, ax=ax, boundaries=bnds, spacing="uniform")
    plt.show()



# Create plot of the total area vs lemniscate

def area_plot(num_areas: int=100, rep: int=REP):
    def power_law(x, a, b, c):
        return a + b * (1/x)** c


    xs = np.arange(0, num_areas)
    ns = np.empty(num_areas)
    area = mand_area
    for i in range(num_areas):
        ns[i] = area(lem=i, rep=rep)

    rng = np.arange(0.1, num_areas, 0.1)
    plaw = [power_law(x, 1.505, 1, 1) for x in rng]

    plt.plot(xs, ns, label="Areas")
    plt.plot(rng, plaw, label="Best-fit")
    plt.title(f"Area vs lemniscate, best-fitted with 1/x")
    plt.xlabel("Lemniscate")
    plt.ylabel("Area of Mandelbrot Set")
    plt.legend()
    plt.show()



# Stratified sampling by bins

@jit
def stratified_mand_area(lem, xbins=10, ybins=10, xmin=-2.0, xmax=2.0, ymin=-2.0, ymax=2.0, rep=REP):
    xrange = xmax - xmin
    yrange = ymax - ymin
    xbinlen = xrange/xbins
    ybinlen = yrange/ybins
    xvals = np.arange(xmin, xmax + xbinlen, xbinlen)
    yvals = np.arange(ymin, ymax + ybinlen, ybinlen)

    binarea = xbinlen * ybinlen

    rand = random.uniform
    area = 0
    for i in range(len(xvals) - 1):
        for j in range(len(yvals) - 1):
            num_in_lem = 0
            for k in range(rep):
                y = rand(xvals[i], xvals[i+1]) + 1.0j*rand(yvals[j], yvals[j+1])
                if in_mandelbrot(y, lem) == 1:
                    num_in_lem += 1

            area += num_in_lem * binarea / rep

    return area

# Compare and plot std. deviations

def deviation_comps(maxlem=5000, delta_lem=100, samples=10, xbins=10, ybins=10, xmin=-2.0, xmax=2.0, ymin=-2.0, ymax=2.0, rep=REP):
    numsteps = int(maxlem/delta_lem)
    lems = np.arange(numsteps)
    stdev_orig = np.empty(numsteps)
    stdev_strat = np.empty(numsteps)

    temp_vals_orig = np.empty(samples)
    temp_vals_strat = np.empty(samples)
    lem = 0
    ar1 = mand_area
    ar2 = stratified_mand_area
    print("Generating data for StDev comparisons...")
    for i in range(numsteps):
        lems[i] = lem
        for j in range(samples):
            temp_vals_orig[j] = ar1(lem, xmin, xmax, ymin, ymax, rep)
            temp_vals_strat[j] = ar2(lem, xbins, ybins, xmin, xmax, ymin, ymax, rep)

        progress_bar(lem+delta_lem, maxlem)
        stdev_orig[i] = stdev(temp_vals_orig)
        stdev_strat[i] = stdev(temp_vals_strat)
        lem += delta_lem

    plt.scatter(lems[1:], stdev_orig[1:], c='red', label="Original")
    plt.scatter(lems[1:], stdev_strat[1:], c='blue', label="Stratified")
    plt.title(f"Standard Deviation vs lemniscate number")
    plt.legend()
    plt.xlabel("Number of lemniscates")
    plt.ylabel("Standard Deviation")
    plt.show()






if __name__ == "__main__":
    print("Homework 3, Ryan French")
    print("Note: generated plots here will have lower resolution in order to lower computational time.\n\n")


    print("Part A: Mandelbrot set image")
    mandelbrot_image(maxlem=256, pts_x=1, pts_y=1, dpi=72)
    print("\n")

    input("Press Enter to continue.")

    print("\nPart B: Area plot, unstratified")
    area_plot(num_areas=300, rep=500)
    print("\n")

    input("Press Enter to continue.")

    print("\nPart C: Comparing standard deviations")
    deviation_comps(maxlem=1000, delta_lem=100, samples=10, xbins=5, ybins=5, rep=500)
    print("\n")

    input("Press Enter to continue.")

    print("\nPart C: Calculation of Mandelbrot Set area")

    lem1 = 5000
    lem2 = 10_000

    avg1 = mand_area(lem1, rep=5_000)
    avg1_strat = stratified_mand_area(lem1, xbins=10, ybins=10, rep=5_000)
    avg2 = mand_area(lem2, rep=5_000)
    avg2_strat = stratified_mand_area(lem2, xbins=10, ybins=10, rep=5_000)

    print(f"With {lem1} lemniscates:")
    print(f"{round(avg1, 5)} (non-strat);\t{round(avg1_strat, 5)} (strat)\n")

    print(f"With {lem2} lemniscates:")
    print(f"{round(avg2, 5)} (non-strat);\t{round(avg2_strat, 5)} (strat)\n")

    print("Compare these to the theoretical area of 1.506484")
