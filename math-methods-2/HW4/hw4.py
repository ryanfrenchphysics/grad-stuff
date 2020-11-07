#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
# PHSX 567 homework assignment 4: Optimal Filtering
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
# Date: 02-28-2020 ##############################################################################

import math
import warnings
import numpy as np
from numpy.fft import fft, ifft, fftshift, fftfreq, irfft
from numpy.testing import suppress_warnings
from numpy import convolve
from scipy import integrate, signal as signal
import scipy.io as sio
from scipy.io.wavfile import read as wavread, write as wavwrite
from numba import jit, NumbaDeprecationWarning, NumbaPendingDeprecationWarning, NumbaWarning
import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt


# If we want to ignore all warnings:
# warnings.filterwarnings("ignore")
with suppress_warnings() as sup:
    sup.filter(np.ComplexWarning)

warnings.simplefilter(
    'ignore', category=NumbaDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaPendingDeprecationWarning)
warnings.simplefilter(
    'ignore', category=NumbaWarning)


music = "Lacumparsita-noisy.wav"
noise = "noise-sample.wav"


# Script globals
SIGMA = 100
PLOT_FLAG = False


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
def read_wav(in_wav):
    # returns sample_rate, data (numpy array)
    return wavread(in_wav)


@jit
def write_wav(out_wav):
    return wavwrite(filename, rate, out_wave)


@jit
def get_all_data():
    rate, data = read_wav(music)
    nrate, ndata = read_wav(noise)
    return (rate, data, nrate, ndata)


@jit
def get_power_spectra(kernel, sigma, plot_flag=PLOT_FLAG):
    r, d, nr, nd = get_all_data()
    # Generate periodogram vals for music & noise
    f, pxx = signal.periodogram(d, r, return_onesided=False)
    nf, npxx = signal.periodogram(nd, nr, return_onesided=False)

    # Do convolution:
    convolved = signal.fftconvolve(pxx, kernel, 'same')
    nconvolved = signal.fftconvolve(npxx, kernel, 'same')

    if plot_flag:
        # Plot all power spectra
        # First, switch indices so that they make sense again (FFTs in numpy switch them up)
        f_ = fftshift(f)
        pxx_ = fftshift(pxx)
        npxx_ = fftshift(npxx)
        convolved_ = fftshift(convolved)
        nconvolved_ = fftshift(nconvolved)
        fig, axs = plt.subplots(2, 2, sharex=True)
        fig.suptitle(
            "Full Power Spectra, Convolved with Gaussian (sigma = " + str(sigma) + ")")
        axs[0, 0].semilogy(f_, pxx_)
        axs[0, 0].set_title("Signal Periodogram")
        axs[0, 1].semilogy(f_, npxx_)
        axs[0, 1].set_title("Noise Periodogram")
        axs[1, 0].semilogy(f_, convolved_)
        axs[1, 0].set_title("Signal Periodogram, convolved")
        axs[1, 1].semilogy(f_, nconvolved_)
        axs[1, 1].set_title("Noise Periodogram, convolved")
        for a in axs.flat:
            a.set_ylim([1e-6, 1e10])
            a.set(xlabel='Frequency (Hz)', ylabel='PSD (V**2/Hz)')
        for a in axs.flat:
            a.label_outer()
        fig.tight_layout()
        plt.show()

        # Plot positive-valued spectra up to Nyquist Frequency
        fig, axs = plt.subplots(2, 2)
        fig.suptitle(
            "Positive Power Spectra, Convolved with Gaussian (sigma = " + str(sigma) + "), Up to f_Nyquist")
        axs[0, 0].semilogy(f_, pxx_)
        axs[0, 0].set_title("Signal Periodogram")
        axs[0, 1].semilogy(f_, npxx_)
        axs[0, 1].set_title("Noise Periodogram")
        axs[1, 0].semilogy(f_, convolved_)
        axs[1, 0].set_title("Signal Periodogram, convolved")
        axs[1, 1].semilogy(f_, nconvolved_)
        axs[1, 1].set_title("Noise Periodogram, convolved")
        for a in axs.flatten():
            a.set_ylim([1e-6, 1e10])
            a.set_xlim([0., r / 2])
            a.set(xlabel='Frequency (Hz)', ylabel='PSD (V**2/Hz)')
        for a in axs.flat:
            a.label_outer()

        fig.tight_layout()
        plt.show()

    return f, pxx, nf, npxx, convolved, nconvolved


@jit
def gaussian_smoothing_test(plot_flag=PLOT_FLAG):
    kernel = signal.gaussian(50000, 10000)

    # Dirac delta function, 500 pts, @ center of these points
    dirac = signal.unit_impulse(50000, 'mid')

    convolved = signal.fftconvolve(dirac, kernel, 'same')
    #convolved = fftshift(convolved)
    xs = np.linspace(0., 1000, 50000)
    if plot_flag:
        fig, ax = plt.subplots(sharex=True)
        ax.plot(xs[1:], convolved[1:], 'r-', label='Convolution')
        ax.plot(xs[1:], kernel[1:], 'b--', label='Original Gaussian')
        ax.plot(xs[1:], dirac[1:], 'g', label='Dirac Delta')
        leg = ax.legend()
        fig.suptitle("Gaussian Smoothing Test")
        ax.set(xlabel='unit...time? idk, you choose',
               ylabel='the dependent variable')
        ax.label_outer
        plt.show()
    return convolved


@jit
def low_pass(sample, noise, frequencies, times, plot_flag=PLOT_FLAG):
    
    filtered = np.real(np.sqrt(
        (sample - noise) /
        sample))

    if plot_flag:
        filtered_t = fftshift(ifft(filtered))
        fig, axs = plt.subplots(2)
        fig.suptitle(
            "Low Pass Filter, in f and t space")
        axs[0].plot(fftshift(frequencies), fftshift(filtered))
        axs[0].set_title("Frequency Space")
        axs[0].set(xlabel='frequency (Hz)')
        axs[1].plot(times, filtered_t)
        axs[1].set_title("Time space")
        axs[1].set(xlabel='time (s)')
        fig.tight_layout()
        plt.show()

    return filtered


@jit
def recover(low_pass, signal_data, times, offset, plot_flag=PLOT_FLAG):
    # Low pass is in frequency domain. Need data in frequency domain
    sig = fft(signal_data)
    recovered_signal = np.real(ifft(sig * low_pass))

    if plot_flag:
        signal = fftshift(recovered_signal)
        num = signal.size
        plt.plot(times, recovered_signal, 'r-')
        plt.title("Recovered Signal")
        plt.xlabel("time (s)")
        plt.ylabel("signal (V)")
        plt.show()

    return recovered_signal


@jit
def main():
    sig = 3
    offset = 0

    rate, data, nrate, ndata = get_all_data()

    timelength = data.shape[0] / rate
    times = np.linspace(0., timelength, data.shape[0])

    kernel = signal.gaussian(data.size, sig)
    f, pxx, nf, npxx, convolved, nconvolved = get_power_spectra(
        kernel, sig, PLOT_FLAG)

    # Test smoothing
    gaussian_smoothing_test(PLOT_FLAG)

    # Apply custom low-pass filter
    filtered = low_pass(convolved, nconvolved, f, times, PLOT_FLAG)

    # Recover the final signal
    recovered_signal = recover(filtered, data, times, offset, PLOT_FLAG)
    # Save this signal into a wav file
    recovered_signal.real.astype(np.int16)
    wavwrite("reconstructed.wav", rate, recovered_signal)


if __name__ == "__main__":
    main()
