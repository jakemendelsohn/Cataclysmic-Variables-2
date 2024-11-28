# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:07:51 2024

@author: jakem
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import lightkurve as lk
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import emcee
from lightkurve import LightCurveCollection
from astropy.io import fits
from scipy.interpolate import interp1d
from ztfquery import lightcurve
from astropy import time
from ztfquery import query
import pandas as pd
from scipy.fft import fft, fftfreq


def sector_data(index):
    search_result = lk.search_lightcurve('04 07 4.080 +18 55 37.20', mission='TESS')
    print(search_result)
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    #sap_lc_cleaned = sap_lc.remove_nans()
    quality_flags = sap_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0  # Keeps only points with a quality flag of 0 (good data)
    
    time = sap_lc.time.value[good_quality_mask]
    flux = sap_lc.flux.value[good_quality_mask]
    flux_error = sap_lc.flux_err.value[good_quality_mask]  
 
    #sap_lc_cleaned = sap_lc.remove_outliers()
    #time = sap_lc_cleaned.time.value
    #flux = sap_lc_cleaned.flux.value
    return time,flux,exptime,flux_error

def frequency_range(time,flux,del_t):
    N = len(time)
    seconds_per_day = 86400
    min_freq = 3/(N*del_t)*seconds_per_day
    max_freq = 1/(2*del_t)*seconds_per_day
    return min_freq,max_freq

def power_law_model(F, A, alpha, B):
    return A * F**alpha + B
    

def LC_model(time,flux,exptime,flux_error):
    min_freq, max_freq = frequency_range(time,flux,exptime)
    # Compute the Lomb-Scargle Periodogram within the specified frequency range
    num_frequency_points = 100000  # You can adjust this based on the desired resolution
    frequency = np.linspace(min_freq, max_freq, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency)
    
    initial_guess = [1e-2, -1, 1e-6]
    params, cov = curve_fit(power_law_model,frequency, power, p0=initial_guess)
    
    A, alpha, B = params

    plt.loglog(frequency, power, label='Data')
    plt.loglog(frequency, power_law_model(frequency, *params), label=f'Fit: A={A:.2e}, alpha={alpha:.2f}, B={B:.2e}')
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Power')
    plt.legend()
    plt.show()
    return A, alpha, B, power, frequency

def simulated_lc(power,frequency, A, alpha, B, num = 10000):
        mean = 0
        std_dev = 1
        params = A,alpha,B
        simulated_freqs = []
        reals = []
        imags = []
        for i in range(1, (num // 2) + 1):  # Only loop over positive frequencies
            random_numbers = np.random.normal(mean, std_dev, 2)
            real = random_numbers[0] * np.sqrt(power_law_model(frequency[i], *params) / 2)
            imag = random_numbers[1] * np.sqrt(power_law_model(frequency[i], *params) / 2)
            simulated_freqs.append(real + imag * 1j)

        # Handle zero-frequency and Nyquist (or last frequency) components
        simulated_freqs = [0] + simulated_freqs  # Add DC component (mean = 0)

        if num % 2 == 0:
            # Even case: add Nyquist frequency (real-only)
            simulated_freqs.append(
                np.random.normal(mean, std_dev) * np.sqrt(power_law_model(frequency[num // 2], *params))
                )

        # Add negative frequencies as conjugates of positive frequencies
        negative_freqs = np.conjugate(simulated_freqs[-2:0:-1])  # Skip DC and Nyquist
        simulated_freqs = np.concatenate((simulated_freqs, negative_freqs))
        
        # Inverse Fourier Transform to generate time series
        time_series = np.fft.ifft(simulated_freqs).real  # Ensure real time series
        
        # Plot the time series
        plt.figure(figsize=(10, 6))
        plt.plot(time_series, label="Simulated Light Curve")
        plt.title("Time Series from Inverse Fourier Transform")
        plt.xlabel("Time Steps")
        plt.ylabel("Amplitude")
        plt.legend()
        plt.show()
        
        return time_series
    
def bootstrap(time_series):
    print(time_series)
        
       
  
        
        
    
        
        
        
time,flux,exptime, flux_error = sector_data(0)
A, alpha, B,power, frequency = LC_model(time,flux,exptime, flux_error)
time_series = simulated_lc(power,frequency,A, alpha,B)
bootstrap(time_series)

