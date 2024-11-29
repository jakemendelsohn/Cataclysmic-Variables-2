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
from scipy.interpolate import interp1d


def sector_data(index):
    search_result = lk.search_lightcurve('RX J2015.6+3711', mission='TESS')
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

def simulated_lc(time,flux,power,frequency, A, alpha, B, num = 100000):
        mean_lc = np.average(flux)
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


        if num % 2 == 0:
            # Even case: add Nyquist frequency (real-only)
            simulated_freqs.append(
                np.random.normal(mean, std_dev) * np.sqrt(power_law_model(frequency[num // 2], *params))
                )

        # Add negative frequencies as conjugates of positive frequencies
        negative_freqs = np.conjugate(simulated_freqs[-2:0:-1])  # Skip DC and Nyquist
        simulated_freqs = np.concatenate((simulated_freqs, negative_freqs))
        simulated_freqs[0] = mean_lc * len(simulated_freqs) #Add DC term to make mean of simulated lc match the original mean
        
        # Inverse Fourier Transform to generate time series
        time_series = np.fft.ifft(simulated_freqs).real  # Ensure real time series
        ift_time = np.linspace(0,len(time_series),num=num)
        
        print(len(time_series))

        
        #Need to interpolate evenly spaced data to achieve same spacing as original data
        interpolator = interp1d(ift_time, time_series, kind='linear', fill_value="extrapolate")
        aligned_flux = interpolator(np.linspace(ift_time[0], ift_time[-1], len(time)))

        #plt.figure(figsize=(10, 6))
        #plt.plot(time, aligned_flux, label="Simulated Light Curve (Aligned)")
        #plt.title("Simulated Light Curve Aligned to Original Time")
        #plt.xlabel("Time")
        #plt.ylabel("Flux")
        #plt.legend()
        #plt.show()
        
        ls = LombScargle(time, aligned_flux)
        power = ls.power(frequency)
        
        aligned_time = time
        
        return aligned_time, aligned_flux,frequency, power
    
def bootstrap(N,time,flux,power, frequency, A, alpha, B):
    power_matrix = np.zeros((N, len(frequency)))
    for i in range(0,N-1):
        aligned_time, aligned_flux, frequency, new_power = simulated_lc(time, flux, power, frequency, A, alpha, B)
        power_matrix[i, :] = new_power
    
    cutoff_powers = []
    for col_idx in range(power_matrix.shape[1]):
        powers = power_matrix[:, col_idx]
        cutoff_power = np.percentile(powers, 99.7)  # 99.7% cutoff
        cutoff_powers.append(cutoff_power)
    
    cutoff_powers = np.array(cutoff_powers)
        
    plt.figure()
    plt.plot(frequency, cutoff_powers*frequency, label="99.7% Cutoff Powers")
    plt.plot(frequency, power*frequency, label = "Original Power")
    plt.xlabel("Frequency")
    plt.ylabel("Cutoff Power")
    plt.title("99.7% Cutoff Power vs Frequency")
    plt.xlim(3,10)
    plt.ylim(0,0.01)
    plt.legend()
    plt.show()
    
        
       
  
        
        
time,flux,exptime, flux_error = sector_data(0)
plt.plot(time,flux, lw=1)
plt.show()
A, alpha, B,power, frequency = LC_model(time,flux,exptime, flux_error)
aligned_time, aligned_flux,frequency, new_power = simulated_lc(time,flux,power,frequency,A, alpha,B)
N = 10
bootstrap(N,time, flux, power, frequency, A, alpha, B)

