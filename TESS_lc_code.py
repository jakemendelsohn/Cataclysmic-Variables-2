# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:17:32 2024

@author: jakem
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import lightkurve as lk


def sector_data():
    search_result = lk.search_lightcurve('EX Hydrae', mission='TESS')
    print(search_result)
    lc = search_result[3].download()
    exptime = search_result.table['exptime'][3]
    sap_lc = lc.SAP_FLUX
    sap_lc_cleaned = sap_lc.remove_nans()
    sap_lc_cleaned.plot()
    plt.figure(figsize=(12, 6))
    plt.show()
    time = sap_lc_cleaned.time.value
    flux = sap_lc_cleaned.flux.value
    periodogram = sap_lc_cleaned.to_periodogram(normalization='amplitude', minimum_frequency=10, maximum_frequency=1000)
    periodogram.plot()
    plt.show()
    return time,flux,exptime
    

def frequency_range(time,flux,del_t):
    N = len(time)
    seconds_per_day = 86400
    min_freq = 3/(N*del_t)*seconds_per_day
    max_freq = 1/(2*del_t)*seconds_per_day
    return min_freq,max_freq

def Lomb_Scargle(time,flux,exptime):
    # Specify the desired frequency range
    min_freq, max_freq = frequency_range(time,flux,exptime)
    # Compute the Lomb-Scargle Periodogram within the specified frequency range
    num_frequency_points = 100000  # You can adjust this based on the desired resolution
    frequency = np.linspace(min_freq*100, max_freq/5, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency)# Manually compute power without autopower
    # Plot the Lomb-Scargle Periodogram
    plt.figure(figsize=(10, 6))
    plt.plot(frequency, power*frequency, 'k', lw=1)
    
    # Set logarithmic scale for frequency and power if needed
    plt.xscale('log')
    #plt.yscale('log')
    plt.yscale('linear')
    
    # Set axis labels
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Power x Frequency')
    
    # Set plot title
    plt.title('Lomb-Scargle Periodogram (Specified Frequency Range)')
    
    
    #Plot peaks#
    peak_frequencies, peak_powers = peak_finder(frequency, power)
    y_vals = np.linspace(0,13,1000)
    for i in range(0,len(peak_powers)):
       x_vals = np.linspace(peak_frequencies[i],peak_frequencies[i],1000)
       plt.plot(x_vals,y_vals,linestyle = ':',label = peak_frequencies[i])


    #plt.legend(title = "Peak Frequencies (c/d)")
    # Show the plot
    plt.show()
    return frequency,power,peak_frequencies,peak_powers

def peak_finder(frequency, power,  height_threshold=0.05, prominence=0.0001):
    y = frequency*power
    peaks, properties = find_peaks(y, height=height_threshold, prominence=prominence)
    
    # Extract the frequencies and powers of the found peaks
    peak_frequencies = frequency[peaks]
    peak_powers = y[peaks]
    
    # Print the results
    print("Found peaks at the following frequencies (c/d) and their corresponding powers:")
    for i in range(len(peaks)):
        print(f"Peak {i + 1}: Frequency = {peak_frequencies[i]:.6f} c/d, Power = {peak_powers[i]:.6e}, Period = {1/peak_frequencies[i]:.6f} d")
        
    
    return peak_frequencies, peak_powers

def peak_classification(frequency,power,peak_frequencies,peak_powers,orbital_period = 0.06823,spin_period = 0.046546):
    natural_orbital_frequency = 1/orbital_period
    natural_spin_period = 1/spin_period
    print(natural_orbital_frequency)
    orbital_frequencies = []
    spin_frequencies = []
    
    for freq in peak_frequencies:
        # Check if it's close to the natural orbital frequency or its harmonics
        #if abs(freq - natural_orbital_frequency) < 0.5:
            #orbital_frequencies.append(freq)
        if abs(freq % natural_orbital_frequency) < 0.1:
            orbital_frequencies.append(freq)  
            
    
    print(orbital_frequencies)


time,flux,exptime = sector_data()
frequency,power,peak_frequencies,peak_powers = Lomb_Scargle(time,flux,exptime)
peak_classification(frequency,power,peak_frequencies,peak_powers)