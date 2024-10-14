# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 12:26:32 2024

@author: jakem
"""

import numpy as np
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt  
import csv
from scipy.signal import find_peaks
import scipy


def lightcurve_plot(time, flux):
    plt.scatter(time, flux, s=1)  # s=1 reduces the point size for performance
    plt.xlabel('Time')
    plt.ylabel('Flux')
    plt.show()
    
def frequency_range(time,flux,file_name):
    N = len(time)
    if file_name.endswith('LC.csv'):
        del_t = 120
    elif file_name.endswith('SC.csv'):
        del_t = 20
    else:
        raise ValueError("The file name does not end with 'LC' or 'SC'.")
    seconds_per_day = 86400
    min_freq = 3/(N*del_t)*seconds_per_day
    max_freq = 1/(2*del_t)*seconds_per_day
    return min_freq,max_freq
    
def Lomb_Scargle(time,flux,file_name):
    # Specify the desired frequency range
    min_freq, max_freq = frequency_range(time,flux,file_name)
    # Compute the Lomb-Scargle Periodogram within the specified frequency range
    num_frequency_points = 1000000  # You can adjust this based on the desired resolution
    frequency = np.linspace(min_freq, max_freq, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency)  # Manually compute power without autopower
    # Plot the Lomb-Scargle Periodogram
    plt.figure(figsize=(10, 6))
    plt.plot(frequency, power*frequency, 'k', lw=1)
    
    # Set logarithmic scale for frequency and power if needed
    plt.xscale('log')
    plt.yscale('log')
    
    # Set axis labels
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Power x Frequency')
    
    # Set plot title
    plt.title('Lomb-Scargle Periodogram (Specified Frequency Range)')
    
    fap = ls.false_alarm_probability(max(power))
    print(f"False Alarm Probability: {fap}")
    
    # Show the plot
    plt.show()
    return frequency,power   

def peak_finder(frequency, power):
    peak_idx = np.argmax(power*frequency)  # Index of the peak power
    peak_frequency = frequency[peak_idx]  # Frequency corresponding to the peak
    peak_period = 1 / peak_frequency
    print("The peak period is", peak_period*60*24, "minutes.")
    print("The peak frequency is", peak_frequency, "c/d.")
    return peak_frequency, peak_period

def fold_data(time, flux, period):
    # Fold the time data by calculating the phase (time modulo period)
    phase = np.mod(time, period)/period
    return phase, flux

def bin_folded_data(phase, flux, num_bins=20):
    bins = np.linspace(0, 1, num_bins+1)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bin_means = np.zeros(num_bins)
    bin_errors = np.zeros(num_bins)
    
    for i in range(num_bins):
        in_bin = (phase >= bins[i]) & (phase < bins[i+1])
        bin_means[i] = np.mean(flux[in_bin])
        bin_errors[i] = np.std(flux[in_bin]) / np.sqrt(np.sum(in_bin))  # Standard error
    
    return bin_centers, bin_means, bin_errors

def plot_folded_data(bin_centers, bin_means, bin_errors):
    plt.figure(figsize=(10, 6))
    # Plot the binned data with error bars
    #plt.errorbar(bin_centers, bin_means, yerr=bin_errors, fmt='k.', label='Binned Data')
    # Plot the binned data with error bars for two cycles
    plt.errorbar(np.concatenate([bin_centers, bin_centers + 1]), 
                 np.concatenate([bin_means, bin_means]), 
                 yerr=np.concatenate([bin_errors, bin_errors]), 
                 fmt='k.', label='Binned Data')
    

# Variables to store the second and third columns
time = np.array([])
flux = np.array([])
# File path to the CSV file
file_path = 'C:/Users/jakem/OneDrive/Documents/Year 4 Project/Data/TICID_350765919_Sector_34_Cadence_LC.csv'
# Open and read the CSV file
with open(file_path, 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header if there is one
    for row in reader:
        if len(row) >= 3:  # Ensure there are at least 3 columns
            time = np.append(time, float(row[1]))  # Index 1 is the second column
            flux = np.append(flux, float(row[2]))  # Index 2 is the third column

print(time)
print(scipy.__version__)
lightcurve_plot(time, flux)
frequency, power = Lomb_Scargle(time,flux,file_path)
peak_frequency, peak_period = peak_finder(frequency, power)
phase, folded_flux = fold_data(time, flux, peak_period)
bin_centers, bin_means, bin_errors = bin_folded_data(phase, folded_flux)
plot_folded_data(bin_centers,bin_means,bin_errors)

plt.figure(figsize=(10, 6))