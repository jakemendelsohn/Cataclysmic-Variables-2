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


def sector_data(index):
    search_result = lk.search_lightcurve('DW Cnc', mission='TESS')
    print(search_result)
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    sap_lc_cleaned = sap_lc.remove_nans()
    #sap_lc_cleaned.plot()
    #plt.figure(figsize=(12, 6))
    #plt.show()
    time = sap_lc_cleaned.time.value
    flux = sap_lc_cleaned.flux.value
    #periodogram = sap_lc_cleaned.to_periodogram(normalization='amplitude', minimum_frequency=10, maximum_frequency=1000)
    #periodogram.plot()
    #plt.show()
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
    frequency = np.linspace(min_freq*10, max_freq/5, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency,normalization = 'model')# Manually compute power without autopower
    # Plot the Lomb-Scargle Periodogram
    plt.figure(figsize=(10, 6))
    plt.plot(frequency, power*frequency, 'k', lw=1)
    
    # Set logarithmic scale for frequency and power if needed
    plt.xscale('log')
    plt.yscale('log')
    #plt.yscale('linear')
    #plt.xscale('linear')
    
    # Set axis labels
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Power x Frequency')
    
    # Set plot title
    plt.title('Lomb-Scargle Periodogram (Specified Frequency Range)')
    
    
    #Plot peaks#
    peak_frequencies, peak_powers = peak_finder(frequency, power)
    #print("Peak Frequencies", peak_frequencies)
    orbital_frequencies, spin_frequencies,new_frequencies, beat_frequencies = peak_classification(frequency,power,peak_frequencies,peak_powers)
    for freq2 in new_frequencies:
        alrm = false_alarm(ls, power,frequency, freq2)
        #print("Freq", freq2, ":", alrm)
    
    y_vals = np.linspace(0,13,1000)
    for freq in orbital_frequencies:
        x_vals = np.linspace(freq,freq,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spin_frequencies:
        x_vals = np.linspace(freq1,freq1,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beat_frequencies:
        x_vals = np.linspace(freq1,freq1,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'black')
    #print("The remaining frequency peaks are", remaining_frequencies)
    x_values = np.linspace(22.47804849975284,22.47804849975284,1000)
    plt.plot(x_values,y_vals,linestyle = ":", color = 'green')
    
    
    plt.plot([], [], linestyle=":", color='blue', label='Orbital Frequencies')  # Add one blue line to the legend
    plt.plot([], [], linestyle=":", color='red', label='Spin Frequencies')      # Add one red line to the legend
    plt.plot([], [], linestyle=":", color='green', label='Remaining Frequencies')
    plt.plot([], [], linestyle=":", color='black', label='Beat Frequencies')
    plt.legend()
    # Show the plot
    plt.show()
    return frequency,power,orbital_frequencies,spin_frequencies,new_frequencies,beat_frequencies

def mulitple_sector_LS(index_list):
    results = {}

    # Iterate over the index list and call the original function
    for idx in index_list:
        time, flux, exptime = sector_data(idx)
        
        # Storing the 3 results as a tuple in a dictionary for easy access
        results[f"var_{idx}_1_2_3"] = (time, flux, exptime)
    times = [value[0] for value in results.values()]
    fluxes = [value[1] for value in results.values()]
    exptimes = [value[2] for value in results.values()]
    
    results2 = {}
    for i in range (0,len(times)):
        frequency,power,orbital_frequencies,spin_frequencies,new_frequencies, beat_frequencies = Lomb_Scargle(times[i],fluxes[i],exptimes[i])
        results2[f"var_{i}_1_2_3_4_5_6"] = (frequency,power,orbital_frequencies,spin_frequencies,new_frequencies,beat_frequencies)

    frequencies = [value[0] for value in results2.values()]
    powers = [value[1] for value in results2.values()]
    orbitals = [value[2] for value in results2.values()]
    spins = [value[3] for value in results2.values()]
    news = [value[4] for value in results2.values()]
    beats = [value[5] for value in results2.values()]
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 12), sharex=True)
    plt.subplots_adjust(hspace=0)  # hspace=0 removes the space between the subplot
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Frequency (c/d)', fontsize = 14)
    fig.text(0, 0.5, 'Power x Frequency', va='center', rotation='vertical', fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    
    #AX1#
    ax1.plot(frequencies[0], powers[0]*frequencies[0], 'k', lw=1)
    y_vals = np.linspace(0,13,1000)
    for freq in orbitals[0]:
        x_vals = np.linspace(freq,freq,1000)
        ax1.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spins[0]:
        x_vals = np.linspace(freq1,freq1,1000)
        ax1.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beats[0]:
        x_vals = np.linspace(freq2,freq2,1000)
        ax1.plot(x_vals,y_vals,linestyle = ":", color = 'black')
    #print("The remaining frequency peaks are", remaining_frequencies)
    x_values = np.linspace(22.47804849975284,22.47804849975284,1000)
    ax1.plot(x_values,y_vals,linestyle = ":", color = 'green')
    #AX2#
    ax2.plot(frequencies[1], powers[1]*frequencies[1], 'k', lw=1)
    y_vals = np.linspace(0,13,1000)
    for freq in orbitals[1]:
        x_vals = np.linspace(freq,freq,1000)
        ax2.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spins[1]:
        x_vals = np.linspace(freq1,freq1,1000)
        ax2.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beats[1]:
        x_vals = np.linspace(freq2,freq2,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'black')
    #print("The remaining frequency peaks are", remaining_frequencies)
    x_values = np.linspace(22.47804849975284,22.47804849975284,1000)
    ax2.plot(x_values,y_vals,linestyle = ":", color = 'green')
    
    
    #Legend#
    plt.plot([], [], linestyle=":", color='blue', label='Orbital Frequencies')  # Add one blue line to the legend
    plt.plot([], [], linestyle=":", color='red', label='Spin Frequencies')      # Add one red line to the legend
    plt.plot([], [], linestyle=":", color='green', label='Remaining Frequencies')
    plt.plot([], [], linestyle=":", color='black', label='Beat Frequencies')
    plt.legend()
    
    
    # Show both plots vertically
    plt.tight_layout()
    plt.show()
    
    return times, fluxes, orbitals, spins, news

def peak_finder(frequency, power,  height_threshold=0.02, prominence=0.01):
    y = frequency*power
    peaks, properties = find_peaks(y, height=height_threshold, prominence=prominence)
    
    # Extract the frequencies and powers of the found peaks
    peak_frequencies = frequency[peaks]
    peak_powers = y[peaks]
    
    return peak_frequencies, peak_powers

def peak_classification(frequency,power,peak_frequencies,peak_powers,tolerance = 0.0015):
    orbital_period = 0.05979267
    spin_period = 0.02679429
    natural_orbital_frequency = 1/orbital_period
    natural_spin_frequency = 1/spin_period
    natural_beat_frequency = abs(natural_spin_frequency-natural_orbital_frequency)
    print("nat", natural_beat_frequency)
    orbital_frequencies = []
    spin_frequencies = []
    beat_frequencies = []
    for freq in peak_frequencies:
        if abs(freq / natural_orbital_frequency - round(freq / natural_orbital_frequency)) < tolerance:
            orbital_frequencies.append(freq)
        if abs(freq / natural_spin_frequency - round(freq / natural_spin_frequency)) < tolerance:
            spin_frequencies.append(freq)
        if abs(freq / natural_beat_frequency - round(freq / natural_beat_frequency)) < tolerance:
            beat_frequencies.append(freq)
    classified_frequencies = set(orbital_frequencies + spin_frequencies + beat_frequencies)
    
    # Filter out the frequencies that are in the classified frequencies set
    remaining_frequencies = [freq for freq in peak_frequencies if freq not in classified_frequencies]
    
    tolerance2 = 0.3  # For example, consider frequencies within 0.1 c/d as "close"
    print("spin",spin_frequencies)
    print("orbital", orbital_frequencies)
    print("beat", beat_frequencies)
    # List to store new frequencies
    new_frequencies = []
    
    # Loop through remaining frequencies
    for freq in remaining_frequencies:
        # Check if this frequency is close to any orbital or spin frequency
        is_known = np.any(np.abs(orbital_frequencies - freq) < tolerance2) or np.any(np.abs(spin_frequencies - freq) < tolerance2) or np.any(np.abs(beat_frequencies - freq) < tolerance2)
        
        # If it's not close to any known frequency, add it to new_frequencies
        if not is_known:
            new_frequencies.append(freq)
            
    # Convert new_frequencies to a numpy array for easier handling (optional)
    new_frequencies = np.array(new_frequencies)
            
            # Output the new frequencies
    print("New frequencies:", new_frequencies)
    return orbital_frequencies, spin_frequencies,new_frequencies, beat_frequencies

def false_alarm(ls,power,frequency,specific_frequency):
    closest_idx = np.argmin(np.abs(frequency - specific_frequency))
    prob = ls.false_alarm_probability(power[closest_idx])
    return prob

def bin_folded_data(phase, flux, num_bins=50):
    bins = np.linspace(0, 1, num_bins+1)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bin_means = np.zeros(num_bins)
    bin_errors = np.zeros(num_bins)
    
    for i in range(num_bins):
        in_bin = (phase >= bins[i]) & (phase < bins[i+1])
        bin_means[i] = np.mean(flux[in_bin])
        bin_errors[i] = np.std(flux[in_bin]) / np.sqrt(np.sum(in_bin))  # Standard error
    
    return bin_centers, bin_means, bin_errors

def phase_fold_binned(time, flux, peak_frequencies):
    peak_periods = 1/peak_frequencies 
    colors = ['r','g','b','y','c']
    
    for i in range(0,len(peak_periods)):
        phase = (time % peak_periods[i]) / peak_periods[i]
        bin_centers, bin_means, bin_errors = bin_folded_data(phase,flux)
        plt.errorbar(np.concatenate([bin_centers, bin_centers +1]), 
                     np.concatenate([bin_means, bin_means]), 
                     yerr=np.concatenate([bin_errors, bin_errors]), 
                     fmt='k.', label=peak_frequencies[i],color = colors[i])
    plt.legend(title = "Peak Frequencies (c/d)")
    plt.xlabel("Phase")
    plt.ylabel("Flux e/s")


indexes = [3,4]
index = 3
times, fluxes, orbitals, spins, news = mulitple_sector_LS(indexes)
peak_frequencies = np.array([21.486123909219167,14.656134673552952,22.4780485])
#phase_fold_binned(times[1], fluxes[1], peak_frequencies)
#time,flux,exptime = sector_data(index)
#frequency,power,peak_frequencies,peak_powers = Lomb_Scargle(time,flux,exptime)
#peak_classification(frequency,power,peak_frequencies,peak_powers)