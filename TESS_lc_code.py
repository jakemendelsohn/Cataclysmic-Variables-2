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


def sector_data(index):
    search_result = lk.search_lightcurve('TIC 15853131', mission='TESS')
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

def stitch_flatten(indeces):
    light_curves = []
    times = []
    fluxes = []
    flux_errors = []
    labels = ["sector 43", "sector 44", "sector 70", "sector 71"]   #adjust accordingly
    search_result = lk.search_lightcurve('TIC 15853131', mission='TESS')
    for i, index in enumerate(indeces):
        lc = search_result[index].download()
        processed_lc = lc.flatten()
        processed_lc = processed_lc
        sap_lc = processed_lc.SAP_FLUX
        quality_flags = sap_lc.quality
        flagged_indices = np.where(quality_flags != 0)[0]
        good_quality_mask = quality_flags == 0
        time = sap_lc.time.value[good_quality_mask]
        flux = sap_lc.flux.value[good_quality_mask]
        mean = np.mean(flux)
        flux = flux/mean
        flux_error = sap_lc.flux_err.value[good_quality_mask]  
        light_curves.append(sap_lc)
        times.append(time)
        fluxes.append(flux)
        flux_errors.append(flux_errors)
        plt.plot(time, flux,lw=1, label=labels[i])
        plt.legend()
    lc_collection = LightCurveCollection(light_curves)
    stitched_lc = lc_collection.stitch()
    
    quality_flags = stitched_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0
    
    flux_stitched = stitched_lc.flux.value[good_quality_mask]
    time_stitched = stitched_lc.time.value[good_quality_mask]
    exptime_stitched = 120 #Adjust Accordingly
    #plt.plot(time,flux,lw = 1)
    plt.xlabel("Time (BJD-2457000, days)")
    plt.ylabel("Normalised Flux")
    plt.show()
    
    return flux_stitched,time_stitched,exptime_stitched
    
def XMM_spec():
    # Load the spectrum data file
    with fits.open('C:/Users/jakem/OneDrive/Documents/Year 4 Project/XMM Data/GUEST18477658/3154/0782060201/PN/P0782060201PNS003SRSPEC0001.FTZ') as hdul:
        data = hdul[1].data  # Adjust this index if necessary
        channels = data['CHANNEL']
        counts = data['COUNTS']
        
            
    # Load the ARF file to get effective area (if available)
    with fits.open('C:/Users/jakem/OneDrive/Documents/Year 4 Project/XMM Data/GUEST18477658/3154/0782060201/PN/P0782060201PNS003SRCARF0001.FTZ') as arf_hdul:
        arf_data = arf_hdul[1].data
        effective_area = arf_data['SPECRESP']  # Effective area for each energy bin
        
 
    energy_min = 0.2  # keV
    energy_max = 10.0  # keV
    energies = np.linspace(energy_min, energy_max, len(channels))

    # Convert counts to flux: flux = counts / effective_area
    # Ensure the lengths of `counts` and `effective_area` match
    interp_func = interp1d(np.arange(len(effective_area)), effective_area, kind='linear', fill_value="extrapolate")
    effective_area_resampled = interp_func(np.linspace(0, len(effective_area) - 1, len(counts)))
    flux = counts / effective_area_resampled
    
    # Plotting the flux vs. energy spectrum
    plt.figure(figsize=(10, 6))
    plt.plot(energies, flux, drawstyle='steps-mid')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Flux (photons/cm²/s/keV)')
    plt.title('Flux vs Energy Spectrum')
    plt.show()

def ZTF_data():
    zquery = query.ZTFQuery()
    lcq = lightcurve.LCQuery.from_position(197.501495, +75.721959, 5)
    data = lcq.download_data()
    lcq = lightcurve.LCQuery(data)
    try:
        data = lcq.download_data()
        if data is not None:
            print("Data downloaded successfully")
            print(data)
        else:
            print("No data returned.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
        
def multiple_LC_plot(index_list):
    results = {}

    # Iterate over the index list and call the original function
    for idx in index_list:
        time, flux, exptime,flux_error = sector_data(idx)
        
        # Storing the 3 results as a tuple in a dictionary for easy access
        results[f"var_{idx}_1_2_3"] = (time, flux, exptime)
    times = [value[0] for value in results.values()]
    fluxes = [value[1] for value in results.values()]
    
    # Plotting
    plt.figure(figsize=(10, 6))
    
    # Use a color map to assign different colors to each sector
    colors = plt.cm.viridis(np.linspace(0, 1, len(index_list)))
    
    # Iterate through each sector and plot with a different color
    for i, (time, flux) in enumerate(zip(times, fluxes)):
        plt.plot(time, flux, lw=1, color=colors[i], label=f'Sector {index_list[i]}')
    
    plt.xlabel('Time (MJD)')
    plt.ylabel('Flux (e/s)')
    plt.title('Light Curve')
    plt.legend(loc='best', fontsize='small')  # Add a legend to identify sectors
    plt.show()
    
    
        
    
    
    

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
    frequency = np.linspace(min_freq, max_freq, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency,normalization = 'model')# Manually compute power without autopower
    # Plot the Lomb-Scargle Periodogram
    plt.figure(figsize=(10, 6))
    plt.plot(frequency, power*frequency, 'k', lw=1)
    
    # Set logarithmic scale for frequency and power if needed
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.yscale('linear')
    #plt.xscale('linear')
    plt.xlim(0,5)
    plt.ylim(0,0.05)
    
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
        
    freq_int_manual2 = [16.56576]
    
    y_vals = np.linspace(0,13,1000)
    for freq in orbital_frequencies:
        x_vals = np.linspace(freq,freq,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spin_frequencies:
        x_vals = np.linspace(freq1,freq1,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beat_frequencies:
        x_vals = np.linspace(freq2,freq2,1000)
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
        time, flux, exptime,flux_error = sector_data(idx)
        
        # Storing the 3 results as a tuple in a dictionary for easy access
        results[f"var_{idx}_1_2_3"] = (time, flux, exptime,flux_error)
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
    
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_xlim(10,30)
    ax1.set_ylim(0,0.2)
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_xlim(10,30)
    ax2.set_ylim(0,0.2)
    ax2.set_xlabel('Frequency (c/d)', fontsize = 14)
    fig.text(0, 0.5, 'Power x Frequency', va='center', rotation='vertical', fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    #Interesting_frequencies#
    freq_int_manual1 = [0.27,0.89,7.567,8.40575]
    freq_int_manual2 = [0.27,0.89,6.68,8.40575]
    #freq_int_manual1 = [6.8267,22.4769, 81.1037,95.761]
    #freq_int_manual2 = [6.8267,22.4769, 81.1037,95.761]
    
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
    for freq3 in freq_int_manual1:
        x_vals = np.linspace(freq3,freq3,1000)
        ax1.plot(x_vals,y_vals,linestyle = ":", color = 'green')
    #print("The remaining frequency peaks are", remaining_frequencies)

    
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
    for freq3 in freq_int_manual2:
        x_vals = np.linspace(freq3,freq3,1000)
        ax2.plot(x_vals,y_vals,linestyle = ":", color = 'green')
    #print("The remaining frequency peaks are", remaining_frequencies)

    
    
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

def Lomb_Scargle_2D(time, flux, exptime, window_size, step_size, min_freq, max_freq):
    # Parameters for the frequency range
    num_frequency_points = 1000  # Increase if higher frequency resolution is needed
    frequencies = np.linspace(min_freq, max_freq, num_frequency_points)

    # Prepare to store the 2D power spectrum (rows: time windows, columns: frequencies)
    power_spectrum_2D = []

    # Sliding window through the data
    time_start = min(time)
    time_end = max(time)
    current_time = time_start

    # Loop through the time windows
    while current_time + window_size <= time_end:
        # Find indices of the data within the current window
        window_mask = (time >= current_time) & (time < current_time + window_size)
        time_window = time[window_mask]
        flux_window = flux[window_mask]

        # Compute the Lomb-Scargle periodogram for the current window
        ls = LombScargle(time_window, flux_window)
        power = ls.power(frequencies)
        power_spectrum_2D.append(power)

        # Move the window
        current_time += step_size

    # Convert the list of power spectra into a 2D array
    power_spectrum_2D = np.array(power_spectrum_2D)

    # Generate the time axis (midpoints of each time window)
    time_axis = np.arange(time_start + window_size / 2, time_end, step_size)

    # Plot the 2D power spectrum
    plt.figure(figsize=(10, 6))
    plt.imshow(power_spectrum_2D.T, aspect='auto', extent=[time_axis[0], time_axis[-1], min_freq, max_freq], origin='lower', cmap='inferno')
    plt.colorbar(label='Power')
    plt.xlabel('Time (BTJD)')
    plt.ylabel('Frequency (cycles/day)')
    plt.title('2D Lomb-Scargle Power Spectrum')
    plt.show()

def peak_finder(frequency, power,  height_threshold=0.01, prominence=0.001):
    y = frequency*power
    peaks, properties = find_peaks(y, height=height_threshold, prominence=prominence)
    
    # Extract the frequencies and powers of the found peaks
    peak_frequencies = frequency[peaks]
    peak_powers = y[peaks]
    
    return peak_frequencies, peak_powers

def peak_classification(frequency,power,peak_frequencies,peak_powers,tolerance = 0.005):
    orbital_period = 0.131
    #orbital_period = 0.068233846
    spin_period = 0.1222
    natural_orbital_frequency = 1/orbital_period
    natural_spin_frequency = 1/spin_period
    natural_beat_frequency = abs(natural_spin_frequency-natural_orbital_frequency)
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
    #print("New frequencies:", new_frequencies)
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
    
    
def multi_periodic_model(t, params):
    A1, f1, phi1, A2, f2, phi2, A3, f3, phi3, A4, f4, phi4, A5, f5, phi5, offset = params
    component1 = A1 * np.sin(2 * np.pi * f1 * t + phi1)  # Spin period
    component2 = A2 * np.sin(2 * np.pi * f2 * t + phi2)  # Orbital period
    component3 = A3 * np.sin(2 * np.pi * f3 * t + phi3)  # Beat frequency
    component4 = A4 * np.sin(2* np.pi * f4 * t + phi4)
    component5 = A5 * np.sin(2* np.pi * f5 * t + phi5)
    return component1 + component2 + component3 + component4 +  offset

def chi_squared(params, time, flux, flux_error):
    """Chi-squared calculation for the model."""
    A1, f1, phi1, A2, f2, phi2, A3, f3, phi3,A4,f4,phi4,A5,f5,phi5, offset = params
    model_flux = multi_periodic_model(time, A1, f1, phi1, A2, f2, phi2, A3, f3, phi3,A4,f4,phi4,A5,f5,phi5, offset)
    print(np.sum(((flux - model_flux) / flux_error) ** 2))
    return np.sum(((flux - model_flux) / flux_error) ** 2)

def log_likelihood(params, time, flux, flux_error):
    model_flux = multi_periodic_model(time, params)
    chi_squared = np.sum(((flux - model_flux) / flux_error) ** 2)
    return -0.5 * chi_squared

def log_prior(params):
    A1, f1, phi1, A2, f2, phi2, A3, f3, phi3, A4, f4, phi4, A5, f5, phi5, offset = params
    if (
        1 <= A1 <= 10
        and 11.9 <= f1 <= 12.1
        and 0 <= phi1 <= 2 * np.pi
        and 2 <= A2 <= 10
        and 1.8 <= f2 <= 1.88
        and 0 <= phi2 <= 2 * np.pi
        and 1 <= A3 <= 10
        and 3.7 <= f3 <= 3.8
        and 0 <= phi3 <= 2 * np.pi
        and 0 <= A4 <= 10
        and 4.5 <= f4 <= 4.7
        and 0 <= phi4 <= 2 * np.pi
        and 0 <= A5 <= 10
        and 5.26 <= f5 <= 5.29
        and 0 <= phi5 <= 2 * np.pi
        and 0 <= offset <= 1250
    ):
        return 0.0  # log(1)
    return -np.inf  # log(0)

def log_probability(params, time, flux, flux_error):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, time, flux, flux_error)

def MCMC_Fit():
    start_time = 2797
    end_time = 2799
    
    subset_indices = (time >= start_time) & (time <= end_time)
    time_subset = time[subset_indices]
    flux_subset = flux[subset_indices]
    flux_error_subset = flux_error[subset_indices]
    # Initial setup for MCMC
    num_params = 16
    num_walkers = 32
    num_steps = 10000
    initial_guesses = [
        0.5,  # A1
        11.95,  # f1
        2*np.pi,  # phi1
                0.5,  # A2
                1.88,  # f2
                2*np.pi,  # phi2
                4,  # A3
                3.75,  # f3
                np.pi,  # phi3
                2,  # A4
                4.6,  # f4
                np.pi,  # phi4
                1,  # A5
                5.27,  # f5
                np.pi,  # phi5
                1175  # offset
                ]
    
    # Initialize the walkers in a Gaussian ball around initial guesses
    pos = initial_guesses + 1 * np.random.randn(num_walkers, num_params)
    
    # Set up the sampler
    sampler = emcee.EnsembleSampler(num_walkers, num_params, log_probability, args=(time_subset, flux_subset, flux_error_subset))
    
    # Run MCMC
    sampler.run_mcmc(pos, num_steps, progress=True)
    
    # Analyze the results
    samples = sampler.get_chain(discard=100, thin=15, flat=True)
    
    # Plotting the results
    import corner
    fig = corner.corner(samples, labels=["A1", "f1", "phi1", "A2", "f2", "phi2", "A3", "f3", "phi3", "A4", "f4", "phi4", "A5", "f5", "phi5", "offset"])
    plt.show()
    
    # Extract median values for each parameter
    params_mcmc = np.median(samples, axis=0)
    param_names = ["A1", "f1", "phi1", "A2", "f2", "phi2", "A3", "f3", "phi3", "A4", "f4", "phi4", "A5", "f5", "phi5", "offset"]
    params_std = np.std(samples, axis=0)  # 1-sigma uncertainty
    
    # Print the results
    print("Final parameter values and uncertainties:")
    for i, name in enumerate(param_names):
        print(f"{name}: {params_mcmc[i]:.4f} ± {params_std[i]:.4f}")
        
    
    
    # Plot the data with the best-fit model
    fitted_flux = multi_periodic_model(time_subset, params_mcmc)
    
    chi_squared = np.sum(((flux_subset - fitted_flux) / flux_error_subset) ** 2)
    reduced_chi_squared = chi_squared / (len(flux_subset) - len(params_mcmc))

    print("Reduced Chi-squared:", reduced_chi_squared)
    
    plt.figure(figsize=(12, 6))
    plt.scatter(time_subset,flux_subset, label = 'observed data', s=5,color = 'black')
    plt.plot(time_subset,flux_subset,lw=1, color = 'blue')
    plt.plot(time_subset, fitted_flux, label="MCMC Fit", color="red")
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.title("Multi-Periodic Model Fit with MCMC")
    plt.legend()
    plt.show()
        
def model_fit(time, flux, flux_error):
    start_time = 2797
    end_time = 2798
    
    subset_indices = (time >= start_time) & (time <= end_time)
    time_subset = time[subset_indices]
    flux_subset = flux[subset_indices]
    flux_error_subset = flux_error[subset_indices]


    # Initial guesses for the parameters
    initial_guesses = [
    10,  # A1: Spin amplitude
    11.997,  # f1: Spin frequency
    -np.pi,    # phi1: Spin phase
    1,   # A2: Orbital amplitude
    1.8, # f2: Orbital frequency
    0,    # phi2: Orbital phase
    5,   # A3: Beat amplitude
    3.76, # f3: Beat frequency
    np.pi,    # phi3: Beat phase
    3, #A4
    4.58, #f4
    -np.pi, #phi4
    2, #A5
    5.28, #f5
    -np.pi, #phi5
    1000  # Offset
]

    # Bounds for the parameters
    bounds = [
    (4,10),      # Bounds for A1
    (11.9,12.1),         # Bounds for f1
    (-np.pi, np.pi), # Bounds for phi1
    (0.5,2),      # Bounds for A2
    (6.5,6.5),       # Bounds for f2
    (-np.pi, np.pi), # Bounds for phi2
    (2,5),       # Bounds for A3
    (3,8),       # Bounds for f3
    (-np.pi, np.pi),# Bounds for phi3
    (2,4),
    (4.5,4.7),
    (-np.pi,np.pi),
    (1,3),
    (5.26,5.29),
    (-np.pi,np.pi),
    (800, 1200)     # Bounds for offset
]

    # Minimize the chi-squared function
    result = minimize(
        chi_squared, initial_guesses, args=(time_subset, flux_subset, flux_error_subset),
        bounds=bounds
    )

    # Extract the fitted parameters
    params = result.x
    fitted_flux = multi_periodic_model(time_subset, *params)
    chi2 = chi_squared(params, time_subset, flux_subset, flux_error_subset)
    reduced_chi2 = chi2 / (len(flux_subset) - len(params))  # Reduced chi-squared

    print("Parameters:", params)
    print("Chi-squared:", chi2)
    print("Reduced Chi-squared:", reduced_chi2)

    # Plot the fit
    plt.figure(figsize=(12, 6))
    plt.plot(time_subset, flux_subset, label='Observed Data', alpha=0.6)
    plt.scatter(time_subset, flux_subset, s=5, color = 'black')
    plt.plot(time_subset, fitted_flux, label='Fitted Model', color='red')
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.legend()
    plt.title("Multi-Periodic Model Fit with Chi-squared Minimization")
    plt.show()

    return params, result.hess_inv, fitted_flux, chi2, reduced_chi2

indexes = [0,1]
indeces = [0,1,2,3]
flux_stitched, time_stitched, exptime_stitched = stitch_flatten(indeces)
Lomb_Scargle(time_stitched, flux_stitched, exptime_stitched)
#XMM_spec()
#ZTF_data()
#multiple_LC_plot(indexes)
#times, fluxes, orbitals, spins, news = mulitple_sector_LS(indexes)
#peak_frequencies = np.array([0.27])
#phase_fold_binned(times[0], fluxes[0], peak_frequencies)


time,flux,exptime,flux_error = sector_data(indexes[1])
print(time)
#frequency,power,peak_frequencies,peak_powers = Lomb_Scargle(time,flux,exptime)
#peak_classification(frequency,power,peak_frequencies,peak_powers)
window_size = 0.3  # Window size in the same units as time (e.g., days or minutes)
step_size = 5  # Step size for sliding the window
min_freq = 5  # Minimum frequency (cycles/day)
max_freq = 40  # Maximum frequency (cycles/day)

# Call the function with your time, flux, and exptime data
#Lomb_Scargle_2D(time, flux, exptime, window_size, step_size, min_freq, max_freq)
#model_fit(time,flux,flux_error)
#MCMC_Fit()