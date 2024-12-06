# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:55:54 2024

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
    return bin_centers,bin_means
    
def phase_fold2D(time, flux, freq1, freq2, num_bins=80,smooth=False):
    freq1 = np.array([freq1])
    freq2 = np.array([freq2])
    
    period1 = 1/freq1
    period2 = 1/freq2
    

    # Calculate the phases for both frequencies
    phase1 = ((time%period1)/period1)
    phase2 = ((time%period2)/period2)
    


    phase1_binned = np.linspace(0,1,24)
    phase2_binned = np.linspace(0,1,24)
    
    print("Phase1 range:", phase1.min(), phase1.max())
    print("Phase2 range:", phase2.min(), phase2.max())
    M = np.zeros([23,23])

    for i in range(0, len(phase1_binned)-1):
        for j in range(0,len(phase2_binned)-1):
            idx = np.where(
            (phase1 > phase1_binned[i]) & (phase1 < phase1_binned[i + 1]) &
            (phase2 > phase2_binned[j]) & (phase2 < phase2_binned[j + 1])
        )[0]
            #print(idx)
            M[i,j] = np.mean(flux[idx])
    print("M",M)
    M_repeated = np.tile(M, (4, 4))  # Repeat the matrix 2x2

    # Create the plot for the repeated matrix
    plt.figure(figsize=(10, 8))
    plt.imshow(
        M_repeated, origin='lower', aspect='auto', cmap='viridis',
        extent=[0, 4, 0, 4]  # Ensure extent spans 0–2
    )
    
    # Add a color bar
    plt.colorbar(label="Flux Value")
    
    # Add labels and title
    plt.xlabel(f"Phase folded on {freq1[0]:.3f} c/d")
    plt.ylabel(f"Phase folded on {freq2[0]:.3f} c/d")
    plt.title("2D Phase Folded Flux")

    # Show the plot
    plt.show()
    
    
    
    
peak_frequencies = np.array([7.633,16.356])
    
time,flux,exptime,flux_error = sector_data(2)
phase_fold_binned(time,flux, peak_frequencies)
phase_fold2D(time,flux,peak_frequencies[0], peak_frequencies[1])





