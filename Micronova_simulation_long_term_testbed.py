# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 13:56:39 2025

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
import os
import glob
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.interpolate import interp1d
from astroquery.vizier import Vizier
from scipy import stats
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from ztfquery import query
from ztfquery import lightcurve
import glob
import random
from sklearn.decomposition import PCA
from scipy.special import erf

def gen_scatter(average_value, time_min,time_max, time_step, scatter_std):
    time = np.arange(time_min, time_max, time_step/(3600*24))
    # Create flux array with normal scatter around the average_value
    flux = np.random.normal(loc=average_value, scale=scatter_std, size=len(time))
    return time, flux

def find_consecutive_bursts(flux, threshold, min_consecutive=3):
    above_threshold = flux > threshold
    consecutive_counts = np.convolve(above_threshold, np.ones(min_consecutive, dtype=int), mode='same')
    burst_indices = np.where(consecutive_counts >= min_consecutive)[0]
    return list(burst_indices)  # Ensure it's a list


def ASASSN_data(cutoff):
    file = 'C:/Users/jakem/OneDrive/Documents/Year 4 Project/Compiled ASASSN Data/IPs_ASAS-SN_data/IPs_ASAS-SN_data/DQ_Her_ASAS-SN_LC.csv'
    df = pd.read_csv(file)
    df = df[df['mag_err'] != 99.99]
    time = df.iloc[:, 0]
    time = time-2457000
    time_mask = time >= cutoff
    flux = df.iloc[:, 7]
    
    time = time[time_mask]
    flux = flux[time_mask]
    
    band = df.iloc[:, 9]
    flux_error = df.iloc[:, 8]
    return time,flux,flux_error
    
    
def ASAS_SN_19bh(index,name):
    search_result = lk.search_lightcurve(name, mission='TESS')
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    quality_flags = sap_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0
    
    time = sap_lc.time.value[good_quality_mask]
    flux = sap_lc.flux.value[good_quality_mask]
    flux_error = sap_lc.flux_err.value[good_quality_mask] 
    
    time_mask = (time >= 2348.5) & (time <= 2355)
    time = time[time_mask]
    flux = flux[time_mask]
    
    
    factor = (3.5*10**(34))/np.max(flux)
    flux = flux*factor

    return time,flux,exptime,flux_error

def synthesize_tess_over_asassn(tess_time, tess_flux,start_time, end_time,sequence_start):
    delta_t = 120
    mask = tess_flux <= 0.7e34
    TESS_quiescence_time = tess_time[mask]
    TESS_quiescence_flux = tess_flux[mask]
    quiescent_flux_av = np.mean(TESS_quiescence_flux)
    threshold = 1.6*quiescent_flux_av
    time,flux = gen_scatter(quiescent_flux_av,start_time,end_time,120.28,0.01e34)
    threshold_line = np.linspace(threshold,threshold,len(time))
    
    #plt.scatter(time,flux,s=1)
    #plt.scatter(tess_time,tess_flux,s=4)
    #plt.xlim(2330,2370)
    
    burst_sequence_start = sequence_start
    shifted_tess_time = tess_time - tess_time[0]

    # The total duration of the TESS light curve
    burst_time_duration = shifted_tess_time[-1] - shifted_tess_time[0]
    
    while burst_sequence_start<time[-1]:
        #mask = (time >= burst_sequence_start) & (time <= burst_sequence_start+burst_time_duration)
        #print(len(time[mask]))
        #flux[mask] = tess_flux
        #burst_sequence_start += 365.25
        #print("hi")
        start_idx = np.argmin(np.abs(time - burst_sequence_start))
        end_idx = start_idx + len(shifted_tess_time)
        if end_idx>len(time):
            burst_sequence_start+=365.25
        else:
            # Now replace exactly that slice of length len(tess_time)
            time[start_idx:end_idx] = burst_sequence_start + shifted_tess_time  
            flux[start_idx:end_idx] = tess_flux
            burst_sequence_start+=365.25
    
    
    matched_flux = np.empty(len(time_A), dtype=float)
    idxs = np.searchsorted(time, time_A)
    # Because searchsorted gives the insertion position, clip boundary
    idxs = np.clip(idxs, 1, len(time)-1)
    
    # Determine whether the left or right neighbor is closer
    left = idxs - 1
    right = idxs
    left_is_closer = np.abs(time[left] - time_A) < np.abs(time[right] - time_A)
    
    # Select whichever side is closer
    closest_idxs = np.where(left_is_closer, left, right)
    
    # Finally, get the flux at those positions
    matched_flux = flux[closest_idxs]

    # Now scatter‐plot the matched points in red
    #plt.scatter(time, flux, s=1, label="Full curve")   # the main data
    #plt.scatter(time_A, matched_flux, color="red", s=8, label="ASASSN cadence")
    #plt.plot(time,threshold_line,linestyle = "--",lw=3)
    #plt.legend(loc="upper left")
    #plt.show()
    
    
    burst_list = find_consecutive_bursts(matched_flux, threshold)
    print(burst_list)
    return burst_list

def saturating_exp(x, y0, A, k, x0):
    # y(x) = y0 + A * (1 - exp(-k*(x - x0)))
    # For x < x0, this might dip below y0 if not carefully constrained,
    # so often we keep x0 <= min(x) or handle that logic. 
    return y0 + A * (1.0 - np.exp(-k*(x - x0)))

    
    


    
time,flux,exptime,flux_error = ASAS_SN_19bh(1, "ASASSN -19bh")

cutoff_times = np.linspace(0,3300,25)
print(cutoff_times)
probabilities = np.zeros(len(cutoff_times))

for j in range(0,len(cutoff_times)):
    
    time_A,flux_A,error_A = ASASSN_data(cutoff_times[j])
    
    
    start_time = time_A.iloc[0]
    end_time = time_A.iloc[-1]
    iterations = 51
    equilibrium = start_time+365
    burst_detection = np.empty(0)
    
    for i in range(0,iterations-1):
        rand_int = random.randint(0, 365)
        sign = random.choice([-1, 1])
        sequence_start =  equilibrium + sign * rand_int
        
        # Synthesize TESS data across that entire range
        burst_list = synthesize_tess_over_asassn(
            time,
            flux,
            start_time,
            end_time,
            sequence_start
            )
        if len(burst_list)>0:
            burst_detection = np.append(burst_detection,1)
        else:
            burst_detection = np.append(burst_detection,0)
        
    count_of_ones = np.sum(burst_detection == 1)
    print("ASAS-SN detects a burst", count_of_ones, "times in", iterations-1, "runs")
    print(count_of_ones/iterations)
    probabilities[j] = count_of_ones/iterations

total_times = (end_time+-269.96982999984175)-cutoff_times

p0 = [0.0, 1.0, 0.001, total_times.min()]  # example

# 4) Fit the function to your data
popt, pcov = curve_fit(saturating_exp, total_times, probabilities, p0=p0)

# 5) Extract the best-fit parameters
y0_opt, A_opt, k_opt, x0_opt = popt
print("Fitted parameters:")
print("y0 =", y0_opt)
print("A  =", A_opt)
print("k  =", k_opt)
print("x0 =", x0_opt)

# 6) Generate a smooth curve for plotting
x_fit = np.linspace(total_times.min(), total_times.max(), 200)
y_fit = saturating_exp(x_fit, *popt)


plt.scatter(total_times,probabilities, s=10)
plt.plot(x_fit, y_fit)
plt.xlabel("Total Time (days)")
plt.ylabel("Probability of Detection")
plt.legend()
plt.show()

    


