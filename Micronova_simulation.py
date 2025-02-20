# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:35:54 2025

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

def ASAS_SN_19bh(index):
    search_result = lk.search_lightcurve('ASASSN -19bh', mission='TESS')
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    quality_flags = sap_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0
    
    time = sap_lc.time.value[good_quality_mask]
    flux = sap_lc.flux.value[good_quality_mask]
    flux_error = sap_lc.flux_err.value[good_quality_mask]  

    return time,flux,exptime,flux_error

def ASAS_SN_delt():
    file = 'C:/Users/jakem/OneDrive/Documents/Year 4 Project/Compiled ASASSN Data/IPs_ASAS-SN_data/IPs_ASAS-SN_data/DQ_Her_ASAS-SN_LC.csv'
    df = pd.read_csv(file)
    df = df[df['mag_err'] != 99.99]
    time = df.iloc[:, 0]
    time = time-2457000
    flux = df.iloc[:, 7]
    band = df.iloc[:, 9]
    flux_error = df.iloc[:, 8]
    
    delta_array = np.diff(time)
        
    delta_sorted = np.sort(delta_array)

    n = len(delta_sorted)
    cdf = np.arange(1, n+1) / float(n)
    
    # Plot the empirical CDF
   # plt.scatter(delta_sorted, cdf,s=1)
    #plt.xlabel('Delta time')
    #plt.ylabel('CDF')
    #plt.show()
    
    return delta_sorted,cdf

def LC_model():
    time,flux,exptime,flux_error = ASAS_SN_19bh(1)
    time_mask = (time >= 2347) & (time <= 2357)
    time = time[time_mask]
    flux = flux[time_mask]
    
    factor = (2*10**(34))/np.max(flux)
    flux = flux*factor
    
    plt.scatter(time,flux,s=1,color = "black")
    
    delta_sorted, cdf = ASAS_SN_delt()
    
    time_array = np.empty(0)
    flux_array = np.empty(0)
    time_current = time[0]
    flux_current = flux[0]
    time_array = np.append(time_array,time_current)
    flux_array = np.append(flux_array,flux_current)
    while time_current < time[-1]:
        random_num = random.uniform(0, 1)
        delta_sample = np.interp(random_num, cdf, delta_sorted)
        time_current+=delta_sample
        time_array = np.append(time_array,time_current)
        idx = np.argmin(np.abs(time - time_current))
        sim_flux = flux[idx]
        flux_array = np.append(flux_array,sim_flux)
    
    plt.scatter(time_array,flux_array,color = "red", s=5)
    plt.show()
    peak_lumi = np.max(flux_array)
    
    threshold = 0.4*10**(34)
    above_threshold = np.where(flux_array > threshold)[0]
    if len(above_threshold) == 0:
        print("No points above threshold — no outburst detected.")
        peak_lumi, integral_low,integral_up,outburst_duration,duration_lower = 0,0,0,0,0
    else:
        t_start = time_array[above_threshold[0]-1]
        t_end = time_array[above_threshold[-1]+1]
        outburst_duration = t_end - t_start
        print("Duration Upper Limit", outburst_duration)
        
        t_start_lower = time_array[above_threshold[0]]
        t_end_lower = time_array[above_threshold[-1]]
        duration_lower = t_end_lower - t_start_lower
        print("Duration Lower Limit", duration_lower)
        #error1 = t_start - filtered_time_g2[above_threshold[0]-1]
        #error2 = filtered_time_g2[above_threshold[-1]+1]-t_end
        #error_tot = error1+error2
        #print("Outburst duration:", outburst_duration, "+-",error_tot, "days")
    
        start = np.where(time_array == t_start)[0]
        end = np.where(time_array == t_end)[0]
        start = start[0]
        end = end[0]
        t_burst  = time_array[start:end]
        lum_burst = flux_array[start:end]
        spd = 86400.0
        time_burst = t_burst * spd
        
        integral_up = np.trapz(lum_burst,time_burst)
        print(f"Energy Upper Limit = {integral_up:.2e} erg")
        
        start_lower = np.where(time_array == t_start_lower)[0]
        end_lower = np.where(time_array == t_end_lower)[0]
        start_lower = start_lower[0]
        end_lower = end_lower[0]
        t_burst = time_array[start_lower:end_lower]
        lum_burst = flux_array[start_lower:end_lower]
        spd = 86400.0
        time_burst = t_burst * spd
        
        integral_low = np.trapz(lum_burst,time_burst)
        print(f"Energy Lower Limit = {integral_low:.2e} erg")
        
    return peak_lumi, integral_low,integral_up,outburst_duration,duration_lower
        
        
def simulation(iterations):
    lum_array = np.empty(0)
    int_up_array = np.empty(0)
    int_low_array = np.empty(0)
    dur_up_array = np.empty(0)
    dur_low_array = np.empty(0)
    for i in range(0,iterations):       
        peak_lumi,integral_low, integral_up,duration_upper,duration_lower = LC_model()
        lum_array = np.append(lum_array,peak_lumi)
        int_up_array = np.append(int_up_array,integral_up)
        int_low_array = np.append(int_low_array,integral_low)
        dur_up_array = np.append(dur_up_array,duration_upper)
        dur_low_array = np.append(dur_low_array,duration_lower)
        
    avg_int_per_run = 0.5 * (int_up_array + int_low_array)
    avg_dur_per_run = 0.5 * (dur_up_array + dur_low_array)

    # (B) Or if you just want the **overall** average across all runs:
    mean_int_overall = np.mean(avg_int_per_run)  # average across all iterations
    mean_dur_overall = np.mean(avg_dur_per_run)
    
    plt.scatter(avg_int_per_run,lum_array)
    plt.xlabel("Average Eenergy (ergs)")
    plt.ylabel("Peak Luminosity (erg/s)")
    

    
#LC_model()
#ASAS_SN_delt()
simulation(100)