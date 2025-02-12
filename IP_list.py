# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 16:20:39 2024

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



def read_files():
    os.chdir('C:/Users/jakem/OneDrive/Documents/Year 4 Project/Compiled ASASSN Data/IPs_ASAS-SN_data/IPs_ASAS-SN_data')
    FileList = glob.glob('*.csv')
    print(f"Found files: {FileList}")  # Ensure this prints as a list of file names
    for file in FileList:
        print("Current file:", file)
    for file in FileList:
            print(f"Processing file: {file}")
            file_path = os.path.join(os.getcwd(), file)
            df = pd.read_csv(file_path)  # Read each file
            df = df[df['mag_err'] != 99.99]
            time = df.iloc[:, 0]
            time = time-2457000
            flux = df.iloc[:, 7]
            band = df.iloc[:, 9]
            indices_V = [i for i in range(len(band)) if band.iloc[i] == 'V']
            indices_g = [h for h in range(len(band)) if band.iloc[h] == 'g']
            
            time_mask_V = []
            flux_mask_V = []
            flux_mask_V2 = []
            time_mask_V2 = []
            time_mask_g = []
            flux_mask_g = []
            flux_mask_g2 = []
            time_mask_g2 = []
            
            
            for j in range(0,len(indices_V)):
                time_mask_V.append(time.iloc[indices_V[j]])
                flux_mask_V.append(flux.iloc[indices_V[j]])
                
            indices_cutoff_V = [f for f in range(len(flux_mask_V)) if flux_mask_V[f]<95]
            
            for i in range(0,len(indices_cutoff_V)):
                time_mask_V2.append(time_mask_V[indices_cutoff_V[i]])
                flux_mask_V2.append(flux_mask_V[indices_cutoff_V[i]])
            
            for k in range(0,len(indices_g)):
                time_mask_g.append(time.iloc[indices_g[k]])
                flux_mask_g.append(flux.iloc[indices_g[k]])
                
            indices_cutoff_g = [f for f in range(len(flux_mask_g)) if flux_mask_g[f]<95]
            
            for i in range(0,len(indices_cutoff_g)):
                time_mask_g2.append(time_mask_g[indices_cutoff_g[i]])
                flux_mask_g2.append(flux_mask_g[indices_cutoff_g[i]])
                
            
            # Set thresholds based on mean and standard deviation
            mean_flux_V = np.mean(flux_mask_V2)
            std_flux_V = np.std(flux_mask_V2)
            threshold_V = mean_flux_V + 3 * std_flux_V  # Adjust threshold as needed
            
            mean_flux_g = np.mean(flux_mask_g2)
            std_flux_g = np.std(flux_mask_g2)
            threshold_g = mean_flux_g + 3 * std_flux_g  # Adjust threshold as needed
            
            if len(flux_mask_g2) > 0:
                burst_indices_g = find_consecutive_bursts(flux_mask_g2, threshold_g)
            else:
                burst_indices_g = []
            if len(flux_mask_V2) > 0:
                burst_indices_V = find_consecutive_bursts(flux_mask_V2, threshold_V)
            else:
                burst_indices_V = []
            
            if len(burst_indices_V) > 0 or len(burst_indices_g) > 0:
                print(f"Bursts found in file: {file}")
                plt.figure(figsize=(10, 6))
                plt.scatter(time_mask_V2, flux_mask_V2, label='V band', s=1)
                plt.scatter(time_mask_g2, flux_mask_g2, label='g band', s=1)
                
                # Highlight bursts
                plt.scatter([time_mask_V2[i] for i in burst_indices_V],[flux_mask_V2[i] for i in burst_indices_V],  label='V bursts',s=20)
                plt.scatter([time_mask_g2[i] for i in burst_indices_g],[flux_mask_g2[i] for i in burst_indices_g],  color='purple',label='g bursts',s=20)
                
                plt.legend()
                plt.title(f'Bursts in {file}')
                plt.xlabel('Time (BJD-2457000)')
                plt.ylabel('Flux (mJy)')
                plt.show()
            
            #plt.scatter(time_mask_V2,flux_mask_V2,s=1,label = "V band")
            #plt.scatter(time_mask_g2, flux_mask_g2, s=1, label ="g band")
            #plt.legend(title = f"{file}")
            #plt.show()
            
def find_consecutive_bursts(flux, threshold, min_consecutive=4):
    above_threshold = flux > threshold
    consecutive_counts = np.convolve(above_threshold, np.ones(min_consecutive, dtype=int), mode='same')
    burst_indices = np.where(consecutive_counts >= min_consecutive)[0]
    return list(burst_indices)  # Ensure it's a list

def specific_file():
    file = 'C:/Users/jakem/OneDrive/Documents/Year 4 Project/Compiled ASASSN Data/IPs_ASAS-SN_data/IPs_ASAS-SN_data/V4743_Sgr_ASAS-SN_LC.csv'
    df = pd.read_csv(file)
    df = df[df['mag_err'] != 99.99]
    time = df.iloc[:, 0]
    time = time-2457000
    flux = df.iloc[:, 7]
    band = df.iloc[:, 9]
    flux_error = df.iloc[:, 8]
    indices_V = [i for i in range(len(band)) if band.iloc[i] == 'V']
    indices_g = [h for h in range(len(band)) if band.iloc[h] == 'g']
    
    time_mask_V = []
    flux_mask_V = []
    flux_error_mask_V = []
    flux_mask_V2 = []
    flux_error_mask_V2 = []
    time_mask_V2 = []
    time_mask_g = []
    flux_mask_g = []
    flux_error_mask_g = []
    flux_mask_g2 = []
    flux_error_mask_g2 = []
    time_mask_g2 = []
    
    
    for j in range(0,len(indices_V)):
        time_mask_V.append(time.iloc[indices_V[j]])
        flux_mask_V.append(flux.iloc[indices_V[j]])
        flux_error_mask_V.append(flux_error.iloc[indices_V[j]])
        
    indices_cutoff_V = [f for f in range(len(flux_mask_V)) if flux_mask_V[f]<95]
    
    for i in range(0,len(indices_cutoff_V)):
        time_mask_V2.append(time_mask_V[indices_cutoff_V[i]])
        flux_mask_V2.append(flux_mask_V[indices_cutoff_V[i]])
        flux_error_mask_V2.append(flux_error_mask_V[indices_cutoff_V[i]])
    
    for k in range(0,len(indices_g)):
        time_mask_g.append(time.iloc[indices_g[k]])
        flux_mask_g.append(flux.iloc[indices_g[k]])
        flux_error_mask_g.append(flux_error.iloc[indices_g[k]])
        
    indices_cutoff_g = [f for f in range(len(flux_mask_g)) if flux_mask_g[f]<95]
    
    for i in range(0,len(indices_cutoff_g)):
        time_mask_g2.append(time_mask_g[indices_cutoff_g[i]])
        flux_mask_g2.append(flux_mask_g[indices_cutoff_g[i]])
        flux_error_mask_g2.append(flux_error_mask_g[indices_cutoff_g[i]])
    
    time_mask_V2       = np.array(time_mask_V2)
    flux_mask_V2       = np.array(flux_mask_V2)
    flux_error_mask_V2 = np.array(flux_error_mask_V2)
    
    time_mask_g2       = np.array(time_mask_g2)
    flux_mask_g2       = np.array(flux_mask_g2)
    flux_error_mask_g2 = np.array(flux_error_mask_g2)
    
    
    plt.scatter(time_mask_g2,flux_mask_g2,color = 'blue', label = "g band",s=1)
    plt.scatter(time_mask_V2, flux_mask_V2, color  = 'orange', label = "V band",s=1)
    plt.legend()
    plt.show()
    
    
    return time_mask_V2,flux_mask_V2, flux_error_mask_V2, time_mask_g2,flux_mask_g2,flux_error_mask_g2

def flux_to_luminosity(flux, distance, band):
    if band == 'V':
        band_centre = 6*10**(14)
        flux_Jy = flux*10**(-3)
        flux_erg_cm2_Hz = flux_Jy/10**(23)
        flux_erg_cm2 = flux_erg_cm2_Hz*band_centre
        lum_erg = flux_erg_cm2*distance
        
      
    if band == 'g':
        band_centre = 6.7*10**(14)
        flux_Jy = flux*10**(-3)
        flux_erg_cm2_Hz = flux_Jy/10**(23)
        flux_erg_cm2 = flux_erg_cm2_Hz*band_centre
        lum_erg = flux_erg_cm2*distance
        
    
    return lum_erg
        
    
def burst_focus(time_mask_V2,flux_mask_V2, flux_error_mask_V2, time_mask_g2,flux_mask_g2,flux_error_mask_g2,distance,band):
    t_lower = 2300
    t_upper = 2500
    
    lum_erg_total_g = flux_to_luminosity(flux_mask_g2, distance, 'g')
    lum_erg_total_V= flux_to_luminosity(flux_mask_V2, distance, 'V')
    
    # Create masks (True/False arrays) indicating which points fall within [t_lower, t_upper]
    mask_V = (time_mask_V2 >= t_lower) & (time_mask_V2 <= t_upper)
    mask_g = (time_mask_g2 >= t_lower) & (time_mask_g2 <= t_upper)
    
    # Apply these masks to time, flux, and flux-error arrays
    filtered_time_V2       = time_mask_V2[mask_V]
    filtered_flux_V2       = flux_mask_V2[mask_V]
    filtered_flux_err_V2   = flux_error_mask_V2[mask_V]
    
    filtered_time_g2       = time_mask_g2[mask_g]
    filtered_flux_g2       = flux_mask_g2[mask_g]
    filtered_flux_err_g2   = flux_error_mask_g2[mask_g]
    
    
    lum_erg_g = flux_to_luminosity(filtered_flux_g2, distance, 'g')
    lum_err_g = flux_to_luminosity(filtered_flux_err_g2, distance, 'g')
    lum_erg_V = flux_to_luminosity(filtered_flux_V2, distance, 'V')
    lum_err_V = flux_to_luminosity(filtered_flux_err_V2, distance, 'V')
    
    

    
    average_g = np.mean(lum_erg_total_g)
    average_V = np.mean(lum_erg_total_V)
    
    
    if band == 'g':
    
        x_average_g = np.linspace(min(filtered_time_g2),max(filtered_time_g2),500)
        #x_average_V = np.linspace(min(filtered_time_V2),max(filtered_time_V2),500)
        y_average_g = np.linspace(average_g,average_g,500)
        #y_average_V = np.linspace(average_V,average_V,500)
        
        
        
        threshold = average_g*1.5
        print("Threshold is", threshold)
        threshold_line = np.linspace(threshold,threshold,len(filtered_time_g2))
        
        above_threshold = np.where(lum_erg_g > threshold)[0]
        if len(above_threshold) == 0:
            print("No points above threshold — no outburst detected.")
        else:
            t_start = filtered_time_g2[above_threshold[0]-1]
            t_end = filtered_time_g2[above_threshold[-1]+1]
            outburst_duration = t_end - t_start
            print("Duration Upper Limit", outburst_duration)
            
            t_start_lower = filtered_time_g2[above_threshold[0]]
            t_end_lower = filtered_time_g2[above_threshold[-1]]
            duration_lower = t_end_lower - t_start_lower
            print("Duration Lower Limit", duration_lower)
            #error1 = t_start - filtered_time_g2[above_threshold[0]-1]
            #error2 = filtered_time_g2[above_threshold[-1]+1]-t_end
            #error_tot = error1+error2
            #print("Outburst duration:", outburst_duration, "+-",error_tot, "days")
        
        start = np.where(filtered_time_g2 == t_start)[0]
        end = np.where(filtered_time_g2 == t_end)[0]
        start = start[0]
        end = end[0]
        t_burst_g = filtered_time_g2[start:end]
        lum_burst_g = lum_erg_g[start:end]
        lum_err_burst_g = lum_err_g[start:end]
        spd = 86400.0
        time_burst_g = t_burst_g * spd
    
        N_sims = 10000
        integrated_values = np.zeros(N_sims)
        lum_err_burst_g = lum_err_burst_g.astype(float)
        
        
        for i in range(N_sims):
            # draw a noisy realization of the luminosities
            L_noisy = (lum_burst_g-average_g) + np.random.normal(loc=0.0, scale=lum_err_burst_g)
            # integrate using trapezoid rule
            integrated_values[i] = np.trapz(L_noisy, time_burst_g)   
        # Mean integrated energy:
        E_mean = np.mean(integrated_values)
        # 1-sigma error (standard deviation of the distribution):
        E_std  = np.std(integrated_values)
        
            
        print(f"Monte Carlo mean energy Upper Limit = {E_mean:.2e} erg ± {E_std:.2e} erg (1σ)")
        
        start_lower = np.where(filtered_time_g2 == t_start_lower)[0]
        end_lower = np.where(filtered_time_g2 == t_end_lower)[0]
        start_lower = start_lower[0]
        end_lower = end_lower[0]
        t_burst_g = filtered_time_g2[start_lower:end_lower]
        lum_burst_g = lum_erg_g[start_lower:end_lower]
        lum_err_burst_g = lum_err_g[start_lower:end_lower]
        spd = 86400.0
        time_burst_g = t_burst_g * spd
    
        N_sims = 10000
        integrated_values = np.zeros(N_sims)
        lum_err_burst_g = lum_err_burst_g.astype(float)
        
        
        for i in range(N_sims):
            # draw a noisy realization of the luminosities
            L_noisy = (lum_burst_g-average_g) + np.random.normal(loc=0.0, scale=lum_err_burst_g)
            # integrate using trapezoid rule
            integrated_values[i] = np.trapz(L_noisy, time_burst_g)
            
        # Mean integrated energy:
        E_mean = np.mean(integrated_values)
        # 1-sigma error (standard deviation of the distribution):
        E_std  = np.std(integrated_values)
            
        print(f"Monte Carlo mean energy Lower Limit = {E_mean:.2e} erg ± {E_std:.2e} erg (1σ)")
        
    
    #start = np.where(filtered_time_V2 == t_start)[0]
    #end = np.where(filtered_time_V2 == t_end)[0]
    #start = start[0]
    #end = end[0]
    #t_burst_g = filtered_time_V2[start:end]
    #lum_burst_g = lum_erg_V[start:end]
    #spd = 86400.0
    #time_burst_g = t_burst_V * spd
    #total_energy = np.trapz(lum_burst_V, time_burst_V)
    #print(f"Total energy of the outburst = {total_energy:.2e} erg")
    
        maximum_g = np.where(lum_erg_g == np.max(lum_erg_g))
        error_mJy = filtered_flux_err_g2[maximum_g]
        error_erg = flux_to_luminosity(error_mJy, distance, 'g')
        #maximum_V = np.where(np.max(lum_erg_V))
        print("The peak luminosity is", np.max(lum_erg_g), "+-", error_erg)
        #print("The peak luminosity is", np.max(lum_erg_V))
        
    
    
        plt.scatter(filtered_time_g2,lum_erg_g,label = "g band",s=1)
        #plt.scatter(filtered_time_V2,lum_erg_V, label = "V band",s=1)
        plt.plot(x_average_g,y_average_g,linestyle = '--',color = 'orange',label = "quiescence")
        plt.plot(filtered_time_g2,threshold_line, label = "Threshold")
        plt.xlabel('Time (BJD-2457000)')
        plt.ylabel('Luminosity (erg/s)')
        #plt.xticks(np.arange(2575, 2601, 2))
        plt.legend()
        plt.show()
    
    if band == "V":
        x_average_V = np.linspace(min(filtered_time_V2),max(filtered_time_V2),500)
        y_average_V = np.linspace(average_V,average_V,500)
        threshold = average_V*1.2
        print("Threshold is", threshold)
        threshold_line = np.linspace(threshold,threshold,len(filtered_time_V2))

        
        above_threshold = np.where(lum_erg_V > threshold)[0]
        if len(above_threshold) == 0:
            print("No points above threshold — no outburst detected.")
        else:
            t_start = filtered_time_V2[above_threshold[0]-1]
            t_end = filtered_time_V2[above_threshold[-1]+1]
            outburst_duration = t_end - t_start
            print("Duration Upper Limit", outburst_duration)
            
            t_start_lower = filtered_time_V2[above_threshold[0]]
            t_end_lower = filtered_time_V2[above_threshold[-1]]
            duration_lower = t_end_lower - t_start_lower
            print("Duration Lower Limit", duration_lower)
            #error1 = t_start - filtered_time_V2[above_threshold[0]-1]
            #error2 = filtered_time_V2[above_threshold[-1]+1]-t_end
            #error_tot = error1+error2
            #print("Outburst duration:", outburst_duration, "+-",error_tot, "days")
        
        start = np.where(filtered_time_V2 == t_start)[0]
        end = np.where(filtered_time_V2 == t_end)[0]
        start = start[0]
        end = end[0]
        t_burst_V = filtered_time_V2[start:end]
        lum_burst_V = lum_erg_V[start:end]
        lum_err_burst_V = lum_err_V[start:end]
        spd = 86400.0
        time_burst_V = t_burst_V * spd
    
        N_sims = 10000
        integrated_values = np.zeros(N_sims)
        lum_err_burst_V = lum_err_burst_V.astype(float)
        
        for i in range(N_sims):
            L_noisy = (lum_burst_V-average_V) + np.random.normal(loc=0.0, scale=lum_err_burst_V)
            integrated_values[i] = np.trapz(L_noisy, time_burst_V)
            
        E_mean = np.mean(integrated_values)
        E_std  = np.std(integrated_values)
            
        print(f"Monte Carlo mean energy = {E_mean:.2e} erg ± {E_std:.2e} erg (1σ)")
        
        start_lower = np.where(filtered_time_V2 == t_start_lower)[0]
        end_lower = np.where(filtered_time_V2 == t_end_lower)[0]
        start_lower = start_lower[0]
        end_lower = end_lower[0]
        t_burst_V = filtered_time_V2[start_lower:end_lower]
        lum_burst_V = lum_erg_V[start_lower:end_lower]
        lum_err_burst_V = lum_err_V[start_lower:end_lower]
        spd = 86400.0
        time_burst_V = t_burst_V * spd
    
        N_sims = 10000
        integrated_values = np.zeros(N_sims)
        lum_err_burst_V = lum_err_burst_V.astype(float)
        
        
        for i in range(N_sims):
            # draw a noisy realization of the luminosities
            L_noisy = (lum_burst_V-average_V) + np.random.normal(loc=0.0, scale=lum_err_burst_V)
            # integrate using trapezoid rule
            integrated_values[i] = np.trapz(L_noisy, time_burst_V)
            
        # Mean integrated energy:
        E_mean = np.mean(integrated_values)
        # 1-sigma error (standard deviation of the distribution):
        E_std  = np.std(integrated_values)
            
        print(f"Monte Carlo mean energy Lower Limit = {E_mean:.2e} erg ± {E_std:.2e} erg (1σ)")
        
        maximum_V = np.where(lum_erg_V == np.max(lum_erg_V))
        error_mJy = filtered_flux_err_V2[maximum_V]
        error_erg = flux_to_luminosity(error_mJy, distance, 'V')
        #maximum_V = np.where(np.max(lum_erg_V))
        print("The peak luminosity is", np.max(lum_erg_V), "+-", error_erg)
        #print("The peak luminosity is", np.max(lum_erg_V))
    
    
    
        plt.scatter(filtered_time_V2,lum_erg_V,label = "V band",s=1)
        #plt.scatter(filtered_time_V2,lum_erg_V, label = "V band",s=1)
        plt.plot(x_average_V,y_average_V,linestyle = '--',color = 'orange',label = "quiescence")
        plt.plot(filtered_time_V2,threshold_line, label = "Threshold")
        plt.xlabel('Time (BJD-2457000)')
        plt.ylabel('Luminosity (erg/s)')
        #plt.xticks(np.arange(2575, 2601, 2))
        plt.legend()
        plt.show()
        
def literature_data():
    df = pd.read_csv('C:/Users/jakem/OneDrive/Documents/Year 4 Project/Classified Optical Outbursts Data.csv')

    # Each column can be converted into a NumPy array like so:
    names = df['Name'].to_numpy()
    types = df['Type'].to_numpy()
    peak_luminosity = df['Peak Luminosity (erg/s)'].to_numpy()
    peak_luminosity_err = df['Peak Luminosity Error (erg/s)'].to_numpy()
    total_energy = df['Total Energy (erg)'].to_numpy()
    total_energy_err = df['Total Energy Error (erg)'].to_numpy()
    duration = df['Duration (d)'].to_numpy()
    duration_err = df['Duration Error (d)'].to_numpy()
    print('Names:', names)
    print('Peak luminosities:', peak_luminosity)
    
    df_dwarf = df[df["Type"] == "Dwarf nova"]
    df_micro = df[df["Type"] == "Micronova"]
    df_donor = df[df["Type"] == "Donor flare"]
    df_gating = df[df["Type"] == "Magnetic Gating"]
    
    # Then extract arrays from each subset. For example:
    names_dwarf = df_dwarf["Name"].to_numpy()
    peak_dwarf = df_dwarf["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_dwarf = df_dwarf['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_dwarf = df_dwarf['Total Energy (erg)'].to_numpy()
    total_en_error_dwarf = df_dwarf['Total Energy Error (erg)'].to_numpy()
    duration_dwarf = df_dwarf['Duration (d)'].to_numpy()
    duration_error_dwarf = df_dwarf['Duration Error (d)'].to_numpy()
    

    
    names_micro = df_micro["Name"].to_numpy()
    peak_micro = df_micro["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_micro = df_micro['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_micro = df_micro['Total Energy (erg)'].to_numpy()
    total_en_error_micro = df_micro['Total Energy Error (erg)'].to_numpy()
    duration_micro = df_micro['Duration (d)'].to_numpy()
    duration_error_micro = df_micro['Duration Error (d)'].to_numpy()

    
    
    names_donor = df_donor["Name"].to_numpy()
    peak_donor = df_donor["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_donor = df_donor['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_donor= df_donor['Total Energy (erg)'].to_numpy()
    total_en_error_donor = df_donor['Total Energy Error (erg)'].to_numpy()
    duration_donor = df_donor['Duration (d)'].to_numpy()
    duration_error_donor = df_donor['Duration Error (d)'].to_numpy()
    
    names_gating = df_gating["Name"].to_numpy()
    peak_gating = df_gating["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_gating = df_gating['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_gating = df_gating['Total Energy (erg)'].to_numpy()
    total_en_error_gating = df_gating['Total Energy Error (erg)'].to_numpy()
    duration_gating = df_gating['Duration (d)'].to_numpy()
    duration_error_gating = df_gating['Duration Error (d)'].to_numpy()
    
    plt.scatter(total_en_micro,peak_micro,label = "Micronovae",color = "blue",s=5)
    plt.scatter(total_en_dwarf,peak_dwarf,label = "Dwarf Novae", color = "red",s=5)
    plt.scatter(total_en_gating,peak_gating, label = "Magnetic Gating", color = "yellow", s=5)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    return (names_micro, peak_micro, peak_error_micro, total_en_micro, total_en_error_micro, duration_micro, duration_error_micro,
        names_donor, peak_donor, peak_error_donor, total_en_donor, total_en_error_donor, duration_donor, duration_error_donor,
        names_gating, peak_gating, peak_error_gating, total_en_gating, total_en_error_gating, duration_gating, duration_error_gating, names_dwarf,peak_dwarf,peak_error_dwarf, total_en_dwarf,total_en_error_dwarf,duration_dwarf,duration_error_dwarf)
    

def new_data():
    df = pd.read_csv('C:/Users/jakem/OneDrive/Documents/Year 4 Project/Outbursts Compiled List (updated durations).xlsx.csv')

    # Each column can be converted into a NumPy array like so:
    names = df['Name'].to_numpy()
    types = df['Type'].to_numpy()
    peak_luminosity = df['Peak Luminosity (erg/s)'].to_numpy()
    peak_luminosity_err = df['Peak Luminosity Error (erg/s)'].to_numpy()
    total_energy_lower = df['Total Energy Lower (erg)'].to_numpy()
    total_energy_upper = df['Total Energy Upper (erg)'].to_numpy()
    duration_lower = df['Duration Lower (d)'].to_numpy()
    duration_upper = df['Duration Upper (d)'].to_numpy()

    
    df_dwarf = df[df["Type"] == "Dwarf nova"]
    df_micro = df[df["Type"] == "Micronova"]
    df_donor = df[df["Type"] == "Donor flare"]
    df_gating = df[df["Type"] == "Magnetic Gating"]
    
    df_dwarf_iron = df_dwarf[df_dwarf['Ironclad?'] == 'Y']
    df_micro_iron = df_micro[df_micro['Ironclad?'] == 'Y']
    
    # Then extract arrays from each subset. For example:
    names_dwarf = df_dwarf["Name"].to_numpy()
    peak_dwarf = df_dwarf["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_dwarf = df_dwarf['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_dwarf_low = df_dwarf['Total Energy Lower (erg)'].to_numpy()
    total_en_dwarf_up = df_dwarf['Total Energy Upper (erg)'].to_numpy()
    duration_dwarf_low = df_dwarf['Duration Lower (d)'].to_numpy()
    duration_dwarf_up = df_dwarf['Duration Upper (d)'].to_numpy()
    P_orb_dwarf = df_dwarf['P_orb (h)'].to_numpy()
    P_spin_dwarf = df_dwarf['P_spin (s)'].to_numpy()

    
    names_micro = df_micro["Name"].to_numpy()
    peak_micro = df_micro["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_micro = df_micro['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_micro_up = df_micro['Total Energy Upper (erg)'].to_numpy()
    total_en_micro_low = df_micro['Total Energy Lower (erg)'].to_numpy()
    duration_micro_up = df_micro['Duration Upper (d)'].to_numpy()
    duration_micro_low = df_micro['Duration Lower (d)'].to_numpy()
    P_orb_micro = df_micro['P_orb (h)'].to_numpy()
    P_spin_micro = df_micro['P_spin (s)'].to_numpy()
    
    names_donor = df_donor["Name"].to_numpy()
    peak_donor = df_donor["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_donor = df_donor['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_donor= df_donor['Total Energy Upper (erg)'].to_numpy()
    total_en_error_donor = df_donor['Total Energy Lower (erg)'].to_numpy()
    duration_donor = df_donor['Duration Upper (d)'].to_numpy()
    duration_error_donor = df_donor['Duration Lower (d)'].to_numpy()
    
    names_gating = df_gating["Name"].to_numpy()
    peak_gating = df_gating["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_gating = df_gating['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_gating = df_gating['Total Energy Upper (erg)'].to_numpy()
    total_en_error_gating = df_gating['Total Energy Lower (erg)'].to_numpy()
    duration_gating = df_gating['Duration Upper (d)'].to_numpy()
    duration_error_gating = df_gating['Duration Lower (d)'].to_numpy()
    P_orb_gating = df_gating['P_orb (h)'].to_numpy()
    P_spin_gating = df_gating['P_spin (s)'].to_numpy()
    
    names_dwarf_iron = df_dwarf_iron["Name"].to_numpy()
    peak_dwarf_iron = df_dwarf_iron["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_dwarf_iron = df_dwarf_iron['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_dwarf_iron_up = df_dwarf_iron['Total Energy Upper (erg)'].to_numpy()
    total_en_dwarf_iron_low = df_dwarf_iron['Total Energy Lower (erg)'].to_numpy()
    duration_dwarf_iron_up = df_dwarf_iron['Duration Upper (d)'].to_numpy()
    duration_dwarf_iron_low = df_dwarf_iron['Duration Lower (d)'].to_numpy()
    P_orb_dwarf_iron = df_dwarf_iron['P_orb (h)'].to_numpy()
    P_spin_dwarf_iron = df_dwarf_iron['P_spin (s)'].to_numpy()
    
    names_micro_iron = df_micro_iron["Name"].to_numpy()
    peak_micro_iron = df_micro_iron["Peak Luminosity (erg/s)"].to_numpy()
    peak_error_micro_iron = df_micro_iron['Peak Luminosity Error (erg/s)'].to_numpy()
    total_en_micro_iron_up = df_micro_iron['Total Energy Upper (erg)'].to_numpy()
    total_en_micro_iron_low = df_micro_iron['Total Energy Lower (erg)'].to_numpy()
    duration_micro_iron_up = df_micro_iron['Duration Upper (d)'].to_numpy()
    duration_micro_iron_low = df_micro_iron['Duration Lower (d)'].to_numpy()
    P_orb_micro_iron = df_micro_iron['P_orb (h)'].to_numpy()
    P_spin_micro_iron = df_micro_iron['P_spin (s)'].to_numpy()
    
    
    
    return (names_micro, peak_micro, peak_error_micro, total_en_micro_up, total_en_micro_low, duration_micro_up, duration_micro_low, P_orb_micro, P_spin_micro,
        names_donor, peak_donor, peak_error_donor, total_en_donor, total_en_error_donor, duration_donor, duration_error_donor,
        names_gating, peak_gating, peak_error_gating, total_en_gating, total_en_error_gating, duration_gating, duration_error_gating, P_orb_gating, P_spin_gating,names_dwarf,peak_dwarf,peak_error_dwarf, total_en_dwarf_up,total_en_dwarf_low,duration_dwarf_up,duration_dwarf_low, P_orb_dwarf, P_spin_dwarf,
        names_micro_iron, peak_micro_iron, peak_error_micro_iron, total_en_micro_iron_up, total_en_micro_iron_low, duration_micro_iron_up, duration_micro_iron_low, P_orb_micro_iron, P_spin_micro_iron,names_dwarf_iron,peak_dwarf_iron,peak_error_dwarf_iron, total_en_dwarf_iron_up,total_en_dwarf_iron_low,duration_dwarf_iron_up,duration_dwarf_iron_low, P_orb_dwarf_iron, P_spin_dwarf_iron)
    
def comparison_plot():
    names_micro1, peak_micro1, peak_error_micro1, total_en_micro1, total_en_error_micro1, duration_micro1, duration_error_micro1,names_donor1, peak_donor1, peak_error_donor1, total_en_donor1, total_en_error_donor1, duration_donor1, duration_error_donor1,names_gating1, peak_gating1, peak_error_gating1, total_en_gating1, total_en_error_gating1, duration_gating1, duration_error_gating1, names_dwarf1,peak_dwarf1,peak_error_dwarf1, total_en_dwarf1,total_en_error_dwarf1,duration_dwarf1,duration_error_dwarf1 = literature_data()
    names_micro, peak_micro, peak_error_micro, total_en_micro_up, total_en_micro_low, duration_micro_up, duration_micro_low,P_orb_micro, P_spin_micro,names_donor, peak_donor, peak_error_donor, total_en_donor, total_en_error_donor, duration_donor, duration_error_donor,names_gating, peak_gating, peak_error_gating, total_en_gating, total_en_error_gating, duration_gating, duration_error_gating, P_orb_gating,P_spin_gating,names_dwarf,peak_dwarf,peak_error_dwarf, total_en_dwarf_up,total_en_dwarf_low,duration_dwarf_up,duration_dwarf_low,P_orb_dwarf,P_spin_dwarf,names_micro_iron, peak_micro_iron, peak_error_micro_iron, total_en_micro_iron_up, total_en_micro_iron_low, duration_micro_iron_up, duration_micro_iron_low, P_orb_micro_iron, P_spin_micro_iron,names_dwarf_iron,peak_dwarf_iron,peak_error_dwarf_iron, total_en_dwarf_iron_up,total_en_dwarf_iron_low,duration_dwarf_iron_up,duration_dwarf_iron_low, P_orb_dwarf_iron, P_spin_dwarf_iron = new_data()
    
    midpoint_en_micro = (total_en_micro_up+total_en_micro_low)/2
    xerr_lower_micro = midpoint_en_micro - total_en_micro_low
    xerr_upper_micro = total_en_micro_up - midpoint_en_micro
    
    midpoint_en_dwarf = (total_en_dwarf_up+total_en_dwarf_low)/2
    xerr_lower_dwarf = midpoint_en_dwarf - total_en_dwarf_low
    xerr_upper_dwarf = total_en_dwarf_up - midpoint_en_dwarf
    

    plt.errorbar(midpoint_en_micro, peak_micro,xerr = [xerr_lower_micro, xerr_upper_micro],yerr = peak_error_micro,fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey', label = 'New IPs/Candidates')
    plt.errorbar(midpoint_en_dwarf, peak_dwarf,xerr = [xerr_lower_dwarf, xerr_upper_dwarf],yerr = peak_error_dwarf,fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey')
        
     
    plt.errorbar(total_en_micro1,peak_micro1,xerr = total_en_error_micro1, yerr= peak_error_micro1,label = "Micronovae",color = "blue",markersize=6, fmt = 'o',elinewidth=0.5, fillstyle = 'none')
    plt.errorbar(total_en_dwarf1, peak_dwarf1,xerr=total_en_error_dwarf1, yerr=peak_error_dwarf1,label="Dwarf Novae", color="red", markersize=6, fmt='o',elinewidth=0.5, fillstyle = 'none')
    plt.errorbar(total_en_gating1, peak_gating1,xerr=total_en_error_gating1, yerr=peak_error_gating1,label="Magnetic Gating", color="green", markersize=6, fmt='o',elinewidth=0.5, fillstyle = 'none')

    plt.xlabel("Total Optical Energy (erg)")
    plt.ylabel("Peak Optical Luminosity (erg/s)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    
    midpoint_dur_micro = (duration_micro_up+duration_micro_low)/2
    xerr_lower_micro_dur = midpoint_dur_micro - duration_micro_low
    xerr_upper_micro_dur = duration_micro_up - midpoint_dur_micro
    
    midpoint_dur_dwarf = (duration_dwarf_up+duration_dwarf_low)/2
    xerr_lower_dwarf_dur = midpoint_dur_dwarf - duration_dwarf_low
    xerr_upper_dwarf_dur = duration_dwarf_up - midpoint_dur_dwarf
    
    plt.errorbar(midpoint_en_micro, midpoint_dur_micro,xerr = [xerr_lower_micro, xerr_upper_micro],yerr = [xerr_lower_micro_dur,xerr_upper_micro_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey', label = 'New IPs/Candidates')
    plt.errorbar(midpoint_en_dwarf, midpoint_dur_dwarf,xerr = [xerr_lower_dwarf, xerr_upper_dwarf],yerr = [xerr_lower_dwarf_dur,xerr_upper_dwarf_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey')
    
    
    plt.errorbar(total_en_micro1, duration_micro1, xerr=total_en_error_micro1, yerr=duration_error_micro1, label="Micronovae", color="blue", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.errorbar(total_en_dwarf1, duration_dwarf1, xerr=total_en_error_dwarf1, yerr=duration_error_dwarf1, label="Dwarf Novae", color="red", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.errorbar(total_en_gating1, duration_gating1, xerr=total_en_error_gating1, yerr=duration_error_gating1, label="Magnetic Gating", color="green", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.xlabel("Total Optical Energy (erg)")
    plt.ylabel("Burst Duration (d)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    plt.errorbar(peak_micro, midpoint_dur_micro,xerr = peak_error_micro,yerr = [xerr_lower_micro_dur,xerr_upper_micro_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey', label = 'New IPs/ Candidates')
    plt.errorbar(peak_dwarf, midpoint_dur_dwarf,xerr = peak_error_dwarf,yerr = [xerr_lower_dwarf_dur,xerr_upper_dwarf_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey')
    
    
    plt.errorbar(peak_micro1, duration_micro1, xerr=peak_error_micro1, yerr=duration_error_micro1, label="Micronovae", color="blue", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.errorbar(peak_dwarf1, duration_dwarf1, xerr=peak_error_dwarf1, yerr=duration_error_dwarf1, label="Dwarf Novae", color="red", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.errorbar(peak_gating1, duration_gating1, xerr=peak_error_gating1, yerr=duration_error_gating1, label="Magnetic Gating", color="green", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.xlabel("Peak Luminosity (erg/s)")
    plt.ylabel("Burst Duration (d)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.show()
    
    
    plt.scatter(P_orb_dwarf,P_spin_dwarf,s=6, color = 'blue', label = "Dwarf Nova")
    plt.scatter(P_orb_micro,P_spin_micro,s=6, color = 'red', label = "Micronova")
    plt.scatter(P_orb_gating,P_spin_gating,s=6, color = 'yellow', label = "Magnetic Gating")
    
    
    P_orb_micro1 = np.array([5.49,6.43])
    P_spin_micro1 = np.array([1909.7,741.6])
    P_orb_gating1 = np.array([1.62,1.41])
    P_spin_gating1 = np.array([2280,2146])
    plt.scatter(P_orb_micro1,P_spin_micro1, s=6, color = 'red', label = "Micronova (old)", marker = '^')
    plt.scatter(P_orb_gating1,P_spin_gating1, s=6, color = 'yellow', label = "Magnetic Gating (old)", marker = '^')
    
    xvals = np.linspace(1,20, 100)

    # Lines:
    yvals_equal = xvals * 3600              # P_spin = P_orb
    yvals_tenth = 0.1 * xvals * 3600        # P_spin = 0.1 * P_orb
        
    plt.scatter(P_orb_dwarf, P_spin_dwarf, s=5, color='blue', label="Dwarf Nova")
    # ... likewise for micronova and magnetic gating ...
    plt.plot(xvals, yvals_equal, '--k',  color = 'black')
    plt.plot(xvals, yvals_tenth, '--g', color = 'black')
    
    # Pick a point in the middle of your x‐range:
    xmid = xvals[len(xvals)//2]
    ymid = yvals_equal[len(xvals)//2]
        
        # Add text right at (xmid, ymid).  No arrowprops needed.
    plt.annotate(
    r"$P_{\mathrm{spin}} = P_{\mathrm{orb}}$",   # LaTeX math‐mode label
    xy=(xmid, ymid),                            # Location in data coords
    xytext=(0, 0),                              # No offset
    textcoords='offset points',                 # Interpret xytext as an offset
    ha='center', va='bottom',                   # Centre the text horizontally
    fontsize=10, color='black'
    )

    # Similarly for P_spin = 0.1 P_orb:
    xmid_01 = xvals[len(xvals)//2]
    ymid_01 = 0.2 * xmid_01 * 3600

    plt.annotate(
    r"$P_{\mathrm{spin}} = 0.1\,P_{\mathrm{orb}}$",
    xy=(xmid_01, ymid_01),
    xytext=(0, 0),
    textcoords='offset points',
    ha='center', va='bottom',
    fontsize = 10, color = 'black')
    
    
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.xlabel("Orbital Period (hrs)")
    plt.ylabel("Spin Period (s)")
    plt.show()
    
    plt.errorbar(P_orb_dwarf, peak_dwarf, yerr = peak_error_dwarf, fmt = 'o',  markersize=5, elinewidth=0.5, capsize=2, capthick=0.5, color = 'blue', label = 'Dwarf Nova')
    plt.errorbar(P_orb_micro, peak_micro, yerr = peak_error_micro, fmt = 'o',  markersize=5, elinewidth=0.5, capsize=2, capthick=0.5, color = 'red', label = "Micronova")
    plt.xlabel("Orbital Period (d)")
    plt.ylabel("Peak Optical Luminosity (erg/s)")
    plt.legend()
    plt.yscale('log')
    
    plt.show()
    
    
    
    ratio_dwarf = P_spin_dwarf/(P_orb_dwarf*3600)
    ratio_micro = P_spin_micro/(P_orb_micro*3600)
    plt.scatter(ratio_dwarf, peak_dwarf, color = 'blue', s=5)
    plt.scatter(ratio_micro, peak_micro, color = 'red', s=5)
    plt.yscale('log')
    plt.legend()
    plt.xlabel("$P_{spin}/P_{orb}$")
    plt.ylabel("Peak Optical Luminosity (erg/s)")
    
    plt.show()
    
    
    midpoint_en_micro = (total_en_micro_iron_up+total_en_micro_iron_low)/2
    xerr_lower_micro = midpoint_en_micro - total_en_micro_iron_low
    xerr_upper_micro = total_en_micro_iron_up - midpoint_en_micro
    
    midpoint_en_dwarf = (total_en_dwarf_iron_up+total_en_dwarf_iron_low)/2
    xerr_lower_dwarf = midpoint_en_dwarf - total_en_dwarf_iron_low
    xerr_upper_dwarf = total_en_dwarf_iron_up - midpoint_en_dwarf
    
    
    plt.errorbar(midpoint_en_micro, peak_micro_iron,xerr = [xerr_lower_micro, xerr_upper_micro],yerr = peak_error_micro_iron,fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey', label = 'Ironclad IPs')
    plt.errorbar(midpoint_en_dwarf, peak_dwarf_iron,xerr = [xerr_lower_dwarf, xerr_upper_dwarf],yerr = peak_error_dwarf_iron,fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey')
     
    plt.errorbar(total_en_micro1,peak_micro1,xerr = total_en_error_micro1, yerr= peak_error_micro1,label = "Micronovae",color = "blue",markersize=6, fmt = 'o',elinewidth=0.5,fillstyle="none")
    plt.errorbar(total_en_dwarf1, peak_dwarf1,xerr=total_en_error_dwarf1, yerr=peak_error_dwarf1,label="Dwarf Novae", color="red", markersize=6, fmt='o',elinewidth=0.5,fillstyle='none')
    plt.errorbar(total_en_gating1, peak_gating1,xerr=total_en_error_gating1, yerr=peak_error_gating1,label="Magnetic Gating", color="yellow", markersize=6, fmt='o',elinewidth=0.5,fillstyle='none')

    plt.xlabel("Total Optical Energy (erg)")
    plt.ylabel("Peak Optical Luminosity (erg/s)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    plt.show()
    
    midpoint_dur_micro = (duration_micro_iron_up+duration_micro_iron_low)/2
    xerr_lower_micro_dur = midpoint_dur_micro - duration_micro_iron_low
    xerr_upper_micro_dur = duration_micro_iron_up - midpoint_dur_micro
    
    midpoint_dur_dwarf = (duration_dwarf_iron_up+duration_dwarf_iron_low)/2
    xerr_lower_dwarf_dur = midpoint_dur_dwarf - duration_dwarf_iron_low
    xerr_upper_dwarf_dur = duration_dwarf_iron_up - midpoint_dur_dwarf
    
    plt.errorbar(midpoint_en_micro, midpoint_dur_micro,xerr = [xerr_lower_micro, xerr_upper_micro],yerr = [xerr_lower_micro_dur,xerr_upper_micro_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey', label = 'Ironclad IPs')
    plt.errorbar(midpoint_en_dwarf, midpoint_dur_dwarf,xerr = [xerr_lower_dwarf, xerr_upper_dwarf],yerr = [xerr_lower_dwarf_dur,xerr_upper_dwarf_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey')
    
    
    plt.errorbar(total_en_micro1, duration_micro1, xerr=total_en_error_micro1, yerr=duration_error_micro1, label="Micronovae", color="blue", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5,fillstyle = 'none')
    plt.errorbar(total_en_dwarf1, duration_dwarf1, xerr=total_en_error_dwarf1, yerr=duration_error_dwarf1, label="Dwarf Novae", color="red", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.errorbar(total_en_gating1, duration_gating1, xerr=total_en_error_gating1, yerr=duration_error_gating1, label="Magnetic Gating", color="yellow", fmt='o', markersize=6, elinewidth=0.5, capsize=2, capthick=0.5, fillstyle = 'none')
    plt.xlabel("Total Optical Energy (erg)")
    plt.ylabel("Burst Duration (d)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    
def ecdf(data):
    x = np.sort(data)
    y = np.arange(1, len(x)+1) / len(x)
    return x, y    

    
def ritter_kolb():
    names_micro1, peak_micro1, peak_error_micro1, total_en_micro1, total_en_error_micro1, duration_micro1, duration_error_micro1,names_donor1, peak_donor1, peak_error_donor1, total_en_donor1, total_en_error_donor1, duration_donor1, duration_error_donor1,names_gating1, peak_gating1, peak_error_gating1, total_en_gating1, total_en_error_gating1, duration_gating1, duration_error_gating1, names_dwarf1,peak_dwarf1,peak_error_dwarf1, total_en_dwarf1,total_en_error_dwarf1,duration_dwarf1,duration_error_dwarf1 = literature_data()
    names_micro, peak_micro, peak_error_micro, total_en_micro_up, total_en_micro_low, duration_micro_up, duration_micro_low,P_orb_micro, P_spin_micro,names_donor, peak_donor, peak_error_donor, total_en_donor, total_en_error_donor, duration_donor, duration_error_donor,names_gating, peak_gating, peak_error_gating, total_en_gating, total_en_error_gating, duration_gating, duration_error_gating, P_orb_gating,P_spin_gating,names_dwarf,peak_dwarf,peak_error_dwarf, total_en_dwarf_up,total_en_dwarf_low,duration_dwarf_up,duration_dwarf_low,P_orb_dwarf,P_spin_dwarf,names_micro_iron, peak_micro_iron, peak_error_micro_iron, total_en_micro_iron_up, total_en_micro_iron_low, duration_micro_iron_up, duration_micro_iron_low, P_orb_micro_iron, P_spin_micro_iron,names_dwarf_iron,peak_dwarf_iron,peak_error_dwarf_iron, total_en_dwarf_iron_up,total_en_dwarf_iron_low,duration_dwarf_iron_up,duration_dwarf_iron_low, P_orb_dwarf_iron, P_spin_dwarf_iron = new_data()
    Vizier.ROW_LIMIT = -1
    
    orbital_periods_ip = np.concatenate((P_orb_micro_iron, P_orb_dwarf_iron))

    # Replace 'J/A+A/...' with the actual catalogue identifier for the Ritter–Kolb catalogue.
    # You can find this identifier on the Vizier page (it might look something like "J/A+A/516/A66" or similar).
    catalog_list = Vizier.get_catalogs("2003A&A...404..301R")  # example identifier
    
    # The catalog is usually the first (or only) table returned.
    catalog_table = catalog_list[0]
    
    # Convert the table to a Pandas DataFrame.
    df = catalog_table.to_pandas()
    
    df['Orb.Per'] = pd.to_numeric(df['Orb.Per'], errors='coerce')

    # Drop rows where Porb is NaN
    df = df.dropna(subset=['Orb.Per'])
    
    # Now extract the orbital_periods as a numpy array
    orbital_periods_df = df['Orb.Per'].values
    

    # If you only want the orbital period column, assuming its name is 'Porb' (adjust if different):
    #orbital_periods_df = df[['Orb.Per']]
    names_df = df[['Name']]
    
    print(orbital_periods_df)
    print(names_df)
    

    ips_x, ips_y = ecdf(orbital_periods_ip)
    cvs_x, cvs_y = ecdf(orbital_periods_df*24)
    
    plt.step(ips_x, ips_y, label='IP Bursts', where='post')
    plt.step(cvs_x, cvs_y, label='CVs', where='post')
    plt.xlabel('Orbital Period')
    plt.ylabel('ECDF')
    plt.legend()
    plt.xlim(0,30)
    plt.show()
    
    ks_stat, p_value = stats.ks_2samp(orbital_periods_ip, orbital_periods_df)
    print("KS statistic:", ks_stat)
    print("p-value:", p_value)
    
    
    star,porb_hrs,pspin_s,porb_pspin = confirmed_IPs()
    ips_x, ips_y = ecdf(orbital_periods_ip)
    ips_conf_x, ips_conf_y = ecdf(porb_hrs)
    
    plt.step(ips_x, ips_y, label='IP Bursts', where='post')
    plt.step(ips_conf_x, ips_conf_y, label='IPs', where='post')
    plt.xlabel('Orbital Period')
    plt.ylabel('ECDF')
    plt.legend()
    plt.xlim(0,30)
    plt.show()
    
    
    

def confirmed_IPs():
    # URL of the page containing the table
    url = 'https://asd.gsfc.nasa.gov/Koji.Mukai/iphome/catalog/omegaomega.html'
    
    # Read all tables on the page; this returns a list of DataFrames
    tables = pd.read_html(url)
    
    # Check how many tables were found
    print(f"Found {len(tables)} tables.")
    
    # Assuming the table you want is the first one:
    df = tables[2]
        
        # Optionally, inspect the DataFrame
    star = df['Star']
    porb_pspin = df['Porb/Pspin']
    porb_hrs = df['Porb (hrs)']
    pspin_s = df['Pspin (s)']
    return star,porb_pspin,pspin_s,porb_pspin
        
        
#read_files()
time_mask_V2,flux_mask_V2, flux_error_mask_V2, time_mask_g2,flux_mask_g2,flux_error_mask_g2 = specific_file()
burst_focus(time_mask_V2,flux_mask_V2, flux_error_mask_V2, time_mask_g2,flux_mask_g2,flux_error_mask_g2,3.777*10**(43), 'g')
#literature_data()
#new_data()
#comparison_plot()
#ritter_kolb()
#confirmed_IPs()
