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
from sklearn.decomposition import PCA

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
    #plt.scatter(delta_sorted, cdf,s=1)
    #plt.xlabel('Delta time')
    #plt.ylabel('CDF')
    #plt.show()
    
    return delta_sorted,cdf

def LC_model(name):
    time,flux,exptime,flux_error = ASAS_SN_19bh(1,name)
    time_mask = (time >= 2347) & (time <= 2357)
    time = time[time_mask]
    flux = flux[time_mask]
    
    factor = (3.5*10**(34))/np.max(flux)
    flux = flux*factor
    
    
    
    delta_sorted, cdf = ASAS_SN_delt()
    
    plt.scatter(time,flux,s=1,color = "black", label = "Original Data")
    
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
    
    plt.scatter(time_array,flux_array,color = "red", s=5,label = "Simulated Points")
    plt.xlabel("Time (BJD-2475000)")
    plt.ylabel("Luminosity (erg/s)")
    plt.legend()
    plt.show()
    peak_lumi = np.max(flux_array)
    
    
    avg_flux_TESS = np.mean(flux)
    print("average",avg_flux_TESS)
    threshold = 1.1*avg_flux_TESS
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
        extracted_set = set(lum_burst)
        remainder = [x for x in flux_array if x not in extracted_set]
        avg = np.mean(remainder)
    
        integral_up = np.trapz(lum_burst-avg,time_burst)
        print(f"Energy Upper Limit = {integral_up:.2e} erg")
        
        start_lower = np.where(time_array == t_start_lower)[0]
        end_lower = np.where(time_array == t_end_lower)[0]
        start_lower = start_lower[0]
        end_lower = end_lower[0]
        t_burst = time_array[start_lower:end_lower]
        lum_burst = flux_array[start_lower:end_lower]
        spd = 86400.0
        time_burst = t_burst * spd
        
        integral_low = np.trapz(lum_burst-avg,time_burst)
        print(f"Energy Lower Limit = {integral_low:.2e} erg")
        
    return peak_lumi, integral_low,integral_up,outburst_duration,duration_lower
        
        
def simulation(iterations,name):
    lum_array = np.empty(0)
    int_up_array = np.empty(0)
    int_low_array = np.empty(0)
    dur_up_array = np.empty(0)
    dur_low_array = np.empty(0)
    for i in range(0,iterations):       
        peak_lumi,integral_low, integral_up,duration_upper,duration_lower = LC_model(name)
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
    
    #plt.scatter(avg_int_per_run,lum_array)
    #plt.xlabel("Average Eenergy (ergs)")
    #plt.ylabel("Peak Luminosity (erg/s)")
    
    return avg_int_per_run,lum_array, int_low_array
    
def simulation_comparison(iterations,name):
    names_micro1, peak_micro1, peak_error_micro1, total_en_micro1, total_en_error_micro1, duration_micro1, duration_error_micro1,names_donor1, peak_donor1, peak_error_donor1, total_en_donor1, total_en_error_donor1, duration_donor1, duration_error_donor1,names_gating1, peak_gating1, peak_error_gating1, total_en_gating1, total_en_error_gating1, duration_gating1, duration_error_gating1, names_dwarf1,peak_dwarf1,peak_error_dwarf1, total_en_dwarf1,total_en_error_dwarf1,duration_dwarf1,duration_error_dwarf1 = literature_data()
    names_micro, peak_micro, peak_error_micro, total_en_micro_up, total_en_micro_low, duration_micro_up, duration_micro_low,P_orb_micro, P_spin_micro,names_donor, peak_donor, peak_error_donor, total_en_donor, total_en_error_donor, duration_donor, duration_error_donor,names_gating, peak_gating, peak_error_gating, total_en_gating, total_en_error_gating, duration_gating, duration_error_gating, P_orb_gating,P_spin_gating,names_dwarf,peak_dwarf,peak_error_dwarf, total_en_dwarf_up,total_en_dwarf_low,duration_dwarf_up,duration_dwarf_low,P_orb_dwarf,P_spin_dwarf,names_micro_iron, peak_micro_iron, peak_error_micro_iron, total_en_micro_iron_up, total_en_micro_iron_low, duration_micro_iron_up, duration_micro_iron_low, P_orb_micro_iron, P_spin_micro_iron,names_dwarf_iron,peak_dwarf_iron,peak_error_dwarf_iron, total_en_dwarf_iron_up,total_en_dwarf_iron_low,duration_dwarf_iron_up,duration_dwarf_iron_low, P_orb_dwarf_iron, P_spin_dwarf_iron = new_data()
    
    energy, lum_array, low_energies = simulation(iterations,name)
    
    midpoint_en_micro = (total_en_micro_up+total_en_micro_low)/2
    xerr_lower_micro = midpoint_en_micro - total_en_micro_low
    xerr_upper_micro = total_en_micro_up - midpoint_en_micro
    
    midpoint_en_dwarf = (total_en_dwarf_up+total_en_dwarf_low)/2
    xerr_lower_dwarf = midpoint_en_dwarf - total_en_dwarf_low
    xerr_upper_dwarf = total_en_dwarf_up - midpoint_en_dwarf
    
    
    midpoint_en_micro_iron = (total_en_micro_iron_up+total_en_micro_iron_low)/2
    xerr_lower_micro_iron = midpoint_en_micro_iron - total_en_micro_iron_low
    xerr_upper_micro_iron = total_en_micro_iron_up - midpoint_en_micro_iron
    
    midpoint_en_dwarf_iron = (total_en_dwarf_iron_up+total_en_dwarf_iron_low)/2
    xerr_lower_dwarf_iron = midpoint_en_dwarf_iron - total_en_dwarf_iron_low
    xerr_upper_dwarf_iron = total_en_dwarf_iron_up - midpoint_en_dwarf_iron
    
    midpoint_dur_micro_iron = (duration_micro_iron_up+duration_micro_iron_low)/2
    xerr_lower_micro_dur = midpoint_dur_micro_iron - duration_micro_iron_low
    xerr_upper_micro_dur = duration_micro_iron_up - midpoint_dur_micro_iron
    
    midpoint_dur_dwarf_iron = (duration_dwarf_iron_up+duration_dwarf_iron_low)/2
    xerr_lower_dwarf_dur = midpoint_dur_dwarf_iron - duration_dwarf_iron_low
    xerr_upper_dwarf_dur = duration_dwarf_iron_up - midpoint_dur_dwarf_iron
    
    central_x = 1.2*10**(39)
    central_y = 3.5*10**(34)
    
    (m1x, m1y), dir1, scale1 = principal_direction(energy, lum_array)
    dx1 = dir1[0] * scale1
    dy1 = dir1[1] * scale1
    #plt.arrow(m1x, m1y, dx1, dy1,
    #color='red', width=0,  # or a small number if you want thickness
    #head_width=0.07*m1y, length_includes_head=True,
    #label='Dist 1 direction'
    #)
    #plt.arrow(m1x, m1y, dx1, dy1,
    #      width=0.02*10**(np.floor(np.log10(np.mean(lum_array)))),
    #      head_width=0.05*abs(dy1),
    #      color='red', alpha=0.8, label='Principal direction')
    
    plt.scatter(energy,lum_array, color = "orange", alpha = 0.3, label = "Simulations")
    plt.scatter(central_x,central_y,color = "red",marker = "x", label = "ASAS-SN 19bh")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Total Optical Energy (erg)")
    plt.ylabel("Peak Optical Luminosity (erg/s)")
    plt.legend()
    plt.show()
    
    
    plt.errorbar(midpoint_en_micro_iron, peak_micro_iron,xerr = [xerr_lower_micro_iron, xerr_upper_micro_iron],yerr = peak_error_micro_iron,fmt = '^', markersize = 6, elinewidth = 0.5,color = 'blue', label = 'Micronova (new)')
    plt.errorbar(midpoint_en_dwarf_iron, peak_dwarf_iron,xerr = [xerr_lower_dwarf_iron, xerr_upper_dwarf_iron],yerr = peak_error_dwarf_iron,fmt = '^', markersize = 6, elinewidth = 0.5,color = 'red', label = "Dwarf nova (new)")
    
    plt.errorbar(total_en_micro1,peak_micro1,xerr = total_en_error_micro1, yerr= peak_error_micro1,label = "Micronovae",color = "blue",markersize=6, fmt = 'o',elinewidth=0.5, fillstyle = 'none')
    plt.errorbar(total_en_dwarf1, peak_dwarf1,xerr=total_en_error_dwarf1, yerr=peak_error_dwarf1,label="Dwarf Novae", color="red", markersize=6, fmt='o',elinewidth=0.5, fillstyle = 'none')
    plt.errorbar(total_en_gating1, peak_gating1,xerr=total_en_error_gating1, yerr=peak_error_gating1,label="Magnetic Gating", color="green", markersize=6, fmt='o',elinewidth=0.5, fillstyle = 'none')


    nonzero_energy = np.count_nonzero(energy)
    print("From", len(energy), "runs, ASAS-SN will detect on average", nonzero_energy, "bursts")
    
    print(np.count_nonzero(lum_array))
    
    plt.scatter(energy,lum_array, color = "orange", alpha = 0.3, label = "Simulations")
    plt.scatter(central_x,central_y,color = "red",marker = "x", label = "ASAS-SN 19bh")
    
    plt.xlabel("Total Optical Energy (erg)")
    plt.ylabel("Peak Optical Luminosity (erg/s)")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()
    
    
    plt.errorbar(midpoint_en_micro_iron, midpoint_dur_micro_iron,xerr = [xerr_lower_micro_iron, xerr_upper_micro_iron],yerr = [xerr_lower_micro_dur,xerr_upper_micro_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey', label = 'Ironclad IPs')
    plt.errorbar(midpoint_en_dwarf_iron, midpoint_dur_dwarf_iron,xerr = [xerr_lower_dwarf_iron, xerr_upper_dwarf_iron],yerr = [xerr_lower_dwarf_dur,xerr_upper_dwarf_dur],fmt = 'o', markersize = 6, elinewidth = 0.5,color = 'grey')
    plt.xscale("log")
    plt.show()
    return m1x, m1y, dir1, scale1
    
def principal_direction(x, y):
    X = np.column_stack((x, y))
    pca = PCA(n_components=2)
    pca.fit(X)

    # The first component is the direction of largest variance
    dir_vector = pca.components_[0]  # shape (2,)
    
    # Mean of the distribution
    mean_x, mean_y = X.mean(axis=0)
    
    # A simple scale: fraction of the data's spread
    # (Tune this to make the arrow look good on your plot)
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    scale_factor = 0.3 * np.sqrt(x_range**2 + y_range**2)
    
    return (mean_x, mean_y), dir_vector, scale_factor
    
#LC_model()
#ASAS_SN_delt()
#simulation(100)
simulation_comparison(200, 'ASASSN -19bh')