# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 11:18:27 2024

@author: jakem
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import lightkurve as lk
from scipy.optimize import curve_fit
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from scipy.stats import linregress


def time_conversion():
    hjd_time = Time('2460395.08940', format='jd', scale='utc',location=EarthLocation(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m))  # example HJD

    # Define the target's sky coordinates (RA and Dec) and observer's location
    ra = '20h15m36.96s'  # Right Ascension of target
    dec = '37d11m23s'  # Declination of target
    star_coords = SkyCoord(ra, dec, frame='icrs')
    
    # Convert HJD to BJD
    bjd_time = hjd_time + hjd_time.light_travel_time(star_coords, kind = 'barycentric')
    
    print("HJD:", hjd_time.jd)
    print("BJD:", bjd_time.jd)
    return bjd_time

def sector_data(index):
    search_result = lk.search_lightcurve('RX J2015.6+3711', mission='TESS')
    print(search_result)
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    #sap_lc_cleaned = sap_lc.remove_nans()
    quality_flags = sap_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    print(f"Data points with quality flags at indices: {flagged_indices}")
    print(len(flagged_indices))
    good_quality_mask = quality_flags == 0  # Keeps only points with a quality flag of 0 (good data)
    
    time = sap_lc.time.value[good_quality_mask]
    flux = sap_lc.flux.value[good_quality_mask]
    #sap_lc_cleaned = sap_lc.remove_outliers()
    #time = sap_lc_cleaned.time.value
    #flux = sap_lc_cleaned.flux.value
    return time,flux,exptime

def ASSASSN_data(time):
    file_path = 'C:/Users/jakem/OneDrive/Documents/Year 4 Project/ASSASSN Data/light_curve_b390fa7d-c639-466f-ad42-1d86a873bed5_BJD.xlsx'  # Replace with the path to your Excel file
    data = pd.read_excel(file_path)

    # Extract the HJD, flux, and flux_err columns into separate variables
    BJD = data['BJD']
    flux = data['flux(mJy)']
    flux_err = data['flux_err']
    
    
    g_band_data = data[(data['Filter'] == 'g') & (data['flux(mJy)'] < 40) & (data['flux(mJy)'] > 1.8)]

    # Save HJD and flux columns in separate variables
    g_band_BJD = g_band_data['BJD']
    g_band_BJD_corrected = []
    for val in g_band_BJD:
        val = val - 2457000
        g_band_BJD_corrected.append(val)
    
    g_band_flux = g_band_data['flux(mJy)']
        
    plt.scatter(g_band_BJD_corrected,g_band_flux,s=1)
    plt.xlabel("Time (BJD-2457000)")
    plt.ylabel("Flux")
    plt.show()
    
    
    closest_start, closest_end = select_section(g_band_BJD_corrected,time[0],time[-1])
    
    g_band_BJD_sector = g_band_BJD_corrected[closest_start:closest_end]
    g_band_flux_sector = g_band_flux[closest_start:closest_end]
        
    
    # Display or use the variables as needed
    print("G Band Time Data (HJD):", g_band_BJD)
    print("G Band Flux Data (flux(mJy)):", g_band_flux)
    
    plt.scatter(g_band_BJD_sector,g_band_flux_sector,s=1)
    plt.xlabel("BJD-2457000")
    plt.ylabel("Flux")
    plt.show()
    return g_band_BJD_sector, g_band_flux_sector
    

def select_section(data, start, end):
    closest_index_start = min(range(len(data)), key=lambda i: abs(data[i] - start))
    closest_index_end = min(range(len(data)), key=lambda i: abs(data[i] - end))
    print("The closest start index is", closest_index_start)
    print("The closest end index is", closest_index_end)
    return closest_index_start,closest_index_end

def calibration(time,flux,g_band_BJD_sector, g_band_flux_sector):
    plt.scatter(time,flux,s=1)
    plt.xlabel("BJD-2457000")
    plt.ylabel("Flux (e/s)")
    plt.show()
    
    ground_time_1st_half = g_band_BJD_sector[:len(g_band_BJD_sector)//2]
    ground_flux_1st_half = g_band_flux_sector[:len(g_band_flux_sector)//2]
    
    closest_time_values = np.array([time[np.abs(time - f).argmin()] for f in ground_time_1st_half])
    closest_indices = [np.abs(time - f).argmin() for f in ground_time_1st_half]
    matched_flux_values = flux[closest_indices]
    print("ground",ground_flux_1st_half)
    print(matched_flux_values)
   # plt.scatter(ground_time_1st_half,ground_flux_1st_half,s=1)
    plt.scatter(matched_flux_values,ground_flux_1st_half)
    plt.xlabel("TESS Flux")
    plt.ylabel("Ground Based Flux")

    # Perform a linear least-squares fit
    slope, intercept, r_value, p_value, std_err = linregress(matched_flux_values, ground_flux_1st_half)

    # Plot the best-fit line
    x_fit = np.array([min(matched_flux_values), max(matched_flux_values)])
    y_fit = slope * x_fit + intercept
    plt.plot(x_fit, y_fit, color="red", label=f"Fit: y = {slope:.2f}x + {intercept:.2f}")

    # Show fit details
    plt.legend()
    plt.show()

    # Print the fit parameters
    print("Slope:", slope)
    print("Intercept:", intercept)
    print("R-squared:", r_value**2)
    
    return slope,intercept, r_value

def calibrated_data(slope,intercept,time,flux):
    flux_calibrated = []
    for fl in flux:
        fl = (fl*slope)+intercept
        flux_calibrated.append(fl)
    plt.scatter(time,flux,s=1)
        
    
    

time,flux,exptime = sector_data(2)
g_band_BJD_sector, g_band_flux_sector = ASSASSN_data(time)
slope,intercept,r_value = calibration(time,flux, g_band_BJD_sector, g_band_flux_sector)
calibrated_data(slope,intercept,time,flux)