# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 16:55:24 2024

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
from matplotlib.animation import FuncAnimation


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

def ATLAS_data(time):
    file_path = 'C:/Users/jakem/OneDrive/Documents/Year 4 Project/ATLAS Data/ATLAS J5016.txt'
    data = pd.read_csv(file_path, delim_whitespace=True, comment='#', header=None)

    # Assign column names to the data if not provided in the file
    data.columns = ["mJD", "m", "dm", "uJy", "duJy", "F", "f", "chi/N", "RA", "Dec", "x", "y", "maj", "min", "phi", "apfit", "mag5sig", "Sky", "Obs"]
    filtered_data = data[data["uJy"] > 0]

    # Extract mJD and uJy columns as separate variables after filtering
    mjd = filtered_data["mJD"]
    flux = filtered_data["uJy"]
    
    # Define the target coordinates (RA and Dec in degrees)
    ra_deg = 303.9042  # Replace with your target RA in degrees
    dec_deg = 37.1897  # Replace with your target Dec in degrees
    target_coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    
    # Define the observatory location (example: Mauna Kea)
    # Replace with your observation location if different
    observatory_location = EarthLocation.of_site('ALMA')  # or specify lat/lon/height
    
    mjd_times = Time(filtered_data["mJD"].values, format='mjd', scale='utc', location=observatory_location)
    bjd_times = mjd_times.tdb + mjd_times.light_travel_time(target_coord)

    # Verify the conversion by comparing MJD and BJD_TDB values
    filtered_data["BJD_TDB"] = bjd_times.value
    plt.scatter(filtered_data["BJD_TDB"], filtered_data["uJy"])
    plt.xlabel("BJD_TDB")
    plt.ylabel("uJy")
    plt.title("Flux (uJy) vs. BJD_TDB")
    plt.show()
    
    filtered_data["BJD_TDB"] = bjd_times.value
    bjd = filtered_data["BJD_TDB"]
    
    plt.scatter(filtered_data["BJD_TDB"], filtered_data["uJy"])
    plt.xlabel("BJD_TDB")
plt.ylabel("uJy")
plt.title("Flux (uJy) vs. BJD_TDB")
plt.show()
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
    
    tess_flux = flux
    assasn_flux = g_band_flux_sector
    
    ground_time_1st_half = g_band_BJD_sector[:len(g_band_BJD_sector)//2]
    ground_flux_1st_half = assasn_flux[:len(assasn_flux)//2]
    
    plt.scatter(ground_time_1st_half, ground_flux_1st_half)
    plt.show()
    
    closest_indices = [np.abs(time - f).argmin() for f in ground_time_1st_half]
    matched_flux_values = tess_flux[closest_indices]
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
    plt.legend()
    plt.show()

    ground_time_2nd_half = g_band_BJD_sector[len(g_band_BJD_sector)//2:]
    ground_flux_2nd_half = assasn_flux[len(assasn_flux)//2:]
    closest_indices2 = [np.abs(time - f).argmin() for f in ground_time_2nd_half]
    matched_flux_values2 = tess_flux[closest_indices2]
    plt.scatter(matched_flux_values2,ground_flux_2nd_half)
    plt.xlabel("TESS Flux")
    plt.ylabel("Ground Based Flux")
    # Perform a linear least-squares fit
    slope2, intercept2, r_value2, p_value, std_err = linregress(matched_flux_values2, ground_flux_2nd_half)
    # Plot the best-fit line
    x_fit2 = np.array([min(matched_flux_values2), max(matched_flux_values2)])
    y_fit2 = slope2 * x_fit2 + intercept2
    plt.plot(x_fit2, y_fit2, color="red", label=f"Fit: y = {slope2:.2f}x + {intercept2:.2f}")
    plt.legend()
    plt.show()
    

    # Print the fit parameters
    print("Slope1:", slope)
    print("Intercept1:", intercept)
    print("R-squared1:", r_value**2)
    
    print("Slope2:", slope2)
    print("Intercept2:", intercept2)
    print("R-squared2:", r_value2**2)
    
    return slope,intercept, r_value, slope2, intercept2,r_value2

def calibrated_data(slope,intercept,slope2,intercept2,time,flux,g_band_BJD_sector,g_band_flux_sector):
    flux_calibrated1 = []
    flux_calibrated2 = []
    mid_index = len(flux) // 2

    # Split the list into two halves
    first_half = flux[:mid_index]
    second_half = flux[mid_index:]
    for fl in first_half:
        fl = (fl*slope)+intercept
        flux_calibrated1.append(fl)
    for f2 in second_half:
        f2 = (f2*slope2)+intercept2
        flux_calibrated2.append(f2)
    
    flux_total = flux_calibrated1+flux_calibrated2
    plt.scatter(time,flux_total,s=1)
    plt.scatter(g_band_BJD_sector, g_band_flux_sector,color = 'black', s=5)
    
def find_combinations(fundamentals, freqs, max_n=5, max_m=5):
    matched_combinations = []
    for i, f1 in enumerate(fundamentals):
        for j, f2 in enumerate(fundamentals):
            if i >= j:  # Avoid duplicate combinations (since f1 and f2 are interchangeable)
                continue
            for n in range(-max_n, max_n + 1):
                for m in range(-max_m, max_m + 1):
                    if n == 0 and m == 0:
                        continue  # skip the trivial case of n=m=0
                    combination_freq = n * f1 + m * f2
                    if combination_freq > 0 and any(np.isclose(combination_freq, f, atol=tolerance) for f in freqs):
                        matched_combinations.append((n, f1, m, f2, combination_freq))
    return matched_combinations

#time,flux,exptime = sector_data(2)
#ATLAS_data(time)
#g_band_BJD_sector, g_band_flux_sector = ASSASSN_data(time)
#slope,intercept,r_value,slope2,intercept2,r_value2 = calibration(time,flux, g_band_BJD_sector, g_band_flux_sector)
#calibrated_data(slope,intercept,slope2, intercept2,time,flux,g_band_BJD_sector,g_band_flux_sector)
all_frequencies = [0.27,0.88,1.788,2.178,6.68,7.633,8.18,8.45,15.266,15.74,16.36,24.54,32.72,25.37,25.85]  # replace with your list of frequencies
fundamental_frequencies = [0.27,0.88,7.63,8.18]
tolerance = 1e-2 # tolerance to account for measurement errors in matching harmonics
# Find matched combinations in the original list
matched_frequencies = find_combinations(fundamental_frequencies, all_frequencies)
print("Matched combinations (n, f1, m, f2, combination frequency):")
for match in matched_frequencies:
    print(match)