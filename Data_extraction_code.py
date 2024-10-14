# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:05:10 2024

@author: jakem
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

def bin_data(x, y, bin_size):
    binned_x = np.mean(x.reshape(-1, bin_size), axis=1)
    binned_y = np.mean(y.reshape(-1, bin_size), axis=1)
    return binned_x, binned_y

def Lomb_Scargle(time,flux,file_name):
    # Specify the desired frequency range
    min_freq, max_freq = 1,1000
    # Compute the Lomb-Scargle Periodogram within the specified frequency range
    num_frequency_points = 1000000  # You can adjust this based on the desired resolution
    frequency = np.linspace(min_freq, max_freq, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency)  # Manually compute power without autopower
    # Plot the Lomb-Scargle Periodogram
    plt.figure(figsize=(10, 6))
    plt.plot(frequency,power, 'k', lw=1)
    # Set logarithmic scale for frequency and power if needed
    plt.xscale('log')
    plt.yscale('log')
    # Set axis labels
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Power x Frequency')
    # Set plot title
    plt.title('Lomb-Scargle Periodogram (Specified Frequency Range)')
    
    peaks, _ = find_peaks(power*frequency, height=0)  # Find all peaks
    sorted_peaks = peaks[np.argsort(power[peaks])][-3:]

    
    # Show the plot
    plt.show()
    return frequency,power  

fits_file = "C:/Users/jakem/OneDrive/Documents/Year 4 Project/TESS Raw Data/MAST_2024-10-14T1021/MAST_2024-10-14T1021/TESS/tess2023069172124-s0063-0000000393471167-0255-s/tess2023069172124-s0063-0000000393471167-0255-s_tp.fits"
print(fits.getdata(fits_file, ext=1).columns)

with fits.open(fits_file, mode="readonly") as hdulist:
    tess_bjds = hdulist[1].data['TIME']
    flux = hdulist[1].data['FLUX']
    print(flux)
    #pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
    
    flux_mean = np.mean(flux, axis=(1, 2))
    qual_flags = hdulist[1].data['QUALITY']
    
    good_data = qual_flags == 0
    tess_bjds_filtered = tess_bjds[good_data]
    flux_mean_filtered = flux_mean[good_data]
    


print(qual_flags)
#Start figure and axis.
fig, ax = plt.subplots()
where_gt0 = np.where(qual_flags > 0)[0]
# Plot the timeseries in black circles.
#ax.plot(tess_bjds, flux_mean, 'ko')
#ax.plot(tess_bjds[where_gt0], flux_mean[where_gt0], 'ro')
ax.plot(tess_bjds_filtered, flux_mean_filtered, 'ko',markersize = 1)

# Bin the data with a bin size of 10 (you can adjust this value)
binned_time, binned_flux = bin_data(tess_bjds_filtered[:len(tess_bjds_filtered)//10*10], flux_mean_filtered[:len(flux_mean_filtered)//10*10], bin_size=10)

# Start figure and axis with larger width
fig, ax = plt.subplots(figsize=(14, 6))  # Increase the width to 14 inches, height to 6 inches

# Plot the smoothed or filtered data (depending on what you've applied)
ax.plot(binned_time, binned_flux, 'k-')

# Optionally add labels
ax.set_xlabel("Time from 2022-10-29 20:00:42.091 (d)")
ax.set_ylabel("Flux (Arbitrary Units)")

# Show the plot
plt.show()

Lomb_Scargle(tess_bjds_filtered,flux_mean_filtered, fits_file)