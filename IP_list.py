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
            time = df.iloc[:, 0]
            flux = df.iloc[:, 7]
            band = df.iloc[:, 9]
            indices_V = [i for i in range(len(band)) if band[i] == 'V']
            indices_g = [h for h in range(len(band)) if band[h] == 'g']
            
            time_mask_V = []
            flux_mask_V = []
            flux_mask_V2 = []
            time_mask_V2 = []
            time_mask_g = []
            flux_mask_g = []
            flux_mask_g2 = []
            time_mask_g2 = []
            
            
            for j in range(0,len(indices_V)):
                time_mask_V.append(time[indices_V[j]])
                flux_mask_V.append(flux[indices_V[j]])
                
            indices_cutoff_V = [f for f in range(len(flux_mask_V)) if flux_mask_V[f]<95]
            
            for i in range(0,len(indices_cutoff_V)):
                time_mask_V2.append(time_mask_V[indices_cutoff_V[i]])
                flux_mask_V2.append(flux_mask_V[indices_cutoff_V[i]])
            
            for k in range(0,len(indices_g)):
                time_mask_g.append(time[indices_g[k]])
                flux_mask_g.append(flux[indices_g[k]])
                
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
                plt.xlabel('Time')
                plt.ylabel('Flux')
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
       
    
read_files()