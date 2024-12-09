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
            
            plt.scatter(time_mask_V2,flux_mask_V2,s=1,label = "V band")
            plt.scatter(time_mask_g2, flux_mask_g2, s=1, label ="g band")
            plt.legend(title = f"{file}")
            plt.show()
       
    
read_files()