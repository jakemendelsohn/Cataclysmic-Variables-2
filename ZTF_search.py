# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:35:53 2025

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

def name_to_coord(target_name): 
    result_table = Simbad.query_object(target_name)
    if result_table is None:
        print(f"Could not resolve object name: {target_name}")
    else:
        ra_hms = result_table["RA"][0]   # in H:M:S
        dec_dms = result_table["DEC"][0] # in D:M:S
            
        # Convert to an Astropy SkyCoord
        coord = SkyCoord(ra_hms, dec_dms, unit=(u.hour, u.deg))
        print(coord)
        
        ra_deg = coord.ra.degree
        dec_deg = coord.dec.degree
    return ra_deg, dec_deg

def ZTF_single(ra,dec):
            
    
    zquery = query.ZTFQuery()
    lcq = lightcurve.LCQuery.from_position(ra,dec, 20)
    ZTF_data = pd.DataFrame({'JD' : lcq.data.mjd+2400000.5, 'Magnitude' : lcq.data.mag, 'Magnitude_Error' : lcq.data.magerr, "Filter" : lcq.data.filtercode})
    #data = lcq.download_data()
    #lcq = lightcurve.LCQuery(data)
    df = ZTF_data
    filter_list = ["zg", "zr", "zi"]
    colour = ["teal", "red", "black"]
    for filt, colour in zip(filter_list, colour):     
        plt.scatter(df.JD[df.Filter == filt]-2457000, df.Magnitude[df.Filter == filt], color=colour, label=filt,s=1)
        plt.gca().invert_yaxis()
        plt.legend() 
        plt.xlabel("Time, BJD")
        plt.ylabel("Magnitude")
        plt.xlim(2450,2700)
        plt.title("ZTF light Curves in Different Bands")
        plt.show()
    try:
        data = lcq.download_data()
        if data is not None:
              print("Data downloaded successfully")
              print(data)
        else:
            print("No data returned.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
def find_files_with_starname(name, folder_path):
    matching_files = []
    for fname in os.listdir(folder_path):
        # Check if the star_name string is part of the filename
        if name in fname:
            matching_files.append(fname)
    return matching_files
        

def ztf_file(name,folder_path):
    
    matching_files = find_files_with_starname(name,folder_path)
    print(matching_files)
    
    for fls in matching_files:
        filepath = folder_path+fls
        df = pd.read_csv(filepath, delimiter='\t')
        filename = os.path.basename(filepath)
        # Split into name and extension: e.g. "somefile_g" and ".csv"
        name, _ = os.path.splitext(filename)
        prefix, suffix = name.split("_ZTF_LC_")
        band = suffix
        name = prefix
        
        time = df.iloc[:, 1]-2458000
        flux = df.iloc[:, 2]
        flux_error = df.iloc[:, 3]
        plt.scatter(time,flux,s=5, label = band)
        plt.ylabel("Flux (mJy)")
        plt.xlabel("Time (JD-2458000)")
        plt.legend(title = "band")
        plt.title(name)
    plt.show()
        
def IP_candidates():
    filename = "C:/Users/jakem/OneDrive/Documents/Year 4 Project/IPs_all_candidates_list.dat"

    star_names = []
    
    with open(filename, "r") as f:
        header = next(f)
        
        for line in f:
            parts = line.strip().split()
            
            name = parts[0] + " " + parts[1]
            star_names.append(name)

    star_names_array = np.array(star_names)
    return star_names_array

def all_candidates():
    star_names_array = IP_candidates()
    for name in star_names_array:
        name_underscored = name.replace(" ", "_")
        print(name_underscored)
        ztf_file(name_underscored,"C:/Users/jakem/OneDrive/Documents/Year 4 Project/Compiled ZTF data/ZTF_IPs_data/")
    
                
                
                
                
                
    
    
    
    
    
    

    
        
#ra,dec =name_to_coord("GI Mon")
#ZTF_single(ra,dec)
#ztf_file(name, "C:/Users/jakem/OneDrive/Documents/Year 4 Project/Compiled ZTF data/ZTF_IPs_data/")
all_candidates()