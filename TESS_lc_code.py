# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:17:32 2024

@author: jakem
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import lightkurve as lk


def sector_data():
    search_result = lk.search_lightcurve('EX Hydrae', mission='TESS')
    print(search_result)
    lc = search_result[3].download()
    sap_lc = lc.SAP_FLUX
    sap_lc_cleaned = sap_lc.remove_nans()
    sap_lc_cleaned.plot()
    plt.figure(figsize=(12, 6))
    plt.show()
    time = sap_lc.time.value
    flux = sap_lc.flux.value
    return time,flux
    

def frequency_range(time,flux,file_name):
    N = len(time)
    if file_name.endswith('LC.csv'):
        del_t = 120
    elif file_name.endswith('SC.csv'):
        del_t = 20
    else:
        raise ValueError("The file name does not end with 'LC' or 'SC'.")
    seconds_per_day = 86400
    min_freq = 3/(N*del_t)*seconds_per_day
    max_freq = 1/(2*del_t)*seconds_per_day
    return min_freq,max_freq
