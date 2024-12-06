# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 16:28:17 2024

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
from ztfquery import lightcurve
from astropy import time
from ztfquery import query
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.interpolate import interp1d

def positive_nsh(P_psh):
    model = 0.003+0.481*P_psh
    error_pos = 0.006+0.501*P_psh
    error_neg = 0+0.461*P_psh
    
    # Plot the model
    plt.figure(figsize=(8, 6))
    plt.plot(P_psh, model, label="Model", color="blue")
    
    # Plot the shaded error region
    plt.fill_between(P_psh, error_neg, error_pos, color="gray", alpha=0.5, label="Error range")
    
    # Add labels and legend
    plt.xlabel("$P_{psh}$")
    plt.ylabel("$\epsilon$")
    plt.title("Model with Error Bounds")
    plt.legend()
    plt.show()
    return model,error_pos,error_neg

def positive_orb(P_orb):
    model = 0 + 0.540*P_orb
    error_pos = 0.0004+0.54026*P_orb
    error_neg = -0.0004+0.53974*P_orb
    return model,error_pos,error_neg

P_psh_vals = np.linspace(1.2,8.4,100)
P_orb_vals = np.linspace(1.2,8.4,100)

positive_nsh(P_psh_vals)