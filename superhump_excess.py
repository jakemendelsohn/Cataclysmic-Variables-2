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

def positive_psh(P_psh):
    model = 0.003+(0.481*P_psh)
    error_pos = 0.006+(0.501*P_psh)
    error_neg = 0+(0.461*P_psh)
    
    # Plot the model
    plt.figure(figsize=(8, 6))
    #plt.plot(P_psh, model, label="Model", color="blue")
    
    # Plot the shaded error region
    #plt.fill_between(P_psh, error_neg, error_pos, color="gray", alpha=0.5, label="Error range")
    
    # Add labels and legend
    plt.xlabel("$P_{psh}$")
    plt.ylabel("$\epsilon$")
    plt.title("Model with Error Bounds")
    plt.legend()

    return model,error_pos,error_neg

def positive_orb(P_orb):
    model = 0 + (0.540*P_orb)
    error_pos = 0.004+0.566*P_orb
    error_neg = -0.004+0.514*P_orb
    return model,error_pos,error_neg

def negative_nsh(P_nsh):
    model = -0.014 - (0.174*P_nsh)
    error_pos = model + (0.005 + 0.006 * P_nsh)  # Shifted above the model
    error_neg = model - (0.005 + 0.006 * P_nsh)
    return model,error_pos,error_neg

def negative_orb(P_orb):
    model = -0.014-(0.165*P_orb)
    error_pos = model + (0.005 + 0.006 * P_orb)  # Shifted above the model
    error_neg = model - (0.005 + 0.006* P_orb)
    return model,error_pos,error_neg

def combined():
    P_psh = np.linspace(0.05, 0.35, 100)  # x-axis range for P_psh
    P_orb = np.linspace(0.05, 0.35, 100)  # x-axis range for P_orb
    P_nsh = np.linspace(0.05, 0.35, 100)
    
    P_orb_actual = 0.1222
    P_psh_actual = 0.132
    P_nsh_actual = 0.118
    e_neg_actual = (P_nsh_actual-P_orb_actual)/P_orb_actual
    e_pos_actual = (P_psh_actual-P_orb_actual)/P_orb_actual
    
    fig, axs = plt.subplots(2, 2, figsize=(14, 14))
    fig.subplots_adjust(hspace=0.17, wspace=0.17)
    
    # Positive psh
    ax1 = axs[0, 0]
    model, error_pos, error_neg = positive_psh(P_psh)
    ax1.plot(P_psh, model, color='magenta', label=r"$\epsilon^{+} - psh$")
    ax1.fill_between(P_psh, error_neg, error_pos, color="lightblue", alpha=0.5, label="Uncertainties")
    ax1.set_title("(a) $\epsilon^{+} - psh$")
    ax1.set_xlabel("psh [d]")
    ax1.set_ylabel(r"$\epsilon^{+}$")
    ax1.scatter(P_psh_actual,e_pos_actual,label = "MGAB V-247")
    ax1.legend()

    
    # Positive orb
    ax2 = axs[0, 1]
    model, error_pos, error_neg = positive_orb(P_orb)
    ax2.plot(P_orb, model, color='magenta', label=r"$\epsilon^{+} - orb$")
    ax2.fill_between(P_orb, error_neg, error_pos, color="peachpuff", alpha=0.5, label="Uncertainties")
    ax2.set_title("(b) $\epsilon^{+} - orb$")
    ax2.set_xlabel("orb [d]")
    ax2.set_ylabel(r"$\epsilon^{+}$")
    ax2.scatter(P_orb_actual,e_pos_actual)
    ax2.legend()

    
    # Negative psh
    ax3 = axs[1, 0]
    model, error_pos, error_neg = negative_nsh(P_nsh)
    ax3.plot(P_psh, model, color='magenta', label=r"$\epsilon^{-} - psh$")
    ax3.fill_between(P_psh, error_neg, error_pos, color="lightblue", alpha=0.5, label="Uncertainties")
    ax3.set_title("(c) $\epsilon^{-} - psh$")
    ax3.set_xlabel("psh [d]")
    ax3.set_ylabel(r"$\epsilon^{-}$")
    ax3.scatter(P_nsh_actual,e_neg_actual)
    ax3.legend()
    
    # Negative orb
    ax4 = axs[1, 1]
    model, error_pos, error_neg = negative_orb(P_orb)
    ax4.plot(P_orb, model, color='magenta', label=r"$\epsilon^{-} - orb$")
    ax4.fill_between(P_orb, error_neg, error_pos, color="lightblue", alpha=0.5, label="Uncertainties")
    ax4.set_title("(d) $\epsilon^{-} - orb$")
    ax4.set_xlabel("orb [d]")
    ax4.set_ylabel(r"$\epsilon^{-}$")
    ax4.scatter(P_orb_actual,e_neg_actual)
    ax4.legend()

    
    plt.show()
    
    
    
P_psh_vals = np.linspace(1.2,8.4,100)
P_orb_vals = np.linspace(1.2,8.4,100)

#positive_psh(P_psh_vals)
combined()