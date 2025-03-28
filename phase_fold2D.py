# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:55:54 2024

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
from scipy.ndimage import gaussian_filter

def sector_data(index):
    search_result = lk.search_lightcurve('04 07 4.080 +18 55 37.20', mission='TESS')
    print(search_result)
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    #sap_lc_cleaned = sap_lc.remove_nans()
    quality_flags = sap_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0  # Keeps only points with a quality flag of 0 (good data)
    
    time = sap_lc.time.value[good_quality_mask]
    flux = sap_lc.flux.value[good_quality_mask]
    flux_error = sap_lc.flux_err.value[good_quality_mask]  
    #sap_lc_cleaned = sap_lc.remove_outliers()
    #time = sap_lc_cleaned.time.value
    #flux = sap_lc_cleaned.flux.value
    return time,flux,exptime,flux_error

def bin_folded_data(phase, flux, num_bins=50):
    bins = np.linspace(0, 1, num_bins+1)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    bin_means = np.zeros(num_bins)
    bin_errors = np.zeros(num_bins)
    
    for i in range(num_bins):
        in_bin = (phase >= bins[i]) & (phase < bins[i+1])
        bin_means[i] = np.mean(flux[in_bin])
        bin_errors[i] = np.std(flux[in_bin]) / np.sqrt(np.sum(in_bin))  # Standard error
    
    return bin_centers, bin_means, bin_errors

def phase_fold_binned(time, flux, peak_frequencies):
    peak_periods = 1/peak_frequencies 
    colors = ['r','g','b','y','c']
    
    for i in range(0,len(peak_periods)):
        phase = (time % peak_periods[i]) / peak_periods[i]
        bin_centers, bin_means, bin_errors = bin_folded_data(phase,flux)
        plt.errorbar(np.concatenate([bin_centers, bin_centers +1]), 
                     np.concatenate([bin_means, bin_means]), 
                     yerr=np.concatenate([bin_errors, bin_errors]), 
                     fmt='k.', label=peak_frequencies[i],color = colors[i])
    plt.legend(title = "Peak Frequencies (c/d)")
    plt.xlabel("Phase")
    plt.ylabel("Flux e/s")
    return bin_centers,bin_means
    
def phase_fold2D(time, flux, freq1, freq2, num_bins=30,smooth=True):
    freq1 = np.array([freq1])
    freq2 = np.array([freq2])
    
    period1 = 1/freq1
    period2 = 1/freq2
    

    # Calculate the phases for both frequencies
    phase1 = ((time%period1)/period1)
    phase2 = ((time%period2)/period2)
    


    phase1_binned = np.linspace(0,1,num_bins+1)
    phase2_binned = np.linspace(0,1,num_bins+1)
    
    print("Phase1 range:", phase1.min(), phase1.max())
    print("Phase2 range:", phase2.min(), phase2.max())
    M = np.zeros([num_bins,num_bins])

    for i in range(0, len(phase1_binned)-1):
        for j in range(0,len(phase2_binned)-1):
            idx = np.where(
            (phase1 > phase1_binned[i]) & (phase1 < phase1_binned[i + 1]) &
            (phase2 > phase2_binned[j]) & (phase2 < phase2_binned[j + 1])
        )[0]
            #print(idx)
            M[i,j] = np.mean(flux[idx])
    print("M",M)
    if smooth:
      # Fill NaNs so the filter doesn't propagate them
      M_filled = np.copy(M)
      nan_mask = np.isnan(M_filled)
      # Replace NaNs with the median (or mean) of non-NaN values
      M_filled[nan_mask] = np.nanmedian(M_filled)
      
      # Apply a Gaussian filter with sigma=1.0 (adjust as needed)
      M_smoothed = gaussian_filter(M_filled, sigma=1.0)
      
      # Restore the NaNs so they stay blank in the plot
      M_smoothed[nan_mask] = np.nan
      
      M = M_smoothed

  # Tile the matrix 2x2 to show multiple cycles in each dimension
    M_repeated = np.tile(M, (2, 2))
  
    # Plot
    plt.figure(figsize=(8,6))
    # extent=[0,2,0,2] means the x-axis runs 0..2 cycles & y-axis 0..2 cycles
    plt.imshow(
      M_repeated,
      origin='lower',
      aspect='auto',
      cmap='magma',
      extent=[0, 2, 0, 2]
  )
    freq1_val = float(freq1)  # or freq1[0]
    freq2_val = float(freq2)
    plt.rcParams.update({
    "font.size": 20,          # Overall font size
    "axes.labelsize": 20,      # Axis label font size
    "xtick.labelsize": 20,     # X tick label font size
    "ytick.labelsize": 20,     # Y tick label font size
    "legend.fontsize": 20      # Legend font size
    })
    plt.colorbar(label="Flux Value")
    plt.xlabel(f"{freq1_val:.3f} c/d Phase")
    plt.ylabel(f"{freq2_val:.3f} c/d Phase")
    #plt.title("2D Phase Folded Flux (Smoothed)" if smooth else "2D Phase Folded Flux")
    plt.show()
  
    return M  # or return M_repeated, etc.
    
def spectrum_data():
    data = np.loadtxt("C:/Users/jakem/OneDrive/Documents/Year 4 Project/TIC15853131_spec_HTC/TIC15853131_spec.dat")


    wavelength = data[:, 0]
    flux       = data[:, 1]
    plt.figure(figsize=(10, 5))
    plt.plot(wavelength, flux * 10**15, color='dodgerblue', linewidth=1)

    # Set axis limits
    plt.ylim(0.05, 0.4)
    plt.xlim(4000, 7500)

    # Label axes
    plt.xlabel(r"Wavelength (Å)", fontsize=20)
    plt.ylabel(r"$F_{\lambda}$ [10$^{-15}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$]", fontsize=20)

    # === Annotated Lines ===
    lines = {
        'Hβ': 4861,
        'HeII–?': 4686,
        'Hα': 6563
    }

    for label, wavelength_val in lines.items():
        plt.axvline(x=wavelength_val, color='r', linestyle=':', linewidth=1)
        plt.text(wavelength_val + 20, 0.38, label, rotation=90, va='top', ha='left', fontsize=10)

    plt.tight_layout()
    plt.show()
    return wavelength, flux
    
c     = 2.99792458e10   # speed of light in cm/s
me    = 9.10938356e-28  # electron mass in g
e_cgs = 4.80320425e-10  # electron charge in statcoulomb
kB    = 1.380649e-16    # Boltzmann's constant in erg/K
hplanck = 6.62607015e-27  # Planck's constant in erg*s

def planck_I_nu(nu, T):
    """
    Planck function for specific intensity in cgs:
      B_nu(T) = (2h nu^3 / c^2) * 1/(exp(h nu / kB T) - 1).
    Returns erg/s/cm^2/Hz/ster
    """
    x = hplanck*nu/(kB*T)
    return (2.0*hplanck*nu**3 / c**2) / (np.exp(x) - 1.0)

def rj_I_nu(nu, T):
    """
    Rayleigh-Jeans approximation to Planck's law in cgs:
      I_RJ = 2 nu^2 k_B T / c^2
    Good for h nu << k_B T.
    """
    return 2.0 * kB * T * (nu**2) / (c**2)

def cyclotron_freq(B):
    """
    Cyclotron (angular) frequency:
      omega_c = e B / (m_e c)
    Returns omega_c in rad/s for B in Gauss.
    """
    return e_cgs * B / (me * c)

def plasma_freq(n_e):
    """
    Plasma (angular) frequency:
      omega_p = sqrt(4 pi n_e e^2 / m_e)
    You may or may not expose n_e as a parameter. 
    """
    return np.sqrt(4.0 * np.pi * n_e * e_cgs**2 / me)


def alpha_mode(nu, B, Te, theta, mode='+'):
    """
    Returns absorption coeff alpha_+ or alpha_- in cgs: cm^-1
    based on the complicated expressions in the references.
    """
    
    if mode=='+':
        n = 1+2
    # TODO: Insert actual multi-harmonic expansions, Bessel functions, etc.
    # For now, we just pretend there's some function alpha(nu).
    # In practice you'd compute the infinite sum in e.g. equation (5).
    alpha_val = 1e-5  # placeholder
    return alpha_val


def j_mode(nu, B, Te, theta, mode='+'):
    """
    Emissivity for the + or - mode in cgs: erg/s/cm^3/ster/Hz
    Again this is a big expression with sums over s, Bessel functions, etc.
    """
    # Use the full integral for j_{±} ~ ∫ n_{±}(p,theta) f(p) ...
    j_val = 1e-10  # placeholder
    return j_val


def cyclotron_spectrum_model(nu_array, B, Te, theta, Lambda,
                             use_planck=True):

    # Choose which “blackbody” to use as the source function
    if use_planck:
        S_nu = planck_I_nu
    else:
        S_nu = rj_I_nu

    # Prepare output array
    I_tot = np.zeros_like(nu_array, dtype=float)

    # Loop over frequencies
    for i, nu in enumerate(nu_array):
        
        # 1) The "single-mode" source intensity (RJ or Planck):
        I_source = S_nu(nu, Te)  # e.g. I_{RJ} or B_nu

        # 2) Compute alpha_+, alpha_- at this frequency
        alpha_p = alpha_mode(nu, B, Te, theta, mode='+')
        alpha_m = alpha_mode(nu, B, Te, theta, mode='-')

        # 3) Form the dimensionless absorption phi_+, phi_- ...
        #    If you wanted your "hat{phi}" to differ from alpha directly,
        #    you could incorporate geometry or dimension, but conceptually:
        phi_p = alpha_p  # or alpha_p * L, etc.
        phi_m = alpha_m

        # 4) Optical depths:
        tau_p = Lambda * phi_p
        tau_m = Lambda * phi_m

        # 5) Emergent intensities:
        I_plus = I_source * (1.0 - np.exp(-tau_p))
        I_minus = I_source * (1.0 - np.exp(-tau_m))

        # 6) Sum of modes
        I_tot[i] = I_plus + I_minus

    return I_tot
    
def fit_wrapper(nu, B, Te, theta, Lambda):
    """Adapter to pass into curve_fit.  Returns model flux at frequencies nu."""
    return cyclotron_spectrum_model(nu, B, Te, theta, Lambda, use_planck=False)



wavelength, flux = spectrum_data()
freqs = c/wavelength

p0 = [1e6, 1e7, np.pi/4, 1.0]  # initial guesses (B, T_e, theta, Lambda)
popt, pcov = curve_fit(fit_wrapper, freqs, flux, p0=p0)

B_fit, Te_fit, theta_fit, Lambda_fit = popt
print("Best-fit B:", B_fit)
print("Best-fit Te:", Te_fit)
print("Best-fit theta:", theta_fit)
print("Best-fit Lambda:", Lambda_fit)


    
peak_frequencies = np.array([8.453,0.27])
    
time,flux,exptime,flux_error = sector_data(1)
phase_fold_binned(time,flux, peak_frequencies)
phase_fold2D(time,flux,peak_frequencies[0], peak_frequencies[1])

spectrum_data()




