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

def Kepler_data(index):
    search_result = lk.search_lightcurve("20:15:35.65390 +37:11:22.3802",radius = 10)
    print(search_result)
    lc = search_result[index].download()
    exptime = search_result.table['exptime'][index]
    sap_lc = lc.SAP_FLUX
    quality_flags = sap_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0  # Keeps only points with a quality flag of 0 (good data)
    
    time = sap_lc.time.value[good_quality_mask]
    flux = sap_lc.flux.value[good_quality_mask]
    print("flux", flux)
    flux_error = sap_lc.flux_err.value[good_quality_mask]  
    plt.figure(figsize = (8,4))
    plt.plot(time,flux,lw=1)
    plt.xlabel("BKJD (days)")
    plt.ylabel("Flux (e/s)")
    plt.title("Kepler K2 Lightcurve")
    #sap_lc_cleaned = sap_lc.remove_outliers()
    #time = sap_lc_cleaned.time.value
    #flux = sap_lc_cleaned.flux.value
    return time,flux,exptime,flux_error

def stitch_flatten(indeces):
    light_curves = []
    times = []
    fluxes = []
    flux_errors = []
    labels = ["sector 41", "sector 55", "sector 82", "sector 71"]   #adjust accordingly
    search_result = lk.search_lightcurve('04 07 4.080 +18 55 37.20', mission='TESS')
    for i, index in enumerate(indeces):
        lc = search_result[index].download()
        processed_lc = lc.flatten()
        processed_lc = processed_lc
        sap_lc = processed_lc.SAP_FLUX
        quality_flags = sap_lc.quality
        flagged_indices = np.where(quality_flags != 0)[0]
        good_quality_mask = quality_flags == 0
        time = sap_lc.time.value[good_quality_mask]
        flux = sap_lc.flux.value[good_quality_mask]
        mean = np.mean(flux)
        flux = flux/mean
        flux_error = sap_lc.flux_err.value[good_quality_mask]  
        light_curves.append(sap_lc)
        times.append(time)
        fluxes.append(flux)
        flux_errors.append(flux_errors)
        plt.plot(time, flux,lw=1, label=labels[i])
        plt.legend()
    lc_collection = LightCurveCollection(light_curves)
    stitched_lc = lc_collection.stitch()
    
    quality_flags = stitched_lc.quality
    flagged_indices = np.where(quality_flags != 0)[0]
    good_quality_mask = quality_flags == 0
    
    flux_stitched = stitched_lc.flux.value[good_quality_mask]
    time_stitched = stitched_lc.time.value[good_quality_mask]
    exptime_stitched = 120 #Adjust Accordingly
    #plt.plot(time_stitched,flux_stitched,lw = 1)
    plt.xlabel("Time (BJD-2457000, days)")
    plt.ylabel("Normalised Flux")
    plt.show()
    
    return flux_stitched,time_stitched,exptime_stitched



def XMM_test():
    with fits.open("C:/Users/jakem/OneDrive/Documents/Year 4 Project/XMM Data/0782060201/3154_0782060201_SCX00000ATS.FIT") as hdul:
        spectrum_data = hdul[1].data
        column_headers = spectrum_data.columns.names
        print("Column Headers:", column_headers)
    
def XMM_time_series():
    with fits.open("C:/Users/jakem/OneDrive/Documents/Year 4 Project/XMM Data/J2015/LightCurve/0744640101/pps/P0744640101PNS003SRCTSR8001.FTZ") as hdul:
        spectrum_data = hdul[1].data
        column_headers = spectrum_data.columns.names
        print("Column Headers:", column_headers)
        time = spectrum_data['TIME']
        time = time-min(time)
        rate = spectrum_data['RATE']  # Source rate
        back = spectrum_data['BACKV']  # Background rate
        err_rate = spectrum_data['ERROR']  # Error in source rate
        err_back = spectrum_data['BACKE']  # Error in background rate
        
            
            
        bin_size = 800  # Desired bin size in seconds
        binned_time, binned_rate, binned_err_rate, binned_back, binned_err_back = rebin_lightcurve(time, rate, err_rate, back, err_back, bin_size)
        binned_source_rate = binned_rate - binned_back
        binned_source_error = np.sqrt(binned_err_rate**2 + binned_err_back**2)
        
        plt.errorbar(binned_time,binned_back,yerr = binned_err_back,fmt = 'o',markersize=2,label = 'Rebinned background data',color = 'black', ecolor = 'red')
        plt.xlabel("Time (s-534962075.616543s)", fontsize=12)
        plt.ylabel("Rate (counts/s)", fontsize=12)
        plt.title(f"Rebinned Background Light Curve (Bin Size = {bin_size}s)", fontsize=12)
        plt.legend(fontsize=12)
        plt.show()
        
        mask = binned_source_rate >= 0  # Keep only points where binned_source_rate >= 0
        binned_time = binned_time[mask]
        binned_source_rate = binned_source_rate[mask]
        binned_source_error = binned_source_error[mask]
        
        # Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(binned_time, binned_source_rate)
        plt.errorbar(binned_time, binned_source_rate, yerr=binned_source_error, fmt='o', markersize=2, label='Rebinned Data', color='black', ecolor='red')
        plt.xlabel("Time (s-534962075.616543s)", fontsize=14)
        plt.ylabel("Rate (counts/s)", fontsize=14)
        plt.title(f"Rebinned Light Curve (Bin Size = {bin_size}s)", fontsize=16)
        plt.legend(fontsize=12)
        #plt.ylim(0,0.7)
        plt.show()
    
    
    
    frequencies = np.logspace(-5, -3, 100000)  # 10^-5 to 1 Hz

    # Lomb-Scargle Periodogram
    ls = LombScargle(binned_time, binned_source_rate)
    power = ls.power(frequencies)
    
    # Plot in log-log scale
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies, power, label="Power Spectrum")
    
    frequencies_to_mark = [1.39e-4, 2.78e-4]
    for freq in frequencies_to_mark:
        plt.axvline(x=freq, color="red", linestyle="--", label=f"f = {freq} Hz")
        
    
    plt.xscale("log")
    plt.xlabel("Frequency (Hz)", fontsize=14)
    plt.ylabel("Power", fontsize=14)
    plt.title("Log-Log Periodogram", fontsize=16)
    plt.legend(fontsize=12)
    plt.show()
    
    return binned_source_rate, binned_time

            
def rebin_lightcurve(time, rate, error, back, err_back, bin_size):
    # Calculate the number of bins
    min_time = np.min(time)
    max_time = np.max(time)
    num_bins = int(np.ceil((max_time - min_time) / bin_size))
    
    # Define the bins
    bin_edges = np.linspace(min_time, max_time, num_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Initialize arrays for binned data
    binned_rate = []
    binned_error = []
    binned_back = []
    binned_err_back = []
    
    # Loop through bins and calculate means
    for i in range(num_bins):
        mask = (time >= bin_edges[i]) & (time < bin_edges[i + 1])
        if np.any(mask):
            binned_rate.append(np.mean(rate[mask]))
            binned_error.append(np.sqrt(np.sum(error[mask]**2)) / np.sum(mask))  # Propagated error
            binned_back.append(np.mean(back[mask]))
            binned_err_back.append(np.sqrt(np.sum(err_back[mask]**2)) / np.sum(mask))  # Propagated error
        else:
            # Fill empty bins with NaNs
            binned_rate.append(np.nan)
            binned_error.append(np.nan)
            binned_back.append(np.nan)
            binned_err_back.append(np.nan)
    
    return bin_centers, np.array(binned_rate), np.array(binned_error), np.array(binned_back), np.array(binned_err_back)


def ZTF_data():
    zquery = query.ZTFQuery()
    lcq = lightcurve.LCQuery.from_position(61.767, +18.927, 5)
    ZTF_data = pd.DataFrame({'JD' : lcq.data.mjd+2400000.5, 'Magnitude' : lcq.data.mag, 'Magnitude_Error' : lcq.data.magerr, "Filter" : lcq.data.filtercode})
    #data = lcq.download_data()
    #lcq = lightcurve.LCQuery(data)
    df = ZTF_data
    filter_list = ["zg", "zr", "zi"]
    colour = ["teal", "red", "black"]
    for filter, colour in zip(filter_list, colour):     
        plt.scatter(df.JD[df.Filter == filter]-2400000.5, df.Magnitude[df.Filter == filter], color=colour, label=filter,s=1)
    plt.gca().invert_yaxis()
    plt.legend() 
    plt.xlabel("Time, BJD")
    plt.ylabel("Magnitude")
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
    
        
def multiple_LC_plot(index_list):
    results = {}

    # Iterate over the index list and call the original function
    for idx in index_list:
        time, flux, exptime,flux_error = sector_data(idx)
        
        # Storing the 3 results as a tuple in a dictionary for easy access
        results[f"var_{idx}_1_2_3"] = (time, flux, exptime)
    times = [value[0] for value in results.values()]
    fluxes = [value[1] for value in results.values()]
    
    # Plotting
    plt.figure(figsize=(10, 6))
    
    # Use a color map to assign different colors to each sector
    colors = plt.cm.viridis(np.linspace(0, 1, len(index_list)))
    
    # Iterate through each sector and plot with a different color
    for i, (time, flux) in enumerate(zip(times, fluxes)):
        plt.plot(time, flux, lw=1, color=colors[i], label=f'Sector {index_list[i]}')
    
    plt.xlabel('Time (MJD)')
    plt.ylabel('Flux (e/s)')
    plt.title('Light Curve')
    plt.legend(loc='best', fontsize='small')  # Add a legend to identify sectors
    plt.show()
    
    
        
    
    
    

def frequency_range(time,flux,del_t):
    N = len(time)
    seconds_per_day = 86400
    min_freq = 3/(N*del_t)*seconds_per_day
    max_freq = 1/(2*del_t)*seconds_per_day
    return min_freq,max_freq

def Lomb_Scargle(time,flux,exptime):
    # Specify the desired frequency range
    min_freq, max_freq = frequency_range(time,flux,exptime)
    # Compute the Lomb-Scargle Periodogram within the specified frequency range
    num_frequency_points = 100000  # You can adjust this based on the desired resolution
    frequency = np.linspace(min_freq, max_freq, num_frequency_points)
    ls = LombScargle(time, flux)
    power = ls.power(frequency,normalization = "model")# Manually compute power without autopower
    print(power)
    # Plot the Lomb-Scargle Periodogram
    plt.figure(figsize=(10, 6))
    plt.plot(frequency, power*frequency, 'k', lw=1)
    
    # Set logarithmic scale for frequency and power if needed
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.yscale('linear')
    #plt.xscale('linear')
    plt.xlim(50,70)
    plt.ylim(0,0.06)
    
    # Set axis labels
    plt.xlabel('Frequency (c/d)')
    plt.ylabel('Power x Frequency')
    
    # Set plot title
    plt.title('Lomb-Scargle Periodogram (Specified Frequency Range)')
    
    
    #Plot peaks#
    peak_frequencies, peak_powers = peak_finder(frequency, power)
    #print("Peak Frequencies", peak_frequencies)
    orbital_frequencies, spin_frequencies,new_frequencies, beat_frequencies = peak_classification(frequency,power,peak_frequencies,peak_powers)
    for freq2 in new_frequencies:
        alrm = false_alarm(ls, power,frequency, freq2)
        #print("Freq", freq2, ":", alrm)
        
    #freq_int_manual2 = [0.27,0.88,1.788,2.178,6.68,8.45,24.54,25.37,25.85]
    #freq_int_manual2 = [2.43,4.08,6.51,7.663,12.23,16.31,20.38]
    freq_int_manual2 = [0.89,6.68,8.455,24.54,25.38,25.86]
    
    y_vals = np.linspace(0,13,1000)
    for freq in orbital_frequencies:
        x_vals = np.linspace(freq,freq,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spin_frequencies:
        x_vals = np.linspace(freq1,freq1,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beat_frequencies:
        x_vals = np.linspace(freq2,freq2,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'black')
    for freq3 in freq_int_manual2:
        x_vals = np.linspace(freq3,freq3,1000)
        plt.plot(x_vals,y_vals,linestyle = ':', color = 'green')
    #print("The remaining frequency peaks are", remaining_frequencies)
    
    
    plt.plot([], [], linestyle=":", color='blue', label='Orbital Frequencies')  # Add one blue line to the legend
    plt.plot([], [], linestyle=":", color='red', label='Spin Frequencies')      # Add one red line to the legend
    plt.plot([], [], linestyle=":", color='green', label='Remaining Frequencies')
    plt.plot([], [], linestyle=":", color='black', label='Beat Frequencies')
    plt.legend()
    # Show the plot
    plt.show()
    return frequency,power,orbital_frequencies,spin_frequencies,new_frequencies,beat_frequencies

def bootstrap_lomb_scargle(time, flux, num_iterations=1000, num_freqs=10000,exptime = 120):
    # Frequency range
    min_freq, max_freq = frequency_range(time, flux, exptime)
    frequencies = np.linspace(min_freq, max_freq, num_freqs)
    power_matrix = np.zeros((num_iterations, len(frequencies)))

    for i in range(num_iterations):
        # Scramble the flux data
        scrambled_flux = np.random.permutation(flux)
        ls = LombScargle(time, scrambled_flux)
        power = ls.power(frequencies, normalization = 'model')
        power_matrix[i, :] = power
    print("Sample power matrix values:\n", power_matrix[:5, :5])
    
    cutoff_powers = []
    for col_idx in range(power_matrix.shape[1]):
        powers = power_matrix[:, col_idx]
        cutoff_power = np.percentile(powers, 99.7)  # 99.7% cutoff
        cutoff_powers.append(cutoff_power)
    cutoff_powers = np.array(cutoff_powers)
    for col_idx in range(min(5, len(frequencies))):  # Plot for the first 5 frequencies
        powers = power_matrix[:, col_idx]
        plt.figure()
        plt.hist(powers, bins=70, alpha=0.7, edgecolor='k')
        plt.axvline(cutoff_powers[col_idx], color='red', linestyle='dashed', label=f"99.7% cutoff: {cutoff_powers[col_idx]:.2f}")
        plt.title(f"Distribution of Powers for Frequency {frequencies[col_idx]:.2f}")
        plt.xlabel("Power")
        plt.ylabel("Number of Occurrences")
        plt.legend()
        plt.show()
        
    ls = LombScargle(time, flux)
    original_power = ls.power(frequencies, normalization = 'model')
    plt.figure()
    plt.plot(frequencies, cutoff_powers, label="99.7% Cutoff Powers")
    plt.plot(frequencies, original_power, label = "Original Power")
    plt.xlabel("Frequency")
    plt.ylabel("Cutoff Power")
    plt.title("99.7% Cutoff Power vs Frequency")
    plt.xlim(0,60)
    plt.ylim(0,0.005)
    plt.legend()
    plt.show()
    
def gaussian(x, a, mu, sigma):
    return a * np.exp(-0.5 * ((x - mu) / sigma)**2)
    

def bootstrap_errors(time, flux, exptime, num_freqs=10000, num_bootstraps=1000, f_min = 8.3, f_max = 8.5):
    # Stack time and flux for easier resampling
    time_flux_pairs = np.column_stack((time, flux))
    min_freq, max_freq = frequency_range(time, flux, exptime)
    frequencies = np.linspace(min_freq, max_freq, num_freqs)
    
    
    frequencies_short = np.linspace(f_min,f_max,num_freqs)

    
    # Exclude frequencies below 1.5
    valid_frequencies = frequencies[(frequencies >= f_min) & (frequencies <= f_max)]

    # Number of data points
    N = len(flux)
    peak_frequencies = []

    for _ in range(num_bootstraps):
        # Resample with replacement
        resampled_indices = np.random.choice(np.arange(N), size=N, replace=True)
        resampled_pairs = time_flux_pairs[resampled_indices]

        # Separate resampled time and flux
        resampled_time = resampled_pairs[:, 0]
        resampled_flux = resampled_pairs[:, 1]

        # Perform Lomb-Scargle periodogram
        ls = LombScargle(resampled_time, resampled_flux)
        power = ls.power(frequencies_short)
        #plt.plot(valid_frequencies,power)
        #plt.show()

        # Find peak frequency (only considering valid frequencies)
        peak_frequency = frequencies_short[np.argmax(power)]
        peak_frequencies.append(peak_frequency)
        
    bin_heights, bin_edges = np.histogram(peak_frequencies, bins=100)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Calculate bin centers
    
    filter_mask = (bin_centers >= f_min) & (bin_centers <= f_max)
    filtered_bin_centers = bin_centers[filter_mask]
    filtered_bin_heights = bin_heights[filter_mask]
    # Fit the Gaussian to the histogram
    popt, pcov = curve_fit(gaussian, filtered_bin_centers, filtered_bin_heights, p0=[max(bin_heights), 8.45, 0.01])
    
    # Extract the fitted parameters
    a_fit, mu_fit, sigma_fit = popt    
        
    #plt.hist(peak_frequencies, bins=1000, alpha=0.5)
    plt.step(filtered_bin_centers, filtered_bin_heights, where="mid", label="Peak Frequency Distribution", color="blue", alpha=0.7)
    x_fit = np.linspace(f_min,f_max, 1500)
    plt.plot(x_fit, gaussian(x_fit, *popt), label=f"Gaussian Fit\n$\mu$={mu_fit:.5f}, $\sigma$={sigma_fit:.5f}", color="red")
    plt.xlim(f_min,f_max)
    #plt.ylim(0,1500)
    plt.xlabel("Peak Frequency")
    plt.ylabel("Count")
    plt.title("Distribution of Peak Frequencies")
    plt.legend()
    plt.show()

    print(f"Fitted Gaussian Parameters: a={a_fit:.5f}, mu={mu_fit:.5f}, sigma={sigma_fit:.5f}")
    return bin_centers,bin_heights, popt

def bootstrap_errors_multiple(indices,f_min = 8.3,f_max = 8.5):
    # Use lists instead of np.array for initial storage
    centres_array = []
    heights_array = []
    popt_array = []

    # Loop through each index in the dataset
    for i in range(len(indices)):
        # Extract data for the current sector
        time, flux, exptime, flux_error = sector_data(indices[i])

        # Get histogram bin centers, heights, and Gaussian fit parameters
        bin_centers, bin_heights, popt = bootstrap_errors(
            time, flux, exptime,num_freqs = 10000,num_bootstraps = 1000
        )

        # Append results to the lists
        centres_array.append(bin_centers)
        heights_array.append(bin_heights)
        popt_array.append(popt)

    # Plotting multiple histograms and Gaussian fits
    fig, axes = plt.subplots(len(indices), 1, figsize=(8, len(indices) * 3), sharex=True)

    for j in range(len(centres_array)):
        ax = axes[j]

        # Plot stepped histogram
        ax.step(
            centres_array[j],
            heights_array[j],
            where="mid",
            label="Histogram",
            color=f"C{j}",
            alpha=0.7,
        )

        # Plot Gaussian fit
        x_fit = np.linspace(f_min, f_max, 1000)
        a_fit, mu_fit, sigma_fit = popt_array[j]
        ax.plot(
            x_fit,
            gaussian(x_fit, a_fit, mu_fit, sigma_fit),
            label=f"Gaussian Fit\n$\\mu={mu_fit:.3f} \\pm {sigma_fit:.3f}$",
            color=f"C{j}",
        )

        # Add labels and grid
        ax.set_ylabel("RMS Power", fontsize=10)
        ax.legend(fontsize=8)

        # Annotate the mean and uncertainty
        ax.text(
            0.95,
            0.85,
            f"$\\mu={mu_fit:.2f}$\n$\\sigma={sigma_fit:.2f}$",
            transform=ax.transAxes,
            fontsize=8,
            ha="right",
            color=f"C{j}",
        )

    # Finalize plot
    axes[-1].set_xlabel("Frequency [d$^{-1}$]", fontsize=10)
    plt.tight_layout()
    plt.show()
    
    stitch_flatten(indeces)
        
        
    
        
def mulitple_sector_LS(index_list):
    results = {}

    # Iterate over the index list and call the original function
    for idx in index_list:
        time, flux, exptime,flux_error = sector_data(idx)
        
        # Storing the 3 results as a tuple in a dictionary for easy access
        results[f"var_{idx}_1_2_3"] = (time, flux, exptime,flux_error)
    times = [value[0] for value in results.values()]
    fluxes = [value[1] for value in results.values()]
    exptimes = [value[2] for value in results.values()]
    
    results2 = {}
    for i in range (0,len(times)):
        frequency,power,orbital_frequencies,spin_frequencies,new_frequencies, beat_frequencies = Lomb_Scargle(times[i],fluxes[i],exptimes[i])
        results2[f"var_{i}_1_2_3_4_5_6"] = (frequency,power,orbital_frequencies,spin_frequencies,new_frequencies,beat_frequencies)

    frequencies = [value[0] for value in results2.values()]
    powers = [value[1] for value in results2.values()]
    orbitals = [value[2] for value in results2.values()]
    spins = [value[3] for value in results2.values()]
    news = [value[4] for value in results2.values()]
    beats = [value[5] for value in results2.values()]
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(10, 12), sharex=True)
    plt.subplots_adjust(hspace=0)  # hspace=0 removes the space between the subplot
    
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    ax1.set_xlim(10,30)
    ax1.set_ylim(0,0.2)
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    ax2.set_xlim(10,30)
    ax2.set_ylim(0,0.2)
    ax2.set_xlabel('Frequency (c/d)', fontsize = 14)
    fig.text(0, 0.5, 'Power x Frequency', va='center', rotation='vertical', fontsize=14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    
    #Interesting_frequencies#
    freq_int_manual1 = [0.27,0.89,7.567,8.40575]
    freq_int_manual2 = [0.27,0.89,6.68,8.40575]
    #freq_int_manual1 = [6.8267,22.4769, 81.1037,95.761]
    #freq_int_manual2 = [6.8267,22.4769, 81.1037,95.761]
    
    #AX1#
    ax1.plot(frequencies[0], powers[0]*frequencies[0], 'k', lw=1)
    y_vals = np.linspace(0,13,1000)
    for freq in orbitals[0]:
       x_vals = np.linspace(freq,freq,1000)
       ax1.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spins[0]:
      x_vals = np.linspace(freq1,freq1,1000)
      ax1.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beats[0]:
      x_vals = np.linspace(freq2,freq2,1000)
      ax1.plot(x_vals,y_vals,linestyle = ":", color = 'black')
    for freq3 in freq_int_manual1:
        x_vals = np.linspace(freq3,freq3,1000)
        ax1.plot(x_vals,y_vals,linestyle = ":", color = 'green')
    #print("The remaining frequency peaks are", remaining_frequencies)

    
    #AX2#
    ax2.plot(frequencies[1], powers[1]*frequencies[1], 'k', lw=1)
    y_vals = np.linspace(0,13,1000)
    for freq in orbitals[1]:
        x_vals = np.linspace(freq,freq,1000)
        ax2.plot(x_vals,y_vals,linestyle = ":", color = 'blue')
    for freq1 in spins[1]:
        x_vals = np.linspace(freq1,freq1,1000)
        ax2.plot(x_vals,y_vals,linestyle = ":", color = 'red')
    for freq2 in beats[1]:
        x_vals = np.linspace(freq2,freq2,1000)
        plt.plot(x_vals,y_vals,linestyle = ":", color = 'black')
    for freq3 in freq_int_manual2:
        x_vals = np.linspace(freq3,freq3,1000)
        ax2.plot(x_vals,y_vals,linestyle = ":", color = 'green')
    #print("The remaining frequency peaks are", remaining_frequencies)

    
    
    #Legend#
    plt.plot([], [], linestyle=":", color='blue', label='Orbital Frequencies')  # Add one blue line to the legend
    plt.plot([], [], linestyle=":", color='red', label='Spin Frequencies')      # Add one red line to the legend
    plt.plot([], [], linestyle=":", color='green', label='Remaining Frequencies')
    plt.plot([], [], linestyle=":", color='black', label='Beat Frequencies')
    plt.legend()
    
    
    # Show both plots vertically
    plt.tight_layout()
    plt.show()
    
    return times, fluxes, orbitals, spins, news

def Lomb_Scargle_2D(time, flux, exptime, window_size, step_size, min_freq, max_freq):
    # Parameters for the frequency range
    num_frequency_points = 1000  # Increase if higher frequency resolution is needed
    frequencies = np.linspace(min_freq, max_freq, num_frequency_points)

    # Prepare to store the 2D power spectrum (rows: time windows, columns: frequencies)
    power_spectrum_2D = []

    # Sliding window through the data
    time_start = min(time)
    time_end = max(time)
    current_time = time_start

    # Loop through the time windows
    while current_time + window_size <= time_end:
        # Find indices of the data within the current window
        window_mask = (time >= current_time) & (time < current_time + window_size)
        time_window = time[window_mask]
        flux_window = flux[window_mask]

        # Compute the Lomb-Scargle periodogram for the current window
        ls = LombScargle(time_window, flux_window)
        power = ls.power(frequencies)
        power_spectrum_2D.append(power)

        # Move the window
        current_time += step_size

    # Convert the list of power spectra into a 2D array
    power_spectrum_2D = np.array(power_spectrum_2D)

    # Generate the time axis (midpoints of each time window)
    time_axis = np.arange(time_start + window_size / 2, time_end, step_size)

    # Plot the 2D power spectrum
    plt.figure(figsize=(10, 6))
    plt.imshow(power_spectrum_2D.T, aspect='auto', extent=[time_axis[0], time_axis[-1], min_freq, max_freq], origin='lower', cmap='inferno')
    plt.colorbar(label='Power')
    plt.xlabel('Time (BTJD)')
    plt.ylabel('Frequency (cycles/day)')
    plt.title('2D Lomb-Scargle Power Spectrum')
    plt.show()

def peak_finder(frequency, power,  height_threshold=0.01, prominence=0.001):
    y = frequency*power
    peaks, properties = find_peaks(y, height=height_threshold, prominence=prominence)
    
    # Extract the frequencies and powers of the found peaks
    peak_frequencies = frequency[peaks]
    peak_powers = y[peaks]
    
    return peak_frequencies, peak_powers

def peak_classification(frequency,power,peak_frequencies,peak_powers,tolerance = 0.005):
    orbital_period = 0.131
    #orbital_period = 0.068233846
    spin_period = 0.122
    natural_orbital_frequency = 1/orbital_period
    natural_spin_frequency = 1/spin_period
    natural_beat_frequency = abs(natural_spin_frequency-natural_orbital_frequency)
    orbital_frequencies = []
    spin_frequencies = []
    beat_frequencies = []
    for freq in peak_frequencies:
        if abs(freq / natural_orbital_frequency - round(freq / natural_orbital_frequency)) < tolerance:
            orbital_frequencies.append(freq)
        if abs(freq / natural_spin_frequency - round(freq / natural_spin_frequency)) < tolerance:
            spin_frequencies.append(freq)
        if abs(freq / natural_beat_frequency - round(freq / natural_beat_frequency)) < tolerance:
            beat_frequencies.append(freq)
    classified_frequencies = set(orbital_frequencies + spin_frequencies + beat_frequencies)
    
    # Filter out the frequencies that are in the classified frequencies set
    remaining_frequencies = [freq for freq in peak_frequencies if freq not in classified_frequencies]
    tolerance2 = 0.3  # For example, consider frequencies within 0.1 c/d as "close"
    print("spin",spin_frequencies)
    print("orbital", orbital_frequencies)
    print("beat", beat_frequencies)
    # List to store new frequencies
    new_frequencies = []
    
    # Loop through remaining frequencies
    for freq in remaining_frequencies:
        # Check if this frequency is close to any orbital or spin frequency
        is_known = np.any(np.abs(orbital_frequencies - freq) < tolerance2) or np.any(np.abs(spin_frequencies - freq) < tolerance2) or np.any(np.abs(beat_frequencies - freq) < tolerance2)
        
        # If it's not close to any known frequency, add it to new_frequencies
        if not is_known:
            new_frequencies.append(freq)
            
    # Convert new_frequencies to a numpy array for easier handling (optional)
    new_frequencies = np.array(new_frequencies)
            
            # Output the new frequencies
    #print("New frequencies:", new_frequencies)
    return orbital_frequencies, spin_frequencies,new_frequencies, beat_frequencies

def false_alarm(ls,power,frequency,specific_frequency):
    closest_idx = np.argmin(np.abs(frequency - specific_frequency))
    prob = ls.false_alarm_probability(power[closest_idx])
    return prob

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
    
    
def multi_periodic_model(t, params):
    A1, f1, phi1, A2, f2, phi2, A3, f3, phi3, A4, f4, phi4, A5, f5, phi5, offset = params
    component1 = A1 * np.sin(2 * np.pi * f1 * t + phi1)  # Spin period
    component2 = A2 * np.sin(2 * np.pi * f2 * t + phi2)  # Orbital period
    component3 = A3 * np.sin(2 * np.pi * f3 * t + phi3)  # Beat frequency
    component4 = A4 * np.sin(2* np.pi * f4 * t + phi4)
    component5 = A5 * np.sin(2* np.pi * f5 * t + phi5)
    return component1 + component2 + component3 + component4 +  offset

def chi_squared(params, time, flux, flux_error):
    """Chi-squared calculation for the model."""
    A1, f1, phi1, A2, f2, phi2, A3, f3, phi3,A4,f4,phi4,A5,f5,phi5, offset = params
    model_flux = multi_periodic_model(time, A1, f1, phi1, A2, f2, phi2, A3, f3, phi3,A4,f4,phi4,A5,f5,phi5, offset)
    print(np.sum(((flux - model_flux) / flux_error) ** 2))
    return np.sum(((flux - model_flux) / flux_error) ** 2)

def log_likelihood(params, time, flux, flux_error):
    model_flux = multi_periodic_model(time, params)
    chi_squared = np.sum(((flux - model_flux) / flux_error) ** 2)
    return -0.5 * chi_squared

def log_prior(params):
    A1, f1, phi1, A2, f2, phi2, A3, f3, phi3, A4, f4, phi4, A5, f5, phi5, offset = params
    if (
        1 <= A1 <= 10
        and 11.9 <= f1 <= 12.1
        and 0 <= phi1 <= 2 * np.pi
        and 2 <= A2 <= 10
        and 1.8 <= f2 <= 1.88
        and 0 <= phi2 <= 2 * np.pi
        and 1 <= A3 <= 10
        and 3.7 <= f3 <= 3.8
        and 0 <= phi3 <= 2 * np.pi
        and 0 <= A4 <= 10
        and 4.5 <= f4 <= 4.7
        and 0 <= phi4 <= 2 * np.pi
        and 0 <= A5 <= 10
        and 5.26 <= f5 <= 5.29
        and 0 <= phi5 <= 2 * np.pi
        and 18 <= offset <= 21
    ):
        return 0.0  # log(1)
    return -np.inf  # log(0)

def log_probability(params, time, flux, flux_error):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, time, flux, flux_error)

def MCMC_Fit(time,flux,flux_error):
    start_time = 2475
    end_time = 2478
    
    subset_indices = (time >= start_time) & (time <= end_time)
    time_subset = time[subset_indices]
    flux_subset = flux[subset_indices]
    flux_error_subset = flux_error[subset_indices]
    # Initial setup for MCMC
    num_params = 16
    num_walkers = 32
    num_steps = 10000
    initial_guesses = [
        0.5,  # A1
        11.95,  # f1
        2*np.pi,  # phi1
                0.5,  # A2
                1.88,  # f2
                2*np.pi,  # phi2
                4,  # A3
                3.75,  # f3
                np.pi,  # phi3
                2,  # A4
                4.6,  # f4
                np.pi,  # phi4
                1,  # A5
                5.27,  # f5
                np.pi,  # phi5
                20 # offset
                ]
    
    # Initialize the walkers in a Gaussian ball around initial guesses
    pos = initial_guesses + 1 * np.random.randn(num_walkers, num_params)
    
    # Set up the sampler
    sampler = emcee.EnsembleSampler(num_walkers, num_params, log_probability, args=(time_subset, flux_subset, flux_error_subset))
    
    # Run MCMC
    sampler.run_mcmc(pos, num_steps, progress=True)
    
    # Analyze the results
    samples = sampler.get_chain(discard=100, thin=15, flat=True)
    
    # Plotting the results
    import corner
    fig = corner.corner(samples, labels=["A1", "f1", "phi1", "A2", "f2", "phi2", "A3", "f3", "phi3", "A4", "f4", "phi4", "A5", "f5", "phi5", "offset"])
    plt.show()
    
    # Extract median values for each parameter
    params_mcmc = np.median(samples, axis=0)
    param_names = ["A1", "f1", "phi1", "A2", "f2", "phi2", "A3", "f3", "phi3", "A4", "f4", "phi4", "A5", "f5", "phi5", "offset"]
    params_std = np.std(samples, axis=0)  # 1-sigma uncertainty
    
    # Print the results
    print("Final parameter values and uncertainties:")
    for i, name in enumerate(param_names):
        print(f"{name}: {params_mcmc[i]:.4f} ± {params_std[i]:.4f}")
        
    
    
    # Plot the data with the best-fit model
    fitted_flux = multi_periodic_model(time_subset, params_mcmc)
    
    chi_squared = np.sum(((flux_subset - fitted_flux) / flux_error_subset) ** 2)
    reduced_chi_squared = chi_squared / (len(flux_subset) - len(params_mcmc))

    print("Reduced Chi-squared:", reduced_chi_squared)
    
    plt.figure(figsize=(12, 6))
    plt.scatter(time_subset,flux_subset, label = 'observed data', s=5,color = 'black')
    plt.plot(time_subset,flux_subset,lw=1, color = 'blue')
    plt.plot(time_subset, fitted_flux, label="MCMC Fit", color="red")
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.title("Multi-Periodic Model Fit with MCMC")
    plt.legend()
    plt.show()
        
def model_fit(time, flux, flux_error):
    start_time = 2797
    end_time = 2798
    
    subset_indices = (time >= start_time) & (time <= end_time)
    time_subset = time[subset_indices]
    flux_subset = flux[subset_indices]
    flux_error_subset = flux_error[subset_indices]


    # Initial guesses for the parameters
    initial_guesses = [
    10,  # A1: Spin amplitude
    11.997,  # f1: Spin frequency
    -np.pi,    # phi1: Spin phase
    1,   # A2: Orbital amplitude
    1.8, # f2: Orbital frequency
    0,    # phi2: Orbital phase
    5,   # A3: Beat amplitude
    3.76, # f3: Beat frequency
    np.pi,    # phi3: Beat phase
    3, #A4
    4.58, #f4
    -np.pi, #phi4
    2, #A5
    5.28, #f5
    -np.pi, #phi5
    1000  # Offset
]

    # Bounds for the parameters
    bounds = [
    (4,10),      # Bounds for A1
    (11.9,12.1),         # Bounds for f1
    (-np.pi, np.pi), # Bounds for phi1
    (0.5,2),      # Bounds for A2
    (6.5,6.5),       # Bounds for f2
    (-np.pi, np.pi), # Bounds for phi2
    (2,5),       # Bounds for A3
    (3,8),       # Bounds for f3
    (-np.pi, np.pi),# Bounds for phi3
    (2,4),
    (4.5,4.7),
    (-np.pi,np.pi),
    (1,3),
    (5.26,5.29),
    (-np.pi,np.pi),
    (800, 1200)     # Bounds for offset
]

    # Minimize the chi-squared function
    result = minimize(
        chi_squared, initial_guesses, args=(time_subset, flux_subset, flux_error_subset),
        bounds=bounds
    )

    # Extract the fitted parameters
    params = result.x
    fitted_flux = multi_periodic_model(time_subset, *params)
    chi2 = chi_squared(params, time_subset, flux_subset, flux_error_subset)
    reduced_chi2 = chi2 / (len(flux_subset) - len(params))  # Reduced chi-squared

    print("Parameters:", params)
    print("Chi-squared:", chi2)
    print("Reduced Chi-squared:", reduced_chi2)

    # Plot the fit
    plt.figure(figsize=(12, 6))
    plt.plot(time_subset, flux_subset, label='Observed Data', alpha=0.6)
    plt.scatter(time_subset, flux_subset, s=5, color = 'black')
    plt.plot(time_subset, fitted_flux, label='Fitted Model', color='red')
    plt.xlabel("Time")
    plt.ylabel("Flux")
    plt.legend()
    plt.title("Multi-Periodic Model Fit with Chi-squared Minimization")
    plt.show()

    return params, result.hess_inv, fitted_flux, chi2, reduced_chi2

indexes = [0,1]
indeces = [0,1,2]

#flux_stitched, time_stitched, exptime_stitched = stitch_flatten(indeces)
#Lomb_Scargle(time_stitched, flux_stitched, exptime_stitched)


#XMM_test()
#flux,time = XMM_time_series()
#ZTF_data()

#time,flux,exptime, flux_error = Kepler_data(indexes[0])
#Lomb_Scargle(time,flux,exptime)


#multiple_LC_plot(indexes)
#times, fluxes, orbitals, spins, news = mulitple_sector_LS(indexes)
#peak_frequencies = np.array([11])
#phase_fold_binned(time, flux, peak_frequencies)


time,flux,exptime,flux_error = sector_data(indeces[1])
#bootstrap_lomb_scargle(time, flux)
bootstrap_errors(time, flux, exptime)
bootstrap_errors_multiple(indeces)
#Lomb_Scargle(time, flux, exptime)

#frequency,power,peak_frequencies,peak_powers = Lomb_Scargle(time,flux,exptime)
#peak_classification(frequency,power,peak_frequencies,peak_powers)
window_size = 0.3  # Window size in the same units as time (e.g., days or minutes)
step_size = 5  # Step size for sliding the window
min_freq = 5  # Minimum frequency (cycles/day)
max_freq = 40  # Maximum frequency (cycles/day)

# Call the function with your time, flux, and exptime data
#Lomb_Scargle_2D(time, flux, exptime, window_size, step_size, min_freq, max_freq)
#model_fit(time,flux,flux_error)
#MCMC_Fit(time,flux,flux_error)