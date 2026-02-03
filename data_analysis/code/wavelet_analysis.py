"""
HPA Axis Nidek Wavelet Analysis

Performs wavelet analysis on HPA axis hormone data (Cortisol/ACTH) to identify
periodic patterns and ultradian rhythms in the time series.

This script uses PyBOAT (Python Biological Oscillations Analysis Toolkit) to:
1. Detrend and normalize hormone signals
2. Compute wavelet spectra
3. Extract dominant frequency ridges
4. Analyze power distributions across participants and time

Usage:
    python data_analysis/code/wavelet_analysis.py --hormone Cortisol
    python data_analysis/code/wavelet_analysis.py --hormone ACTH --dt 10

Note that in the manuscript we used dt 10 for Cortisol (which means the file needs to be created first in data_prep.py)   

Requirements:
    - PyBOAT (pip install pyboat-dtk)
    - data_analysis/options/participant_options.csv with wavelet parameters

Input:
    - data_analysis/data/data_shifted_{hormone}_{dt}T_00_00.csv
    - data_analysis/options/participant_options.csv -- These were the options chosen for each participant.

Output:
    - output/stacked_wavelets_{hormone}.npy
    - output/periods.npy
    - output/times.npy
    - Multiple PNG/PDF plots in output/ directory
"""

import os
import argparse
import numpy as np
import pandas as pd
from pyboat import WAnalyzer, ssg
from pyboat import plotting as pl
from pyboat.core import ar1_powerspec, sliding_window_amplitude
import matplotlib.pyplot as ppl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta, datetime
import seaborn as sns
from scipy.interpolate import interp1d
import matplotlib.dates as matdates
from matplotlib.ticker import MultipleLocator
from scipy.stats import mannwhitneyu, kruskal


def moving_average(data, window_size=45):
    """Apply moving average smoothing to data."""
    return np.convolve(data, np.ones(window_size) / window_size, mode='same')


def convert_timesteps_to_times(ref_time, timesteps):
    """Convert timesteps to time strings."""
    base_time = datetime.strptime(ref_time, "%H:%M")
    times = [base_time + timedelta(minutes=step) for step in timesteps]
    return [time.strftime("%H:%M") for time in times]


def run_pyboat(signal, row_options, cutoff=1440, plot_label="", dt=10, hormone='Cortisol', results_folder='output', save=True):
    """
    Run PyBOAT wavelet analysis on a single participant's signal.
    
    Args:
        signal: Time series signal to analyze
        row_options: DataFrame row with analysis parameters
        cutoff: Detrending cutoff period (minutes)
        plot_label: Label for plots
        dt: Sampling interval (minutes)
        hormone: Hormone name
        results_folder: Output directory
        save: Whether to save plots
    
    Returns:
        tuple: (modulus, ridge, periods)
    """
    min_period = row_options['min_period'].values[0]
    max_period = row_options['max_period'].values[0]
    power_thresh = row_options['power_thresh'].values[0]
    smoothing_wsize = row_options['smoothing_wsize'].values[0]
    window_size = row_options['sliding_window'].values[0]
    
    periods = np.linspace(min_period, max_period, 1000)

    wAn = WAnalyzer(periods, dt, time_unit_label='min')

    # Calculate trend and detrend signal
    trend = wAn.sinc_smooth(signal, T_c=cutoff)
    detrended_signal = signal - trend

    # Plot original signal vs trend
    ppl.figure(figsize=(10, 6))
    ppl.plot(signal, label='Signal', color='blue', alpha=0.6)
    ppl.plot(trend, label='Trend', color='green', alpha=0.8)
    ppl.legend()
    ppl.title(f'Original signal and trend {plot_label}')
    ppl.xlabel('Clocktime')
    ppl.ylabel('Values')
    if save: 
        ppl.savefig(f'{results_folder}/trend_{plot_label}_{hormone}.png', dpi=150, bbox_inches='tight')
        ppl.close()
    
    # Normalize the amplitude with a sliding window
    norm_signal = wAn.normalize_amplitude(detrended_signal, window_size=window_size)

    # Compute wavelet spectrum
    modulus, transform = wAn.compute_spectrum(norm_signal)

    wAn.ax_spec_signal.set_title(f'{plot_label} Signal {hormone}')
    wAn.ax_spec.grid(axis='both', color='white', alpha=0.4)
    
    # Extract ridge
    ridge = wAn.get_maxRidge(power_thresh=power_thresh, smoothing_wsize=smoothing_wsize)
    wAn.draw_Ridge()
    
    # Format x-axis
    ax = ppl.gca()
    start_time = '00:00:00'
    start_time = datetime.strptime(start_time, "%H:%M:%S")
    max_time = 1440
    ax.set_xticks([i for i in range(180, max_time, 360)])
    ax.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(180, max_time, 360)])
    wAn.ax_spec_signal.set_title('')
    
    if save: 
        ppl.savefig(f'{results_folder}/wavelet_{plot_label}_{hormone}.pdf', bbox_inches='tight')
        ppl.savefig(f'{results_folder}/wavelet_{plot_label}_{hormone}.png', dpi=150, bbox_inches='tight')
        ppl.close()
    
    return modulus, ridge, periods


def plot_participant_total_power(data, periods, hormone='Cortisol', results_folder='output', save=True):
    """Plot total wavelet power for each participant."""
    calculated_results = [np.sum(participant, axis=1) for participant in data]

    with ppl.rc_context({
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
    }):
        fig, axs = ppl.subplots(2, 5, figsize=(15, 10), sharey=True)

        for i, total_power in enumerate(calculated_results):
            ax = axs[i // 5, i % 5]
            ax.plot(periods[i], total_power, color='black')
            ax.set_title(f'Participant {i + 1}')
            ax.set_xlabel('Periods')
            if i % 5 == 0:
                ax.set_ylabel('Total wavelet power')
            ax.set_xticks([100, 200, 300])
            ax.tick_params(axis='both', which='both')

        ppl.tight_layout()

        if save:
            ppl.savefig(f'{results_folder}/wavelet_total_power_period_{hormone}.png', dpi=300, bbox_inches='tight')
            ppl.close()

    return calculated_results


def create_heatmaps(stacked_wavelet, periods, times, hormone='Cortisol', results_folder='output'):
    """Create mean and median heatmaps of wavelet power."""
    x_ticks_steps = 36
    y_ticks_steps = 50

    mean_wavelet_power = np.mean(stacked_wavelet, axis=0)
    median_wavelet_power = np.median(stacked_wavelet, axis=0)

    # Mean heatmap
    ppl.figure(figsize=(10, 6))
    ax = sns.heatmap(mean_wavelet_power, cmap='viridis', cbar=True)
    ax.set_xticks(range(0, 144, x_ticks_steps))
    ax.set_xticklabels(times[::x_ticks_steps])
    ax.set_yticks(range(0, 1000, y_ticks_steps))
    ax.set_yticklabels(periods[::y_ticks_steps].round(0))
    ax.invert_yaxis()
    ppl.title('Mean heatmap for wavelet powers')
    ppl.savefig(f'{results_folder}/mean_heatmap_{hormone}.png', dpi=150, bbox_inches='tight')
    ppl.close()

    # Median heatmap
    ppl.figure(figsize=(10, 6))
    ax = sns.heatmap(median_wavelet_power, cmap='viridis', cbar=True)
    ax.set_xticks(range(0, 144, x_ticks_steps))
    ax.set_xticklabels(times[::x_ticks_steps])
    ax.set_yticks(range(0, 1000, y_ticks_steps))
    ax.set_yticklabels(periods[::y_ticks_steps].round(0))
    ax.invert_yaxis()
    ppl.title('Median heatmap for wavelet powers')
    ppl.savefig(f'{results_folder}/median_heatmap_{hormone}.png', dpi=150, bbox_inches='tight')
    ppl.close()

    return mean_wavelet_power, median_wavelet_power


def plot_heatmap_summaries(mean_wavelet_power, median_wavelet_power, periods, hormone='Cortisol', results_folder='output', window_size=45):
    """Plot summary statistics for heatmaps."""
    percentiles = [5, 25, 50, 75, 90]
    
    # Mean heatmap summary
    percentiles_values = np.percentile(mean_wavelet_power, percentiles, axis=1)
    mean_values = np.mean(mean_wavelet_power, axis=1)
    
    percentiles_values_smoothed = [moving_average(pv, window_size) for pv in percentiles_values]
    mean_values_smoothed = moving_average(mean_values, window_size)
    
    ppl.figure(figsize=(10, 6))
    for j, percentile in enumerate(percentiles):
        ppl.plot(periods, percentiles_values_smoothed[j], label=f'{percentile}th percentile')
    ppl.plot(periods, mean_values_smoothed, label='Mean', linestyle='--')
    ppl.title("Summary statistics of stacked heatmap for mean wavelet power")
    ppl.xlabel("Periods")
    ppl.ylabel("Wavelet power")
    ppl.legend(loc='upper right')
    ppl.tight_layout()
    ppl.savefig(f'{results_folder}/mean_heatmap_summary_{hormone}.png', dpi=150, bbox_inches='tight')
    ppl.close()

    # Median heatmap summary
    percentiles_values = np.percentile(median_wavelet_power, percentiles, axis=1)
    mean_values = np.mean(median_wavelet_power, axis=1)
    
    percentiles_values_smoothed = [moving_average(pv, window_size) for pv in percentiles_values]
    mean_values_smoothed = moving_average(mean_values, window_size)
    
    ppl.figure(figsize=(10, 6))
    for j, percentile in enumerate(percentiles):
        ppl.plot(periods, percentiles_values_smoothed[j], label=f'{percentile}th percentile')
    ppl.plot(periods, mean_values_smoothed, label='Mean', linestyle='--')
    ppl.title("Summary statistics of stacked heatmap for median wavelet power")
    ppl.xlabel("Periods")
    ppl.ylabel("Wavelet power")
    ppl.legend(loc='upper right')
    ppl.tight_layout()
    ppl.savefig(f'{results_folder}/median_heatmap_summary_{hormone}.png', dpi=150, bbox_inches='tight')
    ppl.close()


def plot_mean_powers(stacked_wavelet, periods, hormone='Cortisol', results_folder='output'):
    """Plot mean wavelet power for each participant."""
    plt.rcParams.update({'font.size': 24})
    plt.figure(figsize=(10, 8))

    mean_powers = [np.mean(dist, axis=1) for dist in stacked_wavelet]
    
    for i, mean_power in enumerate(mean_powers):
        plt.plot(periods, mean_power, label=f'Participant {i+1}')
    
    plt.xlabel('Period (minutes)')
    plt.grid(True)
    plt.ylabel('Mean wavelet power')
    plt.xlim(0, 300)
    plt.savefig(f"{results_folder}/mean_wavelet_powers_{hormone}.png", dpi=300, bbox_inches="tight")
    plt.close()


def main():
    """Main wavelet analysis pipeline."""
    parser = argparse.ArgumentParser(description='Run wavelet analysis on HPA axis data')
    parser.add_argument('--hormone', default='Cortisol', help='Hormone to analyze (Cortisol or ACTH)')
    parser.add_argument('--dt', type=int, default=10, help='Sampling interval in minutes')
    parser.add_argument('--data_folder', default='data_analysis/data', help='Data directory')
    parser.add_argument('--results_folder', default=None, help='Output directory')
    args = parser.parse_args()

    hormone = args.hormone
    dt = args.dt
    data_folder = args.data_folder
    # Set results_folder to include hormone and dt if not provided
    if args.results_folder is None:
        results_folder = f"data_analysis/output/{hormone}/{dt}min"
    else:
        results_folder = args.results_folder

    # Create output directory
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    print("="*60)
    print(f"HPA Axis Wavelet Analysis - {hormone}")
    print("="*60)
    print(f"Sampling interval: {dt} minutes")
    print(f"Output directory: {results_folder}\n")

    # Set plot aesthetics
    ppl.rcParams.update({'font.size': 18})
    plt.rcParams.update({'font.size': 18})

    # Load participant options
    options = pd.read_csv('data_analysis/options/participant_options.csv')
    print(f"Sliding window stats:")
    print(f"  Std: {options['sliding_window'].std():.2f}")
    print(f"  Min: {options['sliding_window'].min():.2f}")
    print(f"  Max: {options['sliding_window'].max():.2f}")
    print(f"  Mean: {options['sliding_window'].mean():.2f}\n")

    # Load data
    ref_time = '00:00'
    data_file = f"{data_folder}/data_shifted_{hormone}_{dt}min_{ref_time.replace(':', '_')}.csv"
    print(f"Loading data from: {data_file}")
    try:
        test_data = pd.read_csv(data_file)
    except FileNotFoundError:
        print(f"\nERROR: Data file '{data_file}' not found.")
        print("To generate this file, run:")
        print("    python data_analysis/code/data_prep.py")
        print("and ensure the correct interpolation and reference time options are set in data_prep.py.")
        exit(1)

    # Setup time labels
    timesteps = list(range(0, 1430, 10))
    times = convert_timesteps_to_times(ref_time, timesteps)
    tick_labels = ["03:00", "09:00", "15:00", "21:00"]
    tick_positions = [0, 360, 720, 1080]

    # Run wavelet analysis for each participant
    print("\nRunning wavelet analysis for each participant...")
    ppl.ion()
    
    stacked_wavelet = []
    participant_periods = []
    ridges = []

    for test_value in range(1, 11):
        print(f"  Processing participant {test_value}...")
        signal = test_data[test_data['test'] == test_value]['y'].to_list()
        row_options = options.loc[options['participant_id'] == test_value]
        
        modulus, ridge, periods = run_pyboat(
            signal, row_options, 
            plot_label=f'Participant_{test_value}', 
            dt=dt, 
            hormone=hormone,
            results_folder=results_folder
        )
        
        stacked_wavelet.append(modulus)
        participant_periods.append(periods)
        ridges.append(ridge)

    # Save results
    print("\nSaving wavelet analysis results...")
    np.save(f'{results_folder}/stacked_wavelets_{hormone}.npy', stacked_wavelet)
    np.save(f'{results_folder}/periods.npy', periods)
    np.save(f'{results_folder}/times.npy', times)
    print(f"  Saved: {results_folder}/stacked_wavelets_{hormone}.npy")
    print(f"  Saved: {results_folder}/periods.npy")
    print(f"  Saved: {results_folder}/times.npy")

    # Generate plots and analyses
    print("\nGenerating summary plots...")
    plot_participant_total_power(stacked_wavelet, participant_periods, hormone, results_folder)
    
    mean_wavelet_power, median_wavelet_power = create_heatmaps(
        stacked_wavelet, periods, times, hormone, results_folder
    )
    
    plot_heatmap_summaries(
        mean_wavelet_power, median_wavelet_power, periods, hormone, results_folder
    )
    
    plot_mean_powers(stacked_wavelet, periods, hormone, results_folder)

    print("\n" + "="*60)
    print("Wavelet analysis complete!")
    print("="*60)
    print(f"\nAll outputs saved to: {results_folder}/")


if __name__ == "__main__":
    main()
