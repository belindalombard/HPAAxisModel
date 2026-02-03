import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from ddeint import ddeint
import numpy as np
import time
import pandas as pd
import json
from model.code.classes.model_class import Model
import argparse
from model.code.additional_functions.old import model_plotting_functions
from model.code.additional_functions import metrics_calculations as mc
import csv
from model.code.additional_functions import additional_functions as af
from typing import Optional

def save_metrics_csv(metrics, output_dir):
    import csv, os
    csv_file = os.path.join(output_dir, "metrics.csv")
    file_exists = os.path.isfile(csv_file)
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=metrics.keys())
        if not file_exists:
            writer.writeheader()
        writer.writerow(metrics)

def save_numpy_arrays(output_dir, arrays_dict):
    import numpy as np, os
    for name, arr in arrays_dict.items():
        np.save(os.path.join(output_dir, f"{name}.npy"), arr)

def build_output_dir(scenario: str, stressor_parameters: dict) -> str:
    """Return output directory path under stressors/output/<scenario>/<param_string>."""
    param_str = (
        f"mag{stressor_parameters['magnitude']:.2f}_"
        f"dur{stressor_parameters['duration']:.2f}_"
        f"plat{stressor_parameters.get('plateau', 0):.2f}_"
        f"decay{stressor_parameters.get('decay_shape', 0):.2f}_"
        f"beta{stressor_parameters.get('beta', 0):.2f}_"
        f"start{stressor_parameters['start']:.2f}"
    )
    # Determine base directory relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(script_dir, "output", scenario, param_str)
def save_stressor_only_plot(
    output_dir: str,
    stressor_series_full: np.ndarray,
    length_plot: int,
    num_days: int,
    stressor_parameters: dict,
    base_clock: str = "09:00:00",
    filename: str = "stressor_only.png",
    max_minutes: Optional[int] = 2300,
) -> str:
    """Save a simple 1-axis plot of the stressor-only contribution added to CRH.

    Uses the last `length_plot` minutes (local window), capped at 2300 minutes for layout.
    Marks the surgery time and duration (rise + plateau by default unless overridden).
    """
    # Slice to the last window
    start_ix = length_plot * (num_days - 1)
    series_local = np.asarray(stressor_series_full)[start_ix:start_ix + length_plot].astype(float)
    local_len = len(series_local)
    if local_len == 0:
        print("Warning: No stressor-only data available for plotting.")
        return ""

    # Cap to maximum if provided
    if max_minutes is None or max_minutes <= 0:
        plot_upto = local_len
    else:
        plot_upto = min(int(max_minutes), local_len)
    time_points = np.arange(plot_upto, dtype=float)
    series = series_local[:plot_upto]

    # Build the figure
    plt.rcParams.update({"font.size": 11})
    fig = plt.figure(figsize=(8, 2.5))
    ax = plt.gca()
    ax.plot(time_points, series, color="black", label="Stressor-only (CRH add-on)")

    # Stressor timing and duration markers (local coordinates)
    stress_time_local = float(stressor_parameters.get("time_in_scope", 0))
    rise_duration = float(stressor_parameters.get("duration", 0))
    plateau_duration = float(stressor_parameters.get("plateau", 0))
    surgery_override = stressor_parameters.get("surgery_duration") or stressor_parameters.get("surgery_length")
    surgery_duration = float(surgery_override) if surgery_override is not None else (rise_duration + plateau_duration)
    # Determine label based on scenario
    scenario_label = stressor_parameters.get('scenario', None)
    # If scenario is not in stressor_parameters, try to get from global variable
    if scenario_label is None:
        scenario_label = globals().get('scenario', None)
    if scenario_label is not None and str(scenario_label).upper() == 'ACUTE':
        duration_label = "Stressor duration"
        time_label = "Time of stressor"
    else:
        duration_label = "Surgery duration"
        time_label = "Time of surgery"
    if 0 <= stress_time_local <= plot_upto:
        ax.axvline(stress_time_local, color="red", linestyle="--", label=time_label)
        span_start = max(0, stress_time_local)
        span_end = min(stress_time_local + surgery_duration, plot_upto)
        if span_end > span_start:
            ax.axvspan(span_start, span_end, color="red", alpha=0.2, label=duration_label)

    # Labels and ticks
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("Stressor CRH add-on (a.u.)")
    base_dt = datetime.strptime(base_clock, "%H:%M:%S")
    ticks = np.arange(0, (time_points[-1] if len(time_points) else 0) + 1, 360)
    labels = [(base_dt + timedelta(minutes=float(t))).strftime("%H:%M") for t in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)

    # Y bounds with a little headroom
    if len(series):
        ymax = float(np.nanmax(series))
        ax.set_ybound(0, ymax + max(1.0, 0.05 * (ymax if ymax > 0 else 1.0)))

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, loc="upper center", frameon=False, bbox_to_anchor=(0.5, 1.18), ncol=3, fontsize=10)

    fig.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)
    try:
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        if os.path.exists(out_path):
            print(f"Saved stressor-only plot to {os.path.abspath(out_path)}")
    except Exception as e:
        print(f"ERROR saving stressor-only plot to {os.path.abspath(out_path)}: {e}", file=sys.stderr)
        try:
            plt.close(fig)
        except Exception:
            pass
    return out_path if os.path.exists(out_path) else ""
def save_final_style_plot(
    output_dir: str,
    times: np.ndarray,
    simulated_values_original: np.ndarray,
    simulated_values: np.ndarray,
    crh_values: np.ndarray,
    length_plot: int,
    num_days: int,
    stressor_parameters: dict,
    base_clock: str = "09:00:00",
    filename: str = "final_style.png",
    show_legend: bool = True,
    crh_label: str = "Perturbed CRH",
    max_minutes: Optional[int] = 2300,
) -> str:
    """Save a figure styled like first_paper/figure6.py.

    - Use the last 3-day window of the simulation (length_plot minutes).
    - Plot only up to 2300 minutes for layout.
    - Use stressor_parameters['time_in_scope'] for the stressor marker.
    """

    def to_series(arr: np.ndarray) -> np.ndarray:
        arr = np.asarray(arr)
        if arr.ndim == 2 and arr.shape[1] >= 2:
            return arr[:, 1]
        return arr.squeeze()

    # Convert to 1D series
    signal_full = to_series(simulated_values_original)
    pert_full = to_series(simulated_values)
    crh_full = np.asarray(crh_values).squeeze()

    # The simulated series provided are already trimmed to the last `length_plot` minutes.
    # Use them directly, and align CRH by taking the last matching segment.
    signal = signal_full
    pert = pert_full
    if len(crh_full) >= len(signal):
        crh_series = crh_full[-len(signal):]
    else:
        # Fallback: if CRH is shorter (shouldn't happen), trim signals to CRH length
        crh_series = crh_full
        signal = signal[: len(crh_series)]
        pert = pert[: len(crh_series)]

    # Determine time points for the window
    window = int(length_plot)
    local_len = min(window, len(signal), len(pert), len(crh_series))
    if local_len <= 0:
        print("Warning: No data available for final-style plot window.")
        return ""

    time_points = np.arange(local_len, dtype=float)

    # Cap to a maximum window length if provided (default 2300 like the paper)
    if max_minutes is None or max_minutes <= 0:
        plot_upto = local_len
    else:
        plot_upto = min(int(max_minutes), local_len)
    time_points = time_points[:plot_upto]
    signal = signal[:plot_upto]
    pert = pert[:plot_upto]
    crh_series = crh_series[:plot_upto]

    # Build the figure
    plt.rcParams.update({"font.size": 11})
    fig = plt.figure(figsize=(8, 3))
    ax1 = plt.gca()
    ax1.plot(time_points, signal, label="Normal signal", color="orange", linewidth=2)
    ax1.plot(time_points, pert, label="Perturbed signal", color="red")

    ax2 = ax1.twinx()
    ax2.plot(time_points, crh_series, label=crh_label, color="grey")

    # Stressor timing and duration (like figure6.py):
    # - time_in_scope is already in the plotting window coordinates
    # - surgery duration defaults to rise + plateau, unless an explicit surgery duration is provided
    stress_time_local = float(stressor_parameters.get("time_in_scope", 0))
    rise_duration = float(stressor_parameters.get("duration", 0))
    plateau_duration = float(stressor_parameters.get("plateau", 0))
    surgery_override = stressor_parameters.get("surgery_duration")
    if surgery_override is None:
        surgery_override = stressor_parameters.get("surgery_length")
    if surgery_override is None:
        surgery_duration = rise_duration + plateau_duration
    else:
        try:
            surgery_duration = float(surgery_override)
        except Exception:
            surgery_duration = rise_duration + plateau_duration

    if 0 <= stress_time_local <= local_len:
        ax1.axvline(stress_time_local, color="red", linestyle="--", label="Time of surgery")
        span_start = max(0, stress_time_local)
        span_end = min(stress_time_local + surgery_duration, local_len)
        if span_end > span_start:
            ax1.axvspan(span_start, span_end, color="red", alpha=0.2, label="Surgery duration")

    # Labels and legend
    ax1.set_xlabel("Time (hours)")
    ax1.set_ylabel("CORT (nmol/L)")
    ax2.set_ylabel("Circadian drive (CRH)\n(Arbitrary units)")
    ax2.spines["right"].set_color("grey")
    ax2.tick_params(axis="y", colors="grey")
    ax2.yaxis.label.set_color("grey")
    if show_legend:
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(
            handles1 + handles2,
            labels1 + labels2,
            loc="upper center",
            ncol=5,
            frameon=False,
            bbox_to_anchor=(0.5, 1.22),
            fontsize=11,
            columnspacing=0.8,
        )

    # X ticks like the paper: every 6 hours starting from base_clock
    base_dt = datetime.strptime(base_clock, "%H:%M:%S")
    ticks = np.arange(0, (time_points[-1] if len(time_points) else 0) + 1, 360)
    labels = [(base_dt + timedelta(minutes=float(t))).strftime("%H:%M") for t in ticks]
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)

    # Y bounds with a little headroom
    if len(signal) and len(pert):
        ymax = max(np.nanmax(signal), np.nanmax(pert))
        ax1.set_ybound(0, float(ymax) + 50)

    fig.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)

    # Robust save with diagnostics
    try:
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        if os.path.exists(out_path):
            print(f"Saved final-style plot to {os.path.abspath(out_path)}")
        else:
            print(
                f"Attempted to save final-style plot to {os.path.abspath(out_path)} but the file was not found right after saving."
            )
            try:
                entries = os.listdir(output_dir)
                print("Directory listing for output_dir:")
                for name in entries:
                    print(f" - {name}")
            except Exception as e:
                print(f"Could not list directory {os.path.abspath(output_dir)}: {e}")
    except Exception as e:
        print(f"ERROR saving final-style plot to {os.path.abspath(out_path)}: {e}", file=sys.stderr)
        try:
            plt.close(fig)
        except Exception:
            pass
    return out_path if os.path.exists(out_path) else ""


#Change directory - only if needed
# os.chdir('ddemodel')  # Commented out - run from code/ directory

#Configurations
params_config = 'parameters.json'
acth_data_file = f'data/data_shifted_ACTH_1T_9_00.csv'
cort_data_file = f'data/data_shifted_Cortisol_1T_9_00.csv'

def parse_arguments():
    parser = argparse.ArgumentParser(description='Config file')
    parser.add_argument('--config', default=params_config, help='Path parameter config file')
    parser.add_argument('--scenario', default=None, help='HS scenario name: HS1, HS2, HS3, HS4, ACUTE')
    parser.add_argument('--magnitude', type=float, help='Parameter magnitude')
    parser.add_argument('--duration', type=float, help='Parameter duration (rise time in minutes)')
    parser.add_argument('--plateau', type=float, help='Plateau duration in minutes')
    parser.add_argument('--beta', type=float, help='Beta shape parameter')
    parser.add_argument('--decay_shape', type=float, help='Decay shape parameter (kappa)')
    parser.add_argument('--time_in_scope', type=float, help='Parameter time in scope')
    parser.add_argument('--save', action='store_true', help='Save results')
    parser.add_argument('--no_legend', action='store_true', help='Disable legend on saved figure')
    parser.add_argument('--crh_source', choices=['perturbed', 'original'], default='perturbed', help='Which CRH to plot on the right axis')
    parser.add_argument('--plot_days', type=int, help='Override length_of_plot_stressor (days) for this run')
    parser.add_argument('--max_minutes_plot', type=int, help='Override the 2300-minute plotting cap (use a larger number)')
    args = parser.parse_args()
    stressor_parameters = {
        'duration': args.duration, 
        'magnitude': args.magnitude, 
        'time_in_scope': args.time_in_scope,
        'plateau': args.plateau,
        'beta': args.beta,
        'decay_shape': args.decay_shape
    }
    return {'config_file': f'model/config/base/{args.config}'}, stressor_parameters, args.scenario, args.save, args.no_legend, args.crh_source, args.plot_days, args.max_minutes_plot

args, stressor_from_args, scenario_from_args, save, no_legend, crh_source, plot_days_override, max_minutes_plot_override = parse_arguments()

# Choose config file path - adjust for running from code/ directory

config_file = args['config_file']

print(f'Current working directory: {os.getcwd()}')
print(f'Looking for config file at: {config_file}')
print(f'Absolute path: {os.path.abspath(config_file)}')
print(f'File exists: {os.path.exists(config_file)}')

with open(config_file, 'r') as file:
    config = json.load(file)

print(f'Getting options from config file {config_file}.')

SCENARIOS = {
    "HS1": {"duration": 20, "plateau": 160, "beta": 1.5, "decay_shape": 0.01, "time_in_scope": 1440+60+45, "magnitude": 1000},
    "HS2": {"duration": 20, "plateau": 160, "beta": 1.5, "decay_shape": 0.01, "time_in_scope": 1440+45, "magnitude": 1000},
    "HS3": {"duration": 180, "plateau": 0, "beta": 1, "decay_shape": 0.01, "time_in_scope": 1440+60+45, "magnitude": 1000},
    "HS4": {"duration": 6*60,  "plateau": 0, "beta": 0.8, "decay_shape": 0, "surgery_duration": 180, "time_in_scope": 1440+60+45, "magnitude": 1000},
    "ACUTE": {"duration": 120, 'magnitude': 100, 'time_in_scope': 1440}  
}

# Get options from parameter file
parameters = config.get('parameters', {})
fixed_params = config.get('fixed_params', {})
signal = config.get('signal', 'both')
participant_number = config.get('participant', 9)
num_days = config.get('num_days', 4)
reject = config.get('reject', True)
show_observed_values = config.get('show_observed_values', True)
final_value = config.get('final_value', False)
plot_option = config.get('plot_option', 1)
plot_crh = config.get('plot_crh', False)

print(f'Getting stressor options')
length_of_plot_in_days = config.get('length_of_plot_stressor', 3)
if plot_days_override is not None and plot_days_override > 0:
    length_of_plot_in_days = int(plot_days_override)
length_plot = 1440 * length_of_plot_in_days
default_stressor = {
    "magnitude": 50,
    "time_in_scope": 160,
    "duration": 10
}
stressor_parameters_json = config.get('stressor', default_stressor)

# Merge: args > config > scenario > default, with default only filling missing keys
scenario = (scenario_from_args or config.get('stressor', {}).get('scenario', 'HS1')).upper()
profile = SCENARIOS.get(scenario, {})

# Start with command line args, then config, then scenario (profile)
stressor_parameters = {}
for d in (default_stressor, stressor_parameters_json, profile, stressor_from_args):
    for k, v in d.items():
        if v is not None:
            stressor_parameters[k] = v

# Now ensure any missing keys from default_stressor are filled in
for k, v in default_stressor.items():
    stressor_parameters.setdefault(k, v)

# Special handling for kappa/decay_shape
if 'kappa' in stressor_parameters:
    stressor_parameters['decay_shape'] = stressor_parameters['kappa']

stressor_parameters['start'] = stressor_parameters['time_in_scope'] + length_plot * (num_days - 1)
stressor_type = 'HS' if scenario != 'ACUTE' else 'ACUTE'

# After merging scenario/profile into stressor_parameters, ensure 'rho' is None for HS1 so it is auto-calculated:
if scenario == "HS1":
    stressor_parameters['rho'] = None  # force auto-compute based on duration

print(f"Stressor parameters just before model creation:")
print(f"  duration: {stressor_parameters.get('duration')}")
print(f"  plateau: {stressor_parameters.get('plateau')}")
print(f"  start: {stressor_parameters.get('start')}")
print(f"  magnitude: {stressor_parameters.get('magnitude')}")
print(f"  rho: {stressor_parameters.get('rho')}")
print(f"  time_in_scope: {stressor_parameters.get('time_in_scope')}")
print(f"  All stressor_parameters: {stressor_parameters}")
print(f'Setting up model')

dde_model = Model(
    fixed_params=fixed_params,
    suggested_params=parameters,
    signal=signal,
    num_days=num_days,
    reject=reject,
    stressor_parameters=stressor_parameters,
    length_model=length_plot,
    stressor_type=stressor_type
)

timesteps = length_plot * num_days
times = np.linspace(0, timesteps, timesteps)
print(f"  times[0]: {times[0]}, times[-1]: {times[-1]}, len(times): {len(times)}")
parameters = dde_model.suggested_parameters()

print(f"Simulating model values.")
simulated_values = dde_model.simulate(parameters, times)

print(f'Calculating CRH values for stressor signal')
# CRH under stressor (perturbed)
crh_values_full_stressor = [dde_model.crh(t) for t in times]
crh_values = crh_values_full_stressor[length_plot * (num_days - 1):]

# Compute stressor-only contribution by subtracting baseline (without stressor)
print('Calculating stressor-only CRH contribution')
# Temporarily disable stressor to get baseline CRH
stressor_obj_backup = dde_model.stressor_object
dde_model.stressor_object = None
crh_values_full_baseline = [dde_model.crh(t) for t in times]
# Restore stressor object
dde_model.stressor_object = stressor_obj_backup
stressor_only_full = np.asarray(crh_values_full_stressor, dtype=float) - np.asarray(crh_values_full_baseline, dtype=float)

dde_model.stressor_object = None

print(f'Calculating CRH values for original signal')
simulated_values_original = dde_model.simulate(parameters, times)
crh_values_original_full = [dde_model.crh(t) for t in times]
crh_values_original = crh_values_original_full[length_plot * (num_days - 1):]

fig, axes = model_plotting_functions.pertubation_plot(
    times,
    simulated_values,
    num_days,
    config,
    crh_values,
    simulated_values_original=simulated_values_original,
    crh_values_original=crh_values_original,
    length_model=length_plot,
    yaxis=None,
    time_of_stressor=stressor_parameters['time_in_scope']
)

plt.show()
metrics, acth_peaks, acth_peaks_shifted, cort_peaks, cort_peaks_shifted = mc.run_metrics(
    stressor_parameters,
    crh_values_original,
    crh_values,
    simulated_values_original,
    simulated_values
)

print(f'Save is {save}')

if save:
    # Compose output directory once
    output_dir = build_output_dir(scenario, stressor_parameters)

    metrics['duration'] = stressor_parameters['duration']
    metrics['time_in_scope'] = stressor_parameters['time_in_scope']
    metrics['magnitude'] = stressor_parameters['magnitude']
    metrics['start'] = stressor_parameters['start']
    metrics['length_of_plot'] = length_plot

    print(metrics)

    os.makedirs(output_dir, exist_ok=True)

    # 1) Save all arrays
    save_numpy_arrays(output_dir, {
        "simulated_values_original": simulated_values_original,
        "simulated_values_stressor": simulated_values,
        "acth_peaks": acth_peaks,
        "acth_peaks_stressor": acth_peaks_shifted,
        "cort_peaks": cort_peaks,
        "crh_stressor": crh_values,
        "crh": crh_values_original,
        "cort_peaks_shifted": cort_peaks_shifted,
        "stressor_only_full": stressor_only_full,
    })

    # 2) Save metrics
    save_metrics_csv(metrics, output_dir)

    # 3) Save the default simulation plot
    fig.savefig(os.path.join(output_dir, 'simulation.png'))

    # 4) Save the final-style plot
    print("Saving final-style plot...")
    # Choose which CRH series to plot on the right axis
    if crh_source == 'original':
        crh_series_for_plot = crh_values_original_full
        crh_label = 'Original CRH'
    else:
        crh_series_for_plot = crh_values_full_stressor
        crh_label = 'Perturbed CRH'
    # Save the default (paper-style) plot, always capped at 2300 for manuscript consistency
    final_plot_path = save_final_style_plot(
        output_dir,
        times,
        simulated_values_original,
        simulated_values,
        crh_series_for_plot,
        length_plot,
        num_days,
        stressor_parameters,
        crh_label=crh_label,
        show_legend=not no_legend
    )
    if not final_plot_path:
        print(
            f"Warning: Final-style plot was not created. Expected location was {os.path.abspath(os.path.join(output_dir, 'final_style.png'))}"
        )

    # If user requested a longer plot, save a second version with the full window (or up to max_minutes_plot_override)
    if max_minutes_plot_override is not None and max_minutes_plot_override > 2300:
        # Use the full window if max_minutes_plot_override exceeds the available data
        full_window = length_plot
        max_minutes_long = max(max_minutes_plot_override, full_window)
        long_plot_path = save_final_style_plot(
            output_dir,
            times,
            simulated_values_original,
            simulated_values,
            crh_series_for_plot,
            length_plot,
            num_days,
            stressor_parameters,
            crh_label=crh_label,
            show_legend=not no_legend,
            max_minutes=full_window,
            filename="final_style_long.png",
        )
        if not long_plot_path:
            print(
                f"Warning: Long final-style plot was not created. Expected location was {os.path.abspath(os.path.join(output_dir, 'final_style_long.png'))}"
            )

    # 5) Save the stressor-only simple plot
    print("Saving stressor-only plot...")
    stressor_only_path = save_stressor_only_plot(
        output_dir,
        stressor_only_full,
        length_plot,
        num_days,
        stressor_parameters,
        filename="stressor_only.png",
        max_minutes=max_minutes_plot_override if max_minutes_plot_override is not None else 2300,
    )
    if not stressor_only_path:
        print(
            f"Warning: Stressor-only plot was not created. Expected location was {os.path.abspath(os.path.join(output_dir, 'stressor_only.png'))}"
        )
