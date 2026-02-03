import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from scipy.signal import find_peaks
import numpy as np
# filesystem
import os
#from sklearn import metrics
import json
from datetime import datetime, timedelta
from matplotlib.patches import FancyArrowPatch


def get_stressor_directories(magnitude=None, duration=None, time_in_scope=None, directory='../output/ACUTE'):
    """Return list of directories matching the requested magnitude/duration/time_in_scope.

    This is robust to two naming schemes: either simple "<mag>_<dur>_<tim>" or key/value tokens
    like "mag50.00_dur0.00_..._start14556.00". It lists all subdirectories and filters by
    parsing their basenames using parse_directory_name when possible.
    """
    # get all candidate subdirectories
    candidates = glob.glob(os.path.join(directory, '*/'))

    # if no filters requested, return all
    if magnitude is None and duration is None and time_in_scope is None:
        return candidates

    filtered = []
    for d in candidates:
        base = os.path.basename(os.path.normpath(d))
        try:
            mag_i, dur_i, tim_i = parse_directory_name(base)
        except Exception:
            # couldn't parse this directory name, skip it
            continue

        ok = True
        if magnitude is not None:
            if isinstance(magnitude, (list, tuple, np.ndarray)):
                ok = ok and (mag_i in magnitude)
            else:
                ok = ok and (int(mag_i) == int(magnitude))
        if duration is not None:
            if isinstance(duration, (list, tuple, np.ndarray)):
                ok = ok and (dur_i in duration)
            else:
                ok = ok and (int(dur_i) == int(duration))
        if time_in_scope is not None:
            if isinstance(time_in_scope, (list, tuple, np.ndarray)):
                ok = ok and (tim_i in time_in_scope)
            else:
                ok = ok and (int(tim_i) == int(time_in_scope))

        if ok:
            filtered.append(d)

    return filtered


def filter_stressor_data(magnitude=None, duration=None, time_in_scope=None, filename='ddemodel/stressor_output/metrics.csv'): 
    import pandas as pd
    
    data = pd.read_csv(filename)
    
    if magnitude is not None: 
        if isinstance(magnitude, (list, tuple, np.ndarray)): 
            data = data[data['magnitude'].isin(magnitude)]
        else: 
            data = data[data['magnitude'].astype(int) == int(magnitude)]
    
    if duration is not None: 
        if isinstance(duration, (list, tuple, np.ndarray)):
            data = data[data['duration'].isin(duration)]
        else:
            data = data[data['duration'] == duration]
    
    if time_in_scope is not None: 
        if isinstance(time_in_scope, (list, tuple, np.ndarray)):
            data = data[data['time_in_scope'].isin(time_in_scope)]
        else:
            data = data[data['time_in_scope'] == time_in_scope]
    
    return data


def open_image(magnitude, duration, time_in_scope): 
    directory = get_stressor_directories(magnitude, duration, time_in_scope)
    if len(directory) > 1: 
        print(f'Please provide all three parameters')

    if len(directory) == 0: 
        print(f'Please provide valid parameters')
        return 
    image = mpimg.imread(f'{directory[0]}/simulation.png')

    plt.imshow(image)
    plt.axis('off')  
    plt.show()



def get_signals(magnitude, duration, time_in_scope, signal = 'Cortisol'): 
    directory = get_stressor_directories(magnitude, duration, time_in_scope)
    if len(directory) > 1: 
        print(f'Please provide all three parameters')

    if len(directory) == 0: 
        print(f'Please provide valid parameters')
        return 
    #print(directory[0])
    original_signal = np.load(f'{directory[0]}/simulated_values_original.npy')
    pertubated_signal = np.load(f'{directory[0]}/simulated_values_stressor.npy')
    
    if signal == 'Cortisol': 
        num = 1
    else: 
        num = 0

    return original_signal[:, num], pertubated_signal[:, num]


def get_signals_from_directory(directory, signal = 'Cortisol'): 
    
    original_signal = np.load(f'{directory}/simulated_values_original.npy')
    pertubated_signal = np.load(f'{directory}/simulated_values_stressor.npy')
    
    if signal == 'Cortisol': 
        num = 1
    else: 
        num = 0

    return original_signal[:, num], pertubated_signal[:, num]


def parse_directory_name(dir_name):
    """Parse a directory base name and return (magnitude, duration, timing) as ints.

    Supports two formats:
    - simple underscore-separated: "<mag>_<dur>_<tim>" (e.g. "50.00_0.00_14556.00" or "50_0_14556")
    - key-value tokens: "mag50.00_dur0.00_..._start14556.00"

    If a "start" token is present, timing is the absolute index (no modulo). Downstream code can
    derive minute-of-day via timing % 1440 where needed.
    """
    import re
    s = dir_name.lower()
    parts = dir_name.split('_')

    # Case 1: exactly three parts and they look numeric
    if len(parts) == 3:
        try:
            mag = float(parts[0])
            dur = float(parts[1])
            tim = float(parts[2])
            return int(mag), int(dur), int(tim)
        except Exception:
            # fall through to more flexible parsing
            pass

    mag = None
    dur = None
    tim = None

    for token in parts:
        m = re.search(r'[-+]?\d*\.?\d+', token)
        if not m:
            continue
        num = float(m.group(0))
        t = token.lower()
        if t.startswith('mag') and mag is None:
            mag = num
            continue
        if t.startswith('dur') and dur is None:
            dur = num
            continue
        if t.startswith('start') or t.startswith('time') or t.startswith('tim'):
            tim = num
            continue
        # fallback: assign first unknown numeric to mag, then dur, then tim
        if mag is None:
            mag = num
        elif dur is None:
            dur = num
        elif tim is None:
            tim = num

    if mag is None or dur is None or tim is None:
        raise ValueError(f"Could not parse magnitude/duration/timing from '{dir_name}'")

    # If explicit 'start' token was provided keep absolute timing (do not modulo here).

    return int(mag), int(dur), int(tim)

# Calculates the difference in AUC after perturbation for CORT
def calculate_auc_after_perturbation(signal, pert, time_of_stressor, length, fig=None):
    # Ensure we don't index past the end of the arrays.
    start = int(time_of_stressor)
    end = min(len(signal), len(pert), start + int(length))

    sig_slice = signal[start:end]
    pert_slice = pert[start:end]

    # Compute AUCs over the same actual available window
    # Use trapezoid (preferred) and ensure same x-spacing so simple trapezoid is fine
    auc_signal = np.trapezoid(sig_slice)
    auc_pert = np.trapezoid(pert_slice)

    if fig:
        ax1 = fig.axes[0]
        time_range = np.arange(start, end)

        # Only plot if there is data
        if time_range.size > 0 and pert_slice.size > 0:
            ax1.fill_between(time_range, 0, pert_slice, color='red', alpha=0.3, label='AUC Pert')

        if time_range.size > 0 and sig_slice.size > 0:
            ax1.fill_between(time_range, 0, sig_slice, color='orange', alpha=0.5, label='AUC Signal')

        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4)

    return auc_pert - auc_signal


def calculate_peak_diff_after_perturbation(signal, pert, time_of_stressor, fig):
    # Ensure we start within array bounds
    start = int(time_of_stressor)
    if start >= len(signal) or start >= len(pert):
        return np.nan

    sig_seg = signal[start:]
    pert_seg = pert[start:]
    peaks_signal, _ = find_peaks(sig_seg)
    peaks_pert, _ = find_peaks(pert_seg)

    # If no peaks found in either, return NaN and skip plotting
    if len(peaks_signal) == 0 or len(peaks_pert) == 0:
        return np.nan

    peak_signal_time = start + peaks_signal[0]
    peak_pert_time = start + peaks_pert[0]
    peak_diff = float(pert[peak_pert_time] - signal[peak_signal_time])

    if fig:
        ax1 = fig.axes[0]

        ax1.plot(peak_signal_time, signal[peak_signal_time], 'o', label='Signal peak', color='orange', markersize=8, markeredgecolor='black', markeredgewidth=0.5)
        ax1.plot(peak_pert_time, pert[peak_pert_time], 'o', label='Perturbation peak', color='red', markersize=8, markeredgecolor='black', markeredgewidth=0.5)

        ax1.annotate("",
                     (peak_signal_time, signal[peak_signal_time]),
                     textcoords="offset points", xytext=(-15, 10), ha='center', color='orange')

        ax1.annotate("",
                     (peak_pert_time, pert[peak_pert_time]),
                     textcoords="offset points", xytext=(-15, 10), ha='center', color='red')

        ax1.vlines(x=peak_signal_time,
                   ymin=signal[peak_signal_time],
                   ymax=pert[peak_pert_time],
                   color="black", linestyle="-", linewidth=1.5, label="Difference in peaks")

        ax1.hlines(y=signal[peak_signal_time],
                   xmin=peak_signal_time - 20,
                   xmax=peak_signal_time + 20,
                   color="black", linestyle="-", linewidth=1.5)

        ax1.hlines(y=pert[peak_pert_time],
                   xmin=peak_signal_time - 20,
                   xmax=peak_signal_time + 20,
                   color="black", linestyle="-", linewidth=1.5)

        mid_y = (signal[peak_signal_time] + pert[peak_pert_time]) / 2
        ax1.text(peak_signal_time + 50, mid_y, f"Δ = {peak_diff:.2f}",
                 fontsize=10, ha='left', color='black')

        ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4)

    return peak_diff

def plot_time_series(normal_signal, perturbed_signals, time_points, crh_signals = [], stress_times=None, save_dir = None):
    fig, ax1 = plt.subplots(figsize=(10, 5), constrained_layout=True)
    # Plot the last two days (2880 minutes) if available
    DAY = 1440
    WINDOW = 2 * DAY
    orig_len = len(normal_signal)
    # slice normal signal to last 2 days if possible and record start index in original signal
    if orig_len >= WINDOW:
        start_index = orig_len - WINDOW
        normal_signal = normal_signal[-WINDOW:]
    else:
        start_index = 0

    # slice perturbed signals to last 2 days where possible
    new_pert = []
    for ps in perturbed_signals:
        if len(ps) >= WINDOW:
            new_pert.append(ps[-WINDOW:])
        else:
            new_pert.append(ps)
    perturbed_signals = new_pert

    # slice crh signals similarly
    new_crh = []
    for cs in crh_signals:
        if len(cs) >= WINDOW:
            new_crh.append(cs[-WINDOW:])
        else:
            new_crh.append(cs)
    crh_signals = new_crh

    # adjust time_points to match the plotted signal length
    plot_len = len(normal_signal)
    # always use a simple 0..(plot_len-1) minute axis
    time_points = np.arange(0, plot_len)

    ax1.plot(time_points, normal_signal, label="Normal signal", color="orange", linewidth=2, linestyle="-")
    ax1.plot(time_points, perturbed_signals[0], label="Perturbed signal", color="red", linestyle="-")

    line_handles = [ax1.plot([], [], label="Normal signal", color="orange")[0],
                    ax1.plot([], [], label="Perturbed signal", color="red")[0]]

    ax2 = ax1.twinx()
    for i, crh_signal in enumerate(crh_signals):
        ax2.plot(time_points, crh_signal, linestyle="-", label=f"CRH", color="grey")

    crh_handle = ax2.plot([], [], label="Circadian drive (CRH)", color="grey")[0]

    if stress_times is not None:
        # convert stress times to plotted-axis coordinates (account for slicing)
        for stress_time in stress_times:
            plot_x = stress_time - start_index
            if 0 <= plot_x < len(normal_signal):
                ax1.axvline(plot_x, color='red', label='Time of stressor', linestyle='--')

    # Determine time_of_stressor relative to the (possibly sliced) arrays for calculations
    timing_in_plot = None
    if stress_times is not None and len(stress_times) > 0:
        # Normalize absolute stress time to the 2-day plotted window
        abs_time = int(stress_times[0])
        timing_in_plot = abs_time - start_index
        # If outside the 2-day slice, wrap into [0, WINDOW)
        if timing_in_plot < 0 or timing_in_plot >= WINDOW:
            timing_in_plot = abs_time % WINDOW

    if timing_in_plot is not None:
        # Compute AUC over one day from the stressor point within the 2-day plot window
        calculate_auc_after_perturbation(normal_signal, perturbed_signals[0], timing_in_plot,  DAY, fig)
        calculate_peak_diff_after_perturbation(normal_signal, perturbed_signals[0], timing_in_plot, fig)


    auc_handle_signal = ax1.fill_between([], [], [], color='orange', alpha=0.8, label="AUC Signal")
    auc_handle_pert = ax1.fill_between([], [], [], color='red', alpha=0.3, label="AUC Pert")

    peak_handle_signal = ax1.plot([], [], 'o', label='Signal peak', color='orange')[0]
    peak_handle_pert = ax1.plot([], [], 'o', label='Perturbation peak', color='red')[0]
    diff_handle = ax1.plot([], [], 'k-', label="Difference in peaks")[0]

    handles = [line_handles[0], line_handles[1], auc_handle_signal, auc_handle_pert, peak_handle_signal, peak_handle_pert, diff_handle, crh_handle]
    labels = ["Normal signal", "Perturbed signal", "AUC Signal", "AUC Pert", "Signal peak", "Perturbation peak", "Difference in peaks", "Circadian drive (CRH)"]

    ax1.legend(handles=handles, labels=labels, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4)
    ax1.set_xlabel("Time (hours)", fontsize=12)
    ax1.set_ylabel("Cortisol (nmol/L)", fontsize=12)
    ax2.set_ylabel("Circadian drive (CRH)", fontsize=12)
    ax2.spines['right'].set_color('grey')    
    ax2.tick_params(axis='y', colors='grey')
    ax2.yaxis.label.set_color('grey')
    plt.grid()

    start_time='09:00:00'
    start_time = datetime.strptime(start_time, "%H:%M:%S")
    # set xticks dynamically based on plotted length
    max_minutes = len(time_points)
    step = 360  # every 6 hours
    ax1.set_xticks([i for i in range(0, max_minutes + 1, step)])
    ax1.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(0, max_minutes+1, step)])


    if save_dir: 
        plt.savefig(save_dir)
    else: 
        return fig



all_directories = get_stressor_directories()
times = np.linspace(0, 1440*3, 1440*3)
'''
for i in range(len(all_directories)): 
    print(f'Now doing {all_directories[i]}')
    base = all_directories[i].split('\\')[-2]
    magnitude, duration, timing_abs = parse_directory_name(base)

    signal, pert = get_signals_from_directory(all_directories[i])

    crh_signal = np.load(f'{all_directories[i]}/crh_stressor.npy')

    # Build 2-day slices and normalize timing to that window for calculations
    DAY = 1440
    WINDOW = 2 * DAY
    start_index_full = max(0, len(signal) - WINDOW)
    signal_2d = signal[start_index_full:]
    pert_2d = pert[start_index_full:]
    crh_2d = crh_signal[start_index_full:]
    timing_in_plot = timing_abs - start_index_full
    if timing_in_plot < 0 or timing_in_plot >= WINDOW:
        timing_in_plot = timing_abs % WINDOW

    fig = plot_time_series(signal, [pert], times, stress_times = [timing_abs], crh_signals=[crh_signal], save_dir = f'{all_directories[i]}/phase_response.png' )
    metrics = {}
    metrics['difference_in_auc'] = float(calculate_auc_after_perturbation(signal_2d, pert_2d, timing_in_plot,  DAY, fig))
    metrics['peak_diff'] = float(calculate_peak_diff_after_perturbation(signal_2d, pert_2d, timing_in_plot, fig))

    # Add contextual phase labels to metrics.json for tracking
    try:
        timing_mod = int(timing_abs % 1440)
        with open('ddemodel/code/stressor_analysis/points_to_consider.json', "r") as file:
            points_data = json.load(file)
        ultradian_desc = None
        ultradian_phase_rad = None
        for entry in points_data:
            for point in entry.get("ultradian_points", []):
                if int(point.get("timing", -1)) == timing_mod:
                    ultradian_desc = point.get('describe')
                    ultradian_phase_rad = point.get('phase')
                    if ultradian_desc is None: 
                        ultradian_desc = 'midpoint_falling' #terrible, hardcoding.      
                    break
            if ultradian_desc:
                break

        # Hardcode: if this specific stressor time, force midpoint_falling
        if int(timing_abs) == 1440*11 + 4 or timing_mod == 4:
            ultradian_desc = 'midpoint_falling'
            ultradian_phase_rad = float(np.pi)

        # Circadian label from timing ranges (minute-of-day)
        if 50 < timing_mod < 500:
            circadian_desc = 'falling'
        elif 500 < timing_mod < 1000:
            circadian_desc = 'nadir'
        elif 980 < timing_mod < 1250:
            circadian_desc = 'rising'
        else:
            circadian_desc = 'peak'

        if ultradian_desc is not None:
            metrics['ultradian_phase_describe'] = ultradian_desc
        else:
            metrics['ultradian_phase_describe'] = 'midpoint_falling'
        if ultradian_phase_rad is not None:
            metrics['ultradian_phase_radians'] = ultradian_phase_rad

        metrics['circadian_phase_describe'] = circadian_desc
    except Exception:
        # if lookup fails, skip adding labels
        pass
    
    with open(f'{all_directories[i]}/metrics.json', "w") as f: 
        json.dump(metrics, f)

    plt.savefig(f'{all_directories[i]}/phase_response.pdf')
    plt.tight_layout()
    plt.close(fig)


#combined_metrics = pd.DataFrame(columns=["magnitude"
#                        , "duration"
#                        , "time_in_scope"
#                        , "ultradian_phase_radians"
#                        , "ultradian_phase_describe"
#                        , "circadian_phase_describe", "difference_in_auc", "peak_diif"])

'''
combined_metrics = []

for i in range(len(all_directories)): 
    base = all_directories[i].split('\\')[-2]
    magnitude, duration, timing_abs = parse_directory_name(base)
    with open('ddemodel/code/stressor_analysis/points_to_consider.json', "r") as file:
        data = json.load(file)
    
    phase = None
    point = None
    for entry in data:
        for point in entry["ultradian_points"]:
            # JSON timing stored as minute-of-day; compare with timing modulo
            if int(point["timing"]) == int(timing_abs % 1440):
                phase = point['phase']
                describe = point['describe']
    # Hardcode: if this specific stressor time, force midpoint_falling
    if int(timing_abs) == 1440*11 + 4 or int(timing_abs % 1440) == 4:
        describe = 'midpoint_falling'
        phase = float(3*np.pi/2)

    with open(f'{all_directories[i]}/metrics.json', "r") as file:
        data = json.load(file)
    
    difference_in_auc = data['difference_in_auc']
    peak_diff = data['peak_diff']

    timing_mod = int(timing_abs % 1440)
    if 50 < timing_mod < 500: 
        circadian_phase = 'falling'
    elif 500 < timing_mod < 1000: 
        circadian_phase = 'nadir'
    elif 980 < timing_mod < 1250: 
        circadian_phase = 'rising'
    else: 
        circadian_phase = 'peak'


    append_data = {
            "magnitude": magnitude,
            "duration": duration, 
            "time_in_scope": timing_mod, 
            "ultradian_phase_radians": phase,
            "ultradian_phase_describe": describe,
            "circadian_phase_describe": circadian_phase,
            "difference_in_auc": difference_in_auc,
            "peak_diff": peak_diff
    }
    combined_metrics.append(append_data)

metrics_df = pd.DataFrame(combined_metrics)

print(metrics_df)
import os
results_path = 'ddemodel/stressor_output/results.csv'
if os.path.exists(results_path):
    try:
        existing = pd.read_csv(results_path)
        # Combine preserving existing rows (don't redo entries). Use magnitude,duration,time_in_scope as key.
        combined_all = pd.concat([existing, metrics_df], ignore_index=True)
        combined_all = combined_all.drop_duplicates(subset=['magnitude', 'duration', 'time_in_scope'], keep='first')
        combined_all.to_csv(results_path, index=False)
    except Exception:
        # If anything goes wrong reading existing file, fall back to writing new file
        metrics_df.to_csv(results_path, index=False)
else:
    metrics_df.to_csv(results_path, index=False)
