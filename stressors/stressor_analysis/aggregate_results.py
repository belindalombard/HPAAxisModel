"""Aggregate ACUTE sweep results into results.csv for heatmap generation.

Reads each per-run directory in output/ACUTE/, computes difference_in_auc and
peak_diff from the saved .npy signals, adds circadian/ultradian phase labels
from points_to_consider.json, and writes ../output/ACUTE/results.csv.

Usage:
    cd stressors/stressor_analysis
    python aggregate_results.py

    # Or specify a different ACUTE directory:
    python aggregate_results.py --acute-dir ../output/ACUTE
"""

import os
import sys
import json
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Aggregate ACUTE sweep results')
parser.add_argument('--acute-dir', type=str,
                    default=os.path.join(_SCRIPT_DIR, '../output/ACUTE'),
                    help='Path to ACUTE output directory')
args = parser.parse_args()

acute_dir = os.path.abspath(args.acute_dir)
results_path = os.path.join(acute_dir, 'results.csv')
points_path = os.path.join(_SCRIPT_DIR, 'points_to_consider.json')

# ---------------------------------------------------------------------------
# Load points_to_consider.json for phase labels
# ---------------------------------------------------------------------------
with open(points_path) as f:
    points_data = json.load(f)

# Build a flat lookup: timing_mod → (ultradian_describe, ultradian_phase_rad)
phase_lookup = {}
for entry in points_data:
    for pt in entry.get('ultradian_points', []):
        phase_lookup[int(pt['timing'])] = (pt['describe'], pt['phase'])


def get_circadian_desc(timing_mod):
    if 50 < timing_mod < 500:
        return 'falling'
    elif 500 < timing_mod < 1000:
        return 'nadir'
    elif 980 < timing_mod < 1250:
        return 'rising'
    else:
        return 'peak'


# ---------------------------------------------------------------------------
# Parse directory name: mag50.00_dur20.00_plat0.00_decay0.00_beta0.00_start14616.00
# ---------------------------------------------------------------------------
import re

def parse_dir_name(name):
    nums = {}
    for token in name.split('_'):
        m = re.match(r'([a-z]+)([\d.]+)', token)
        if m:
            nums[m.group(1)] = float(m.group(2))
    mag = int(nums.get('mag', 0))
    dur = int(nums.get('dur', 0))
    start = int(nums.get('start', 0))
    return mag, dur, start


# ---------------------------------------------------------------------------
# Metric calculations (same logic as stressor_analysis.py)
# ---------------------------------------------------------------------------
DAY = 1440

def calc_difference_in_auc(signal, pert, time_of_stressor):
    start = int(time_of_stressor)
    end = min(len(signal), len(pert), start + DAY)
    return float(np.trapezoid(pert[start:end]) - np.trapezoid(signal[start:end]))


def calc_peak_diff(signal, pert, time_of_stressor):
    start = int(time_of_stressor)
    if start >= len(signal) or start >= len(pert):
        return np.nan
    peaks_sig, _ = find_peaks(signal[start:])
    peaks_pert, _ = find_peaks(pert[start:])
    if len(peaks_sig) == 0 or len(peaks_pert) == 0:
        return np.nan
    return float(pert[start + peaks_pert[0]] - signal[start + peaks_sig[0]])


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
candidates = sorted(glob.glob(os.path.join(acute_dir, '*/')))
print(f'Found {len(candidates)} directories in {acute_dir}')

rows = []
errors = []

for d in candidates:
    base = os.path.basename(os.path.normpath(d))
    try:
        mag, dur, start = parse_dir_name(base)
    except Exception:
        continue  # skip non-run dirs

    npy_orig = os.path.join(d, 'simulated_values_original.npy')
    npy_pert = os.path.join(d, 'simulated_values_stressor.npy')

    if not os.path.exists(npy_orig) or not os.path.exists(npy_pert):
        errors.append(f'Missing npy: {base}')
        continue

    try:
        orig = np.load(npy_orig)
        pert = np.load(npy_pert)
    except Exception as e:
        errors.append(f'Load error {base}: {e}')
        continue

    cort_orig = orig[:, 1]
    cort_pert = pert[:, 1]

    # Timing relative to the 2-day plot window (same as stressor_analysis.py)
    WINDOW = 2 * DAY
    start_index = max(0, len(cort_orig) - WINDOW)
    cort_2d = cort_orig[start_index:]
    pert_2d = cort_pert[start_index:]
    timing_in_plot = start - start_index
    if timing_in_plot < 0 or timing_in_plot >= WINDOW:
        timing_in_plot = start % WINDOW

    diff_auc = calc_difference_in_auc(cort_2d, pert_2d, timing_in_plot)
    peak_diff = calc_peak_diff(cort_2d, pert_2d, timing_in_plot)

    timing_mod = start % DAY
    ultradian_desc, ultradian_rad = phase_lookup.get(timing_mod, (None, None))
    circadian_desc = get_circadian_desc(timing_mod)

    rows.append({
        'magnitude': mag,
        'duration': dur,
        'time_in_scope': timing_mod,
        'ultradian_phase_radians': ultradian_rad,
        'ultradian_phase_describe': ultradian_desc,
        'circadian_phase_describe': circadian_desc,
        'difference_in_auc': diff_auc,
        'peak_diff': peak_diff,
    })

df = pd.DataFrame(rows)
df = df.drop_duplicates(subset=['magnitude', 'duration', 'time_in_scope'])
df.to_csv(results_path, index=False)

print(f'Written {len(df)} rows to {results_path}')
if errors:
    print(f'\n{len(errors)} errors:')
    for e in errors:
        print(f'  {e}')
