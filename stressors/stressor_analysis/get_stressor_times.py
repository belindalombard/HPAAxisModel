from scipy.signal import hilbert
import stressor_analysis as sa
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.signal import find_peaks
import argparse


def find_peaks_and_throughs(signal): 
    peaks, _ = find_peaks(signal, height=None, distance=None)
    troughs, _ = find_peaks(-signal, height=None, distance=None)
    return peaks, troughs


def calculate_ultradian_period(signal, time_array):

    peaks, _ = find_peaks(signal, height=None, distance=None)
    
    peak_times = time_array[peaks]
    periods = np.diff(peak_times) 
    
    ultradian_period = np.mean(periods)
    return ultradian_period

# Calculate phase of signal ultradian rhythm
def calculate_phase(signal, point_time, time_ref): 
    phase_length = len(signal)
    phase_radians = ((point_time-time_ref)/phase_length)*2*np.pi
    return phase_radians

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Analyze stressor timing and phase')
parser.add_argument('--directory', type=str, 
                    default='../output/ACUTE/mag50.00_dur20.00_plat0.00_decay0.00_beta0.00_start14616.00/',
                    help='Path to simulation directory')
parser.add_argument('--scenario', type=str, default='ACUTE',
                    help='Scenario subdirectory (e.g., ACUTE, HS1, HS2, HS3, HS4)')
args = parser.parse_args()

# Obtaining signals
print(f'################')
print(f'Obtaining signals')

directory = args.directory
signal, pert = sa.get_signals_from_directory(directory, 'Cortisol')

print(f'Signal is {signal}')
#times = np.linspace(0, 7200, 7200)



print(f'################')
print(f'Find peaks and troughs')
print(f'Signal length {len(signal)}')
signal_chopped = signal[0:1440*2]
peaks, troughs = find_peaks_and_throughs(signal_chopped)
print(f'Peaks: {peaks+1440}')
print(f'Troughs {troughs+1440}')


times_chopped = np.arange(0, 2880)

plt.plot(times_chopped, signal_chopped, 'b-')
for peak in peaks: 
    plt.plot(times_chopped[peak], signal_chopped[peak], 'ro')

for trough in troughs: 
    plt.plot(times_chopped[trough], signal_chopped[trough], 'go')


plt.grid()
plt.show()

#throughs we're interested in: 
selected_throughs = [1, 4, 6, 7] 

ultradian_period = calculate_ultradian_period(signal_chopped, times_chopped)
print(f'Ultradian period is: {ultradian_period}')

# Get CRH values
crh_values = np.load(f'{directory}/crh.npy')

through_index = np.argmin(crh_values)

plt.show()

fig, ax1 = plt.subplots(figsize=(8, 4), constrained_layout=True)
ax1.plot(times_chopped, signal_chopped, 'r-', label='Cortisol')

points_dictionary = []

for i in range(len(selected_throughs)): 
    # Selected through is i
    # Extract through position from throughs array
    through_idx = selected_throughs[i]
    if through_idx >= len(troughs):
        # skip if requested index doesn't exist
        continue
    through = troughs[through_idx]
    # Use next trough if it exists, otherwise fall back to the end of the chopped window
    if through_idx + 1 < len(troughs):
        follow_through = troughs[through_idx + 1]
    else:
        follow_through = len(signal_chopped) - 1
    # guard in case indices are equal or reversed
    if follow_through <= through:
        follow_through = min(len(signal_chopped) - 1, through + 1)
    ultradian_cycle = signal_chopped[through:follow_through]


    points_dictionary.append({'ultradian_points': []})

    # Find corresponding peak (first peak after the selected trough)
    find_peak = None
    for peak in peaks:
        if peak > through:
            find_peak = peak
            break
    if find_peak is None:
        # no peak after this trough; skip this ultradian segment
        continue

    midpoint = ((find_peak - through)//2) + through
    midpoint_falling = (((follow_through - find_peak))//2)+find_peak

    points_dictionary[-1]['ultradian_points'].append({'timing': int(find_peak), 'phase': float(np.pi), 'describe': 'peak'})
    ax1.plot([find_peak, find_peak], [0, signal[find_peak]], color='g', linestyle='--')

    points_dictionary[-1]['ultradian_points'].append({'timing': int(midpoint), 'phase': float(np.pi/2), 'describe': 'midpoint_rising'})
    ax1.plot([midpoint, midpoint], [0, signal[midpoint]], color='g', linestyle='--')

    points_dictionary[-1]['ultradian_points'].append({'timing': int(through), 'phase': 0, 'describe': 'through'})
    ax1.plot([through, through], [0, signal[through]], color='g', linestyle='--')


    points_dictionary[-1]['ultradian_points'].append({'timing': int(midpoint_falling), 'phase': float(3*np.pi/2), 'describe': 'midpoint_falling'})
    ax1.plot([midpoint_falling, midpoint_falling], [0, signal[midpoint_falling]], color='g', linestyle='--')
    
    #ax1.plot([follow_through, follow_through], [0, signal[follow_through]], color='r', linestyle='--')

    #for point in points_dictionary[-1]['ultradian_points']: 
    #    ax1.plot(point['timing'], signal[point['timing']], 'ro')
        #radians = calculate_phase(ultradian_cycle, point, through)
    #    ax1.text(point['timing'], signal[point['timing']]+10, str(round(point['phase'], 2)), ha='center')    

ax1.set_ylim(0, None)

with open("points_to_consider.json", "w") as outfile: 
    json.dump(points_dictionary, outfile)

    
ax2 = ax1.twinx()
ax2.plot(times_chopped, crh_values[0:2880], label='CRH Drive', color='black')
ax2.set_ylabel('CRH Drive')
ax1.set_xlabel('Time (minutes)')
ax1.set_ylabel('Cortisol')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
fig.legend(lines1 + lines2, labels1 + labels2, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=2, fontsize='small', frameon=False)

plt.savefig('../output/stressor_fig.pdf')
plt.savefig('../output/stressor_fig.png')

plt.show()
 



'''
def extract_phases(signal, time_array, ultradian_period, circadian_period=1440):
    ultradian_phase = (2 * np.pi * (time_array % ultradian_period) / ultradian_period) % (2 * np.pi)
    circadian_phase = (2 * np.pi * (time_array % circadian_period) / circadian_period) % (2 * np.pi)
    return ultradian_phase, circadian_phase







ultradian_period = calculate_ultradian_period(signal, times)

ultradian_phases, circadian_phases = extract_phases(signal, times, ultradian_period)


def extract_ultradian_cycle(signal, times, start_time, ultradian_period):
    start_idx = np.abs(times - start_time).argmin() 
    end_time = start_time + ultradian_period
    extracted_indices = (times >= start_time) & (times <= end_time)
    extracted_signal = signal[extracted_indices]
    extracted_times = times[extracted_indices]
    
    return extracted_times, extracted_signal


# Cut off 24hr signal

'''