import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.signal import find_peaks
import seaborn as sns
from scipy.signal import find_peaks
import numpy as np
from datetime import datetime, timedelta

#from sklearn import metrics
import json
from matplotlib.patches import FancyArrowPatch


def get_stressor_directories(magnitude=None, duration=None, time_in_scope=None, directory='../output/HS'): 
    if magnitude != None: 
        directory_name = f'{magnitude:.2f}_'
    else: 
        directory_name = f'*_'
    
    if duration != None: 
        directory_name = f'{directory_name}{duration:.2f}'
    else: 
        directory_name = f'{directory_name}*_'
    
    if time_in_scope != None: 
        directory_name = f'{directory_name}_{time_in_scope:.2f}'
    else: 
        directory_name = f'{directory_name}*'

    return glob.glob(f'{directory}/{directory_name}/')



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


def plot_time_series(normal_signal, perturbed_signals, time_points, crh_signals = [], stress_times=None, save_dir = None):
    fig, ax1 = plt.subplots(figsize=(8, 4), constrained_layout=True)
    
    normal_signal = normal_signal[:2300]
    perturbed_signals = perturbed_signals[:2300]
    time_points = time_points[:2300]
    
    max = len(time_points)
    

    ax1.plot(time_points[:max], normal_signal[:max], label="Normal signal", color="blue", linewidth=2, linestyle="-")

    for i, perturbed_signal in enumerate(perturbed_signals):
        ax1.plot(time_points[:max], perturbed_signal[:max], linestyle="-", label=f"Perturbed signal", color="red")


    ax2 = ax1.twinx()

    for i, crh_signal in enumerate(crh_signals):
        #print(crh_signal)
        ax2.plot(time_points[:max], crh_signal[:max], linestyle="-", label=f"CRH", color="grey")

    if stress_times is not None: 
        for stress_time in stress_times: 
            ax1.axvline(stress_time, color='red', label='Time of stressor', linestyle='--')
            start_time = datetime.strptime("09:00:00", "%H:%M:%S")  
            stress_time_label = (start_time + timedelta(minutes=int(stress_time))).strftime("%H:%M")
            ax1.text(stress_time, ax1.get_ylim()[1] * 0.95, stress_time_label, color="red", fontsize=10, ha="right")
            ax1.axvspan(stress_time, stress_time + duration, color='red', alpha=0.2, label='Surgery duration' if i == 0 else None)


    ax1.set_xlabel("Time (hours)", fontsize=12)
    ax1.set_ylabel("Cortisol (nmol/L)", fontsize=12)
    ax2.set_ylabel("Circadian drive (CRH)", fontsize=12)
    ax2.spines['right'].set_color('grey')    
    ax2.tick_params(axis='y', colors='grey')
    ax2.yaxis.label.set_color('grey')
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    ax1.legend(handles1 + handles2, labels1 + labels2, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=5)

    #ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4)
    plt.grid()


    start_time='09:00:00'
    start_time = datetime.strptime(start_time, "%H:%M:%S")
    ax1.set_xticks([i for i in range(0, max + 1, 360)])
    ax1.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(0, max+1, 360)])


    if save_dir: 
        plt.savefig(save_dir)
    else: 
        return fig


all_directories = get_stressor_directories()
times = np.linspace(0, 1440*3, 1440*3)

combined_metrics = []

for i in range(len(all_directories)): 
    print(f'Now doing {all_directories[i]}')
    magnitude, duration, timing = all_directories[i].split('\\')[-2].split('_')
    magnitude, duration, timing = map(int, (float(magnitude), float(duration), float(timing)))

    signal, pert = get_signals_from_directory(all_directories[i])
    signal = signal[:1440*3]
    pert = pert[:1440*3]

    crh_signal = np.load(f'{all_directories[i]}/crh_stressor.npy')
    crh_signal = crh_signal[:1440*3]

    fig = plot_time_series(signal, [pert], times, stress_times = [timing], crh_signals=[crh_signal])

    cut_pert = pert[timing:timing+720]


    peaks, _ = find_peaks(cut_pert)
    print(peaks)
    number_of_peaks = len(peaks)


    plt.savefig(f'{all_directories[i]}/phase_response.png')
    plt.tight_layout()
    plt.close(fig)


    append_data = {
            "magnitude": magnitude,
            "duration": duration, 
            "time_in_scope": timing,
            "num_peaks": number_of_peaks
    }
    combined_metrics.append(append_data)

metrics_df = pd.DataFrame(combined_metrics)

print(metrics_df)
metrics_df.to_csv('ddemodel/stressor_output/results_heart_surgery.csv', index=False)
