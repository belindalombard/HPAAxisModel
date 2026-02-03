import imageio, glob, os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from io import BytesIO
from PIL import Image
import argparse

def plot_time_series(normal_signal, perturbed_signals, time_points, crh_signals = [], stress_times=None, save_dir = None, ymax = None, ymax_crh=None):
    fig, ax1 = plt.subplots(figsize=(8, 4), constrained_layout=True)
    
    normal_signal = normal_signal[:2300]
    perturbed_signals = perturbed_signals[:2300]
    time_points = time_points[:2300]
    
    max = len(time_points)
    

    ax1.plot(time_points[:max], normal_signal[:max], label="Normal signal", color="orange", linewidth=2, linestyle="-")

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
            if ymax: 
                position = ymax/2
            else: 
                position = ax1.get_ylim()[1] * 0.95
            ax1.text(stress_time, position , stress_time_label, color="red", fontsize=10, ha="right")
            ax1.axvspan(stress_time, stress_time + duration, color='red', alpha=0.2, label='Surgery duration' if i == 0 else None)
    if ymax: 
        ax1.set_ybound([0, ymax])
    if ymax_crh: 
        ax2.set_ybound([0, ymax_crh])
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


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate animation of stressor simulations')
parser.add_argument('--scenario', type=str, default='ACUTE',
                    help='Scenario subdirectory (e.g., ACUTE, HS1, HS2, HS3, HS4). Default: ACUTE')
parser.add_argument('--base-dir', type=str, default=None,
                    help='Base directory for simulations. Default: ../output/{scenario}')
parser.add_argument('--output-dir', type=str, default=None,
                    help='Output directory for animations. Default: ../output/animations')
args = parser.parse_args()

# Set paths based on scenario
base_directory = args.base_dir if args.base_dir else f'../output/{args.scenario}'
output_directory = args.output_dir if args.output_dir else '../output/animations'

def get_stressor_directories(magnitude=None, duration=None, time_in_scope=None, directory=None): 
    if directory is None:
        directory = base_directory
    
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

directories = get_stressor_directories(directory='../output/ACUTE', duration=240)
directories.sort(key=lambda x: (float(x.split('\\')[-2].split('_')[0]), float(x.split('\\')[-2].split('_')[2])))

images = []
times = np.linspace(0, 1440*3, 1440*3)

max_yaxis = 0
max_yaxis_crh = 0
for i, dir in enumerate(directories): 
    pert = np.load(f'{dir}/simulated_values_stressor.npy')[:, 1]
    crh_signal = np.load(f'{dir}/crh_stressor.npy')
    if np.max(pert) > max_yaxis: 
        max_yaxis = np.max(pert)
    if np.max(crh_signal) > max_yaxis_crh: 
        max_yaxis_crh = np.max(crh_signal)

max_yaxis += 20
max_yaxis_crh += 20


for i, dir in enumerate(directories): 
    print(f'Doing {i+1}/{len(directories)}.')
    magnitude, duration, timing = dir.split('\\')[-2].split('_')
    magnitude, duration, timing = map(int, (float(magnitude), float(duration), float(timing)))

    signal = np.load(f'{dir}/simulated_values_original.npy')[:, 1]
    pert = np.load(f'{dir}/simulated_values_stressor.npy')[:, 1]

    signal = signal[:1440*3]
    pert = pert[:1440*3]

    crh_signal = np.load(f'{dir}/crh_stressor.npy')
    crh_signal = crh_signal[:1440*3]

    fig = plot_time_series(signal, [pert], times, stress_times = [timing], crh_signals=[crh_signal], ymax=max_yaxis, ymax_crh=max_yaxis_crh)
    
    buf = BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    img = Image.open(buf)
    
    images.append(img)
    plt.close(fig)

# Ensure output directory exists
os.makedirs(output_directory, exist_ok=True)

output_path_mp4 = f'{output_directory}/output_animation_{args.scenario}.mp4'
writer = imageio.get_writer(output_path_mp4, fps=5)
for im in images:
    writer.append_data(np.array(im))
writer.close()

output_path_gif = f'{output_directory}/output_animation_{args.scenario}.gif'
imageio.mimsave(output_path_gif, images, fps=7)