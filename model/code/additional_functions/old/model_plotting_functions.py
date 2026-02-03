import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os 

def seperate_axes_plot(times
                       , simulated_values
                       , num_days=6
                       , participant_number = 0
                       , show_observed_values=False
                       , acth_data_file=''
                       , cort_data_file=''
                       , config = {}
                       , sleep_start = None
                       , sleep_end = None
                       ): 
    x_axis = [x - 1440*(num_days-1) for x in times[((num_days-1)*1440):]]
    print(x_axis)

    fig, ax = plt.subplots(2, 1, figsize=(5, 8))

    if show_observed_values is True:
        if acth_data_file == '' or cort_data_file == '': 
            print(f'Your datafiles are not defined. Can not plot observed values.')
            return

        acth_data = pd.read_csv(acth_data_file)
        cort_data = pd.read_csv(cort_data_file)

        observed_values = [] 
        observed_values.append(acth_data.loc[acth_data['test'] == participant_number, 'y'].tolist())
        observed_values.append(cort_data.loc[cort_data['test'] == participant_number, 'y'].tolist())
        observed_values = np.array(observed_values)
        observed_values = observed_values.T
    plt.rcParams.update({'font.size': 10})


    #fig, ax = plt.subplots(2, figsize=(5, 4))
    if show_observed_values is True: 
        ax[0].plot(x_axis, observed_values[:, 0], 'b-', label='Actual values ACTH')

    ax[0].plot(x_axis, simulated_values[:, 0], 'b--', label='Simulated values ACTH')
    ax[0].set_ylabel('ACTH (nmol/L)')
    #ax[0].legend(prop={'size': 8})
    ax[0].yaxis.set_label_coords(-0.12, 0.5)

    if show_observed_values is True: 
        ax[1].plot(x_axis, observed_values[:, 1], 'g-', label='Actual values CORT')

    ax[1].plot(x_axis, simulated_values[:, 1], 'g--', label='Simulated values CORT')
    ax[1].set_ylabel('Cortisol (nmol/L)')
    ax[1].set_xlabel('Time')
    ax[1].yaxis.set_label_coords(-0.12, 0.5)
    

    x_axis_time_of_day = [(x + 540) % 1440 for x in x_axis] 
    x_axis_time_of_day = [f"{int(time // 60):02d}:{int(time % 60):02d}" for time in x_axis_time_of_day]  
    x_axis_time_of_day.append("09:00")
    ax[0].set_xticks(np.arange(0, len(x_axis_time_of_day), step=1440 // 6))  
    ax[0].tick_params(axis='x', labelbottom=False)


    ax[1].set_xticks(np.arange(0, len(x_axis_time_of_day), step=1440 // 6))  
    ax[1].set_xticklabels(x_axis_time_of_day[::1440 // 6])  
    
    ax[0].grid(True)
    ax[1].grid(True)
    
    if sleep_start is not None and sleep_end is not None:
        sleep_start_min = sleep_start % 1440
        sleep_end_min = sleep_end % 1440
        estimated_days = int(np.ceil(len(x_axis) / 1440))
        for day in range(estimated_days):
            base = day * 1440
            if sleep_start_min < sleep_end_min:
                start = base + sleep_start_min
                end = base + sleep_end_min
                for axis in ax:
                    axis.axvspan(start, end, color='grey', alpha=0.3)
            else:
                for axis in ax:
                    axis.axvspan(base + sleep_start_min, base + 1440, color='grey', alpha=0.3)
                    axis.axvspan(base, base + sleep_end_min, color='grey', alpha=0.3)

    plt.title(config.get('plot_title', f'Participant {participant_number}'))    
    return fig 

def pertubation_plot(times
                       , simulated_values
                       , num_days=6
                       , config = {}
                       , crh_values = []
                       , start_time = '09:00:00'
                       , lines = '-'
                       , simulated_values_original = []
                       , crh_values_original = []
                       , length_model =1440
                       , yaxis = None
                       , time_of_stressor = None
                       ): 

    x_axis = [x - length_model*(num_days-1) for x in times[((num_days-1)*length_model):]]  

    fig, ax1 = plt.subplots(figsize=(9, 4), constrained_layout=True)
    
    ax1.plot(x_axis, simulated_values[:, 0], f'b{lines}', label='Stressor ACTH')
    ax1.set_ylabel('ACTH values', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.plot(x_axis, simulated_values[:, 1], f'g{lines}', label='Stressor CORT')
    ax2.set_ylabel('Cortisol values', color='green')
    ax2.tick_params(axis='y', labelcolor='green')

    ax3 = ax1.twinx()
    ax3.spines['right'].set_position(('outward', 60))  
    ax3.plot(x_axis, crh_values, f'k{lines}', label='Stressor CRH', color='black')

    if time_of_stressor: 
        ax3.axvline(time_of_stressor, color='red', label='Time of stressor', linestyle='--')

    ax3.set_ylabel('CRH values', color='black')
    ax3.tick_params(axis='y', labelcolor='black')

    if len(simulated_values_original) > 0:
        lines = '--'
        ax1.plot(x_axis, simulated_values_original[:, 0], f'b{lines}', label='Original ACTH')
        ax2.plot(x_axis, simulated_values_original[:, 1], f'g{lines}', label='Original CORT')
        if len(crh_values_original) > 0: 
            ax3.plot(x_axis, crh_values_original, f'k{lines}', label='Original CRH', color='gray')

    ax1.set_xlabel('Time')

    start_time = datetime.strptime(start_time, "%H:%M:%S")
    ax1.set_xticks([i for i in range(0, length_model+1, 360)])
    ax1.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(0, length_model+1, 360)])
    
    # Grouped legends
    signal_labels = ['Stressor ACTH', 'Original ACTH', 'Stressor CORT', 'Original CORT']
    crh_labels = ['Stressor CRH', 'Original CRH']
    time_of_stressor_labels = ['Time of stressor']
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines3, labels3 = ax3.get_legend_handles_labels()

    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center',
               bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize='small', frameon=False)
    ax3.legend(lines3, labels3, loc='upper center',
               bbox_to_anchor=(0.5, -0.2), ncol=1, fontsize='small', frameon=False)
    
    if time_of_stressor:
        ax3.legend(time_of_stressor_labels, loc='upper center', 
                   bbox_to_anchor=(0.5, -0.25), ncol=1, fontsize='small', frameon=False)

    if yaxis: 
        ax1.set_ybound(0, yaxis[0])
        ax2.set_ybound(0, yaxis[1])
        ax3.set_ybound(0, yaxis[2])

    plt.title(config.get('plot_title', 'Pertubation plot'))
    # plt.savefig('output/high_res_model.png', dpi=600)  # Commented out - caller saves to appropriate location
    #plt.show()
    return fig, (ax1, ax2, ax3)


def plot_simulation_with_crh(times
                       , simulated_values
                       , num_days=6
                       , participant_number=0
                       , show_observed_values=False
                       , acth_data_file=''
                       , cort_data_file=''
                       , config={}
                       , length_model=1440
                       , start_time='09:00:00'
                       , crh=[]
                       , lines_1='-'
                       , lines_2='--'
                       , outside_bounds=None
                       , sleep_start = None
                       , sleep_end = None
                       ): 

    x_axis = [x - length_model * (num_days - 1) for x in times[((num_days - 1) * length_model):]]

    fig, ax1 = plt.subplots(figsize=(8, 4), constrained_layout=True)

    ax1.plot(x_axis, simulated_values[:, 0], f'b{lines_2}', label='Simulated ACTH', color='blue')
    ax1.set_ylabel('ACTH (pmol/L)', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.plot(x_axis, simulated_values[:, 1], f'g{lines_2}', label='Simulated CORT', color='red')
    ax2.set_ylabel('CORT (nmol/L)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    if outside_bounds:
        ax1_min = max(min(simulated_values[:, 0]) - outside_bounds[0], 0)
        ax1_max = max(simulated_values[:, 0]) + outside_bounds[0]
        ax1.set_ylim(ax1_min, ax1_max)

        ax2_min = max(min(simulated_values[:, 1]) - outside_bounds[1], 0)
        ax2_max = max(simulated_values[:, 1]) + outside_bounds[1]
        ax2.set_ylim(ax2_min, ax2_max)

    if len(crh) > 0: 
        ax3 = ax1.twinx()
        ax3.spines['right'].set_position(('outward', 70))  
        ax3.plot(x_axis, crh, f'k{lines_2}', label='CRH drive', color='black')
        ax3.set_ylabel('CRH (arbritary)', color='black')
        ax3.tick_params(axis='y', labelcolor='black')

    if show_observed_values:
        if acth_data_file == '' or cort_data_file == '':
            print(f"Your data files are not defined. Cannot plot observed values.")
            return

        acth_data = pd.read_csv(acth_data_file)
        cort_data = pd.read_csv(cort_data_file)

        observed_values = []
        observed_values.append(acth_data.loc[acth_data['test'] == participant_number, 'y'].tolist())
        observed_values.append(cort_data.loc[cort_data['test'] == participant_number, 'y'].tolist())
        observed_values = np.array(observed_values).T

        ax1.plot(x_axis, observed_values[:, 0], f'b{lines_1}', label='Observed ACTH', color='blue')
        ax2.plot(x_axis, observed_values[:, 1], f'g{lines_1}', label='Observed CORT', color='red')

    start_time = datetime.strptime(start_time, "%H:%M:%S")
    ax1.set_xticks([i for i in range(0, length_model + 1, 360)])
    ax1.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(0, length_model + 1, 360)])
    ax1.set_xlabel('Time')

    if sleep_start and sleep_end:
        num_minutes = len(x_axis)
        estimated_days = int(np.ceil(num_minutes / 1440))

        for day in range(estimated_days):
            base = x_axis[0] + day * 1440
            if sleep_start < sleep_end:
                ax1.axvspan(base + sleep_start, base + sleep_end, color='grey', alpha=0.3)
            else:
                ax1.axvspan(base + sleep_start, base + 1440, color='grey', alpha=0.3)
                ax1.axvspan(base, base + sleep_end, color='grey', alpha=0.3)


    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    if len(crh)>0:
        lines3, labels3 = ax3.get_legend_handles_labels()
        ax1.legend(lines1 + lines2 + lines3, labels1 + labels2 + labels3, loc='upper center', 
            bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize='small', frameon=False)
    else: 
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center', 
            bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize='small', frameon=False)
        
    plt.title(config.get('plot_title', f'Participant {participant_number}'))
    #plt.savefig('output/high_res.pdf')
    #plt.savefig('output/high_res_combined_plot.png', dpi=600)
    #plt.show()

def plot_crh_only(times
                       , length_model=1440
                       , start_time='09:00:00'
                       , crh=[]
                        , num_days = 6
                       ): 

    x_axis = [x - length_model * (num_days - 1) for x in times[((num_days - 1) * length_model):]]

    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)


    ax.plot(x_axis, crh, f'k-', label='CRH drive', color='black')
    ax.set_ylabel('CRH Drive', color='black')
    ax.tick_params(axis='y', labelcolor='black')


    start_time = datetime.strptime(start_time, "%H:%M:%S")
    ax.set_xticks([i for i in range(0, length_model + 1, 360)])
    ax.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(0, length_model + 1, 360)])
    ax.set_xlabel('Time')


    #plt.title(config.get('plot_title', f'Participant {participant_number}'))
    plt.savefig('output/crh_plot.png', dpi=600)
    #plt.show()
