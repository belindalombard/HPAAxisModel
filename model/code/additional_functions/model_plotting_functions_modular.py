import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from datetime import datetime, timedelta
from scipy.interpolate import interp1d

plt.rcParams.update({
    'font.size': 12,         
    'axes.titlesize': 12,        
    'axes.labelsize': 12,        
    'xtick.labelsize': 10,        
    'ytick.labelsize': 10,     
    'legend.fontsize': 10,       
})

class BasePlot:
    def __init__(self, times, simulated_values, num_days=6, length_model=1440, show_all_days=False):
        self.times = times
        self.simulated_values = simulated_values
        self.num_days = num_days
        self.length_model = length_model
        self.show_all_days = show_all_days  # Show all days instead of just last day
        self.x_axis = self._calculate_x_axis()

    def _calculate_x_axis(self):
        if self.show_all_days:
            # Show all days from start to end - use full time array
            return list(self.times)
        else:
            # Show only the last day (original behavior)
            return [x - self.length_model * (self.num_days - 1) for x in self.times[((self.num_days - 1) * self.length_model):]]

    def _set_xticks(self, ax, start_time):
        start_dt = datetime.strptime(start_time, "%H:%M:%S")
        tick_range = range(0, self.length_model + 1, 360)
        tick_labels = [(start_dt + timedelta(minutes=i)).strftime("%H:%M") for i in tick_range]
        ax.set_xticks(tick_range)
        ax.set_xticklabels(tick_labels)
        ax.set_xlabel('Time of day')

    def _load_observed_values(self, participant_number, acth_file, cort_file):
        acth = pd.read_csv(acth_file).loc[lambda df: df['test'] == participant_number, 'y'].tolist()
        cort = pd.read_csv(cort_file).loc[lambda df: df['test'] == participant_number, 'y'].tolist()
        return np.array([acth, cort]).T



class SeparateAxesPlot(BasePlot):
    def plot(self, participant_number=0, show_observed_values=False, acth_data_file='', show_title=True, cort_data_file='', config={}, start_time='09:00', crh=None):
        has_crh = crh is not None
        
        # Adjust figure size if showing all days
        figsize = (15, 11) 
        fig, ax = plt.subplots(2, 1, figsize=figsize)

        if self.show_all_days:
            plot_data = self.simulated_values
            plot_times = self.times
        else:
            plot_data = self.simulated_values[-self.length_model:]
            plot_times = self.times[-self.length_model:]

        if show_observed_values and acth_data_file and cort_data_file:
            observed = self._load_observed_values(participant_number, acth_data_file, cort_data_file)
            # Note: observed data is only for last day, so only plot if not showing all days
            # Also check that observed data is not empty (e.g., for mean/median participants)
            if not self.show_all_days and len(observed) > 0:
                ax[0].plot(self.x_axis, observed[:, 0], linestyle='--', color='blue', label='Observed ACTH')
                ax[1].plot(self.x_axis, observed[:, 1], linestyle='--', color='red', label='Observed CORT')

        ax[0].plot(self.x_axis, plot_data[:, 0], linestyle='-', color='blue', label='Simulated ACTH')
        ax[1].plot(self.x_axis, plot_data[:, 1], linestyle='-', color='red', label='Simulated CORT')

        ax[0].set_ylabel('ACTH (nmol/L)')
        ax[1].set_ylabel('CORT (nmol/L)')
        ax[-1].set_xlabel('Time (minutes)' if self.show_all_days else 'Time')

  
        # Handle time axis ticks differently for all days vs single day
        if self.show_all_days:
            # For multiple days, show time-of-day markers (00:00, 06:00, 12:00, 18:00)
            base_time = datetime.strptime(start_time, "%H:%M")
            # Create ticks every 6 hours (360 minutes)
            tick_interval = 360
            ticks = np.arange(0, len(self.x_axis), step=tick_interval)
            tick_labels = []
            for tick in ticks:
                minutes_from_start = self.x_axis[int(tick)] if tick < len(self.x_axis) else self.x_axis[-1]
                current_time = base_time + timedelta(minutes=int(minutes_from_start))
                # Format as HH:MM with day number
                day_num = int(minutes_from_start / 1440) + 1
                time_str = current_time.strftime("%H:%M")
                tick_labels.append(f'D{day_num}\n{time_str}')
        else:
            tick_interval = 240
            ticks = np.arange(0, len(self.x_axis), step=tick_interval)
            base_time = datetime.strptime(start_time, "%H:%M")
            tick_labels = [(base_time + timedelta(minutes=int(self.x_axis[int(tick)]))).strftime("%H:%M") for tick in ticks]

        for axis in ax:
            axis.set_xticks(ticks)
            axis.set_xticklabels(tick_labels, rotation=45, fontsize=8 if self.show_all_days else 10)
            axis.grid(True)

        for axis in ax[:-1]:
            axis.tick_params(axis='x', labelbottom=False)
             
        if show_title: 
            fig.suptitle(config.get('plot_title', f'Participant {participant_number}'))
        plt.tight_layout()
        return fig
        


class CRHPlot(BasePlot):
    def plot(self, crh, start_time='09:00:00', 
             crh_suppression=None, annotate=True, crh_unsuppressed=None):
        fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)


        if crh_unsuppressed is not None:
            ax.plot(self.x_axis, crh, '-', label='Hypothalamic drive (modulated)', color='black')
            ax.plot(self.x_axis, crh_unsuppressed, '--', label='Hypothalamic drive (baseline)', color='gray', alpha=0.8)
        else: 
            ax.plot(self.x_axis, crh, '-', label='CRH drive', color='black')


        ax.set_ylabel('Hypothalamic drive ($\Phi$)', color='black')
        ax.tick_params(axis='y', labelcolor='black')
        self._set_xticks(ax, start_time)

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
        #plt.savefig('output/crh_suppression_overlay.png', dpi=600)
        return fig



class PerturbationPlot(BasePlot):
    def plot(self, crh_values, start_time='09:00:00', lines='-', simulated_values_original=None, crh_values_original=None, yaxis=None, time_of_stressor=None):
        fig, ax1 = plt.subplots(figsize=(9, 4), constrained_layout=True)
        ax1.plot(self.x_axis, self.simulated_values[:, 0], linestyle='-', color='blue', label='Stressor ACTH')
        ax1.set_ylabel('ACTH values', color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')

        ax2 = ax1.twinx()
        ax2.plot(self.x_axis, self.simulated_values[:, 1], linestyle='-', color='red', label='Stressor CORT')
        ax2.set_ylabel('CORT values', color='red')
        ax2.tick_params(axis='y', labelcolor='red')

        ax3 = ax1.twinx()
        ax3.spines['right'].set_position(('outward', 60))
        ax3.plot(self.x_axis, crh_values, f'k{lines}', label='Stressor CRH', color='black')

        if time_of_stressor:
            ax3.axvline(time_of_stressor, color='red', label='Time of stressor', linestyle='--')

        if simulated_values_original is not None:
            ax1.plot(self.x_axis, simulated_values_original[:, 0], linestyle='--', color='blue', label='Original ACTH')
            ax2.plot(self.x_axis, simulated_values_original[:, 1], linestyle='--', color='red', label='Original CORT')
            if crh_values_original is not None:
                ax3.plot(self.x_axis, crh_values_original, 'k--', label='Original CRH', color='gray')

        self._set_xticks(ax1, start_time)

        if yaxis:
            ax1.set_ybound(0, yaxis[0])
            ax2.set_ybound(0, yaxis[1])
            ax3.set_ybound(0, yaxis[2])

        fig.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize='small', frameon=False)
        plt.title('Perturbation plot')
        plt.savefig('output/high_res_model.png', dpi=600)
        return fig

class CombinedSimulationPlot(BasePlot):
    def plot(self, crh=None, show_observed_values=False, acth_data_file='', cort_data_file='', config={}, start_time='09:00:00',
         lines_1='-', lines_2='--', outside_bounds=None, participant_number=0,
         simulated_values_baseline=None, crh_baseline=None):
        fig, ax1 = plt.subplots(figsize=(8, 4), constrained_layout=True)
        ax1.plot(self.x_axis, self.simulated_values[:, 0], linestyle='-', color='blue', label='Simulated ACTH')
        ax1.set_ylabel('ACTH (pmol/L)', color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')

        ax2 = ax1.twinx()
        ax2.plot(self.x_axis, self.simulated_values[:, 1], linestyle='-', color='red', label='Simulated CORT')
        ax2.set_ylabel('CORT (nmol/L)', color='red')
        ax2.tick_params(axis='y', labelcolor='red')

        if simulated_values_baseline is not None:
            ax1.plot(self.x_axis, simulated_values_baseline[:, 0], linestyle=':', color='blue', label='Baseline ACTH')
            ax2.plot(self.x_axis, simulated_values_baseline[:, 1], linestyle=':', color='red', label='Baseline CORT')


        if outside_bounds:
            ax1.set_ylim(max(min(self.simulated_values[:, 0]) - outside_bounds[0], 0), max(self.simulated_values[:, 0]) + outside_bounds[0])
            ax2.set_ylim(max(min(self.simulated_values[:, 1]) - outside_bounds[1], 0), max(self.simulated_values[:, 1]) + outside_bounds[1])

        if crh is not None:
            ax3 = ax1.twinx()
            ax3.spines['right'].set_position(('outward', 70))
            ax3.plot(self.x_axis, crh[-len(self.x_axis):], f'k{lines_2}', label='CRH drive', color='black')
            ax3.set_ylabel('CRH (arbritary)', color='black')
            ax3.tick_params(axis='y', labelcolor='black')
            if crh_baseline is not None:
                ax3.plot(self.x_axis, crh_baseline, 'k:', label='Baseline CRH', color='gray')

        if show_observed_values and acth_data_file and cort_data_file:
            observed = self._load_observed_values(participant_number, acth_data_file, cort_data_file)
            # Only plot if observed data exists (not empty for mean/median)
            if len(observed) > 0:
                ax1.plot(self.x_axis, observed[:, 0], f'b{lines_1}', label='Observed ACTH', color='blue')
                ax2.plot(self.x_axis, observed[:, 1], f'g{lines_1}', label='Observed CORT', color='red')

        self._set_xticks(ax1, start_time)
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        legend_handles = handles1 + handles2
        legend_labels = labels1 + labels2

        if crh is not None:
            handles3, labels3 = ax3.get_legend_handles_labels()
            legend_handles += handles3
            legend_labels += labels3

        ax1.legend(legend_handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=False)
        plt.title(config.get('plot_title', f'Participant {participant_number}'))
        plt.tight_layout()
        return fig

