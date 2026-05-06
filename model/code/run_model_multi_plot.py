"""
Multi-Plot Runner for DDE Model
================================
Runs simulation once, generates multiple plot types.
More efficient than running run_model.py multiple times.

"""

import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add the repo root (HPAAxisModel/) to path so that absolute imports like
# `from model.code.classes.X import ...` resolve correctly regardless of cwd.
current_dir = os.path.dirname(os.path.abspath(__file__))   # model/code/
repo_root = os.path.dirname(os.path.dirname(current_dir))  # HPAAxisModel/
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from model.code.classes.model_class import Model
import model.code.additional_functions.model_plotting_functions_modular as mpf

# PROJECT_ROOT is model/code/ (used for relative data paths inside configs)
PROJECT_ROOT = current_dir

def run_simulation_and_plot(
    config_file,
    plot_options=[1],
    output_dir=None
):
    """
    Run simulation once and generate multiple plot types.
    
    Args:
        config_file: Path to config JSON file
        plot_options: List of plot types to generate  e.g. [1, 3, 4, 5]
        output_dir: Directory to save plots    
    Returns:
        dict with 'success', 'plot_files', 'simulated_values'
    """
    
    # Load config
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    # Get parameters from config
    parameters = config.get('parameters', {})
    fixed_params = config.get('fixed_params', {})
    signal = config.get('signal', 'both')
    participant_number = config.get('participant', 2)
    num_days = config.get('num_days', 4)
    reject = config.get('reject', True)
    show_observed_values = config.get('show_observed_values', True)
    plot_crh = config.get('plot_crh', False)
    
    # Data files (relative to repo root)
    acth_data_file = os.path.join(repo_root, 'data_analysis', 'data', 'data_shifted_ACTH_1min_09_00.csv')
    cort_data_file = os.path.join(repo_root, 'data_analysis', 'data', 'data_shifted_Cortisol_1min_09_00.csv')
    
    try:
        participant_number_int = int(participant_number)
    except (ValueError, TypeError):
        participant_number_int = participant_number
    
    print(f'Config: {config_file}')
    print(f'Participant: {participant_number_int}, Days: {num_days}')

    dde_model = Model(
            fixed_params=fixed_params,
            suggested_params=parameters,
            signal=signal,
            num_days=num_days,
            reject=reject
        )
 
    # Create time array
    timesteps = 1440 * num_days
    times = np.arange(0, timesteps, 1)

    # RUN SIMULATION.
    print(f"Running simulation...")
    parameters = dde_model.suggested_parameters()
    simulated_values = dde_model.simulate(parameters, times)
    print(f"Simulation complete! Signal (ACTH) length: {len(simulated_values[:, 0])}")
    simulated_values_full = simulated_values
    crh_values = None
    needs_crh = plot_crh 
    if needs_crh:
        print(f'Calculating CRH values...')
        crh_values = [dde_model.crh(t) for t in times]

    plot_files = []
    plot_names = {
        1: "Separate Axes",
        3: "Combined Plot", 
        4: "Combined with Bounds",
        5: "CRH Plot"
    }
    
    for plot_opt in plot_options:
        print(f"Generating plot {plot_opt}: {plot_names.get(plot_opt, 'Unknown')}...")
        
        # Determine if showing all days
        show_all_days = False
 
        # Generate the plot
        if plot_opt == 1:
            plot_kwargs = {
                'times': times, 
                'simulated_values': simulated_values_full if show_all_days else simulated_values,
                'num_days': num_days,
                'show_all_days': show_all_days
            }

            plotter = mpf.SeparateAxesPlot(**plot_kwargs)
            plotter.plot(
                participant_number=participant_number_int,
                show_observed_values=show_observed_values,
                acth_data_file=acth_data_file,
                cort_data_file=cort_data_file,
                config=config,
                start_time='09:00'
            )
            print(f'Showing observed values {show_observed_values}')
        
        elif plot_opt == 2:
            show_observed_values = False
            lines_2 = '-' if not show_observed_values else '--'
            outside_bounds = None if plot_opt == 3 else (10, 40)
            
            plot_kwargs = {
                'times': times,
                'simulated_values': simulated_values_full if show_all_days else simulated_values,
                'num_days': num_days,
                'show_all_days': show_all_days
            }

            plotter = mpf.CombinedSimulationPlot(**plot_kwargs)
            plotter.plot(
                participant_number=participant_number_int,
                show_observed_values=show_observed_values,
                acth_data_file=acth_data_file,
                cort_data_file=cort_data_file,
                config=config,
                crh=crh_values,
                lines_2=lines_2,
                outside_bounds=outside_bounds
            )
        
        # Save plot
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            plot_file = os.path.join(output_dir, f'plot_{plot_opt}.png')
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plot_files.append(plot_file)
            print(f"   Saved: {plot_file}")
        
        plt.close('all')  # Close to free memory
    
    return {
        'success': True,
        'plot_files': plot_files,
        'simulated_values': simulated_values_full if show_all_days else simulated_values,
        'num_plots': len(plot_files)
    }

