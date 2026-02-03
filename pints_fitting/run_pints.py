#!/usr/bin/env python
# coding: utf-8
"""
HPA Axis Model - Parameter Fitting with PINTS
=============================================
This script fits the HPA axis base model parameters to observed ACTH and Cortisol data
using the PINTS optimization framework. Participant to fit is obtained from the config 
file. 

Usage:
    python pints_fitting/run_pints.py --config parameters.json
    python pints_fitting/run_pints.py --config parameters_9.json --phase_align 1
"""


import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
import pints
from ddeint import ddeint
import math
import pandas as pd
import imageio
import json
import argparse
import time
from model.code.classes.model_class import Model
import model.code.classes.error_function as ef
from model.code.classes.custom_logger import CustomLogger

# === Important configurations for optimizer === #
MAX_ITERATIONS = 1500  

def _paths():
    """Compute and return frequently used project paths.

    Returns a dict with: project_root, ACTH_DATA_FILE, CORT_DATA_FILE,
    OUTPUT_DIR.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..'))
    acth = os.path.join(project_root, 'model', 'data', 'data_shifted_ACTH_1min_09_00.csv')
    cort = os.path.join(project_root, 'model', 'data', 'data_shifted_Cortisol_1min_09_00.csv')
    out_dir = os.path.join(project_root, 'pints_fitting', 'output')
    return {
        'project_root': project_root,
        'ACTH_DATA_FILE': acth,
        'CORT_DATA_FILE': cort,
        'OUTPUT_DIR': out_dir,
    }

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Fit HPA Axis model parameters using PINTS')
    parser.add_argument('--animation', type=int, default=0, help='Generate animation (1) or not (0)')
    parser.add_argument('--phase_align', type=int, default=0, help='Align model output to data by shifting before calculating error (1 or 0)')
    parser.add_argument('--config', type=str, default='parameters.json', help='Configuration file name')

    args = parser.parse_args()
    return {
        'phase_align': args.phase_align, 
        'config_file': f'model/config/base/{args.config}', 
        'animation': args.animation
    }

# === Helper functions ===
def create_next_output_path(base_dir, participant_number):
    """Create the next output path for a given participant number."""
    os.makedirs(base_dir, exist_ok=True)
    existing = (
        int(d.split('.')[-1])
        for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d))
        and d.startswith(f"{participant_number}.")
        and d.split('.')[-1].isdigit()
    )
    next_seq = max(existing, default=0) + 1
    dir_name = f"{participant_number}.{next_seq}"
    participant_dir = os.path.join(base_dir, dir_name)
    os.makedirs(participant_dir, exist_ok=True)
    return participant_dir

def save_parameters(output_path, fitted_params, fixed_params, others, signal, participant_number):
    """Save fitted and fixed parameters to a text file."""
    file_path = output_path + '.txt'
    with open(file_path, 'w') as file:
        file.write(f'{signal} fitted\n')
        file.write(f'Participant {participant_number}\n\n')
        file.write('Fitted Parameters:\n')
        for param, value in fitted_params.items():
            file.write(f'  {param}={value}\n')
        file.write('\nFixed Parameters:\n')
        for param, value in fixed_params.items():
            file.write(f'  {param}={value}\n')
        file.write('\nOther Settings:\n')
        for param, value in others.items():
            file.write(f'  {param}={value}\n')

def main():
    paths = _paths()
    project_root = paths['project_root']
    print(f"Project root: {project_root}")

    # === Get command line arguments and set up accordingly ===
    config_information = parse_arguments()
    phase_align_flag = int(config_information.get('phase_align', 0))
    
    print(f"{'='*60}")
    print(f"HPA Axis Model - Parameter Fitting")
    print(f"{'='*60}")
    print(f"Phase alignment: {'ON' if phase_align_flag == 1 else 'OFF'}")
    print(f"{'='*60}\n")

    # === DATA AND CONFIGURATION FILES ===
    ACTH_DATA_FILE = paths['ACTH_DATA_FILE']
    CORT_DATA_FILE = paths['CORT_DATA_FILE']
    OUTPUT_DIR = paths['OUTPUT_DIR']

    # Get config file path
    config_file = os.path.join(project_root, config_information['config_file'])
    print(f'Configuration file: {config_file}')

    # === Get JSON configuration information ===
    with open(config_file, 'r') as file:
        config = json.load(file)

    # Extract from config file 
    parameters = config.get('parameters', {})
    fixed_params = config.get('fixed_params', {})
    signal = config.get('signal', 'Both')
    participant_number = config.get('participant', 9)
    reject = config.get('reject', True)
    # Ensure participant_number is always an int for data selection
    participant_number_int = int(participant_number) if str(participant_number).isdigit() else participant_number
    num_days = config.get('num_days', 8)

    print(f"Participant: {participant_number}")
    print(f"Signal: {signal}")
    print(f"Simulation days: {num_days}")

    # === Set up output information ===
    output_path = create_next_output_path(OUTPUT_DIR, participant_number)
    print(f"Output directory: {os.path.abspath(OUTPUT_DIR)}")
    print(f"Output path: {os.path.abspath(output_path)}\n")
    
    # Ensure output_path directory exists
    output_path_dir = os.path.abspath(output_path)
    if not os.path.exists(output_path_dir):
        os.makedirs(output_path_dir, exist_ok=True)

    # Save parameters that are being fitted
    save_parameters(
        os.path.join(output_path, 'parameters'), 
        fitted_params=parameters, 
        fixed_params=fixed_params, 
        others={"num_days": num_days}, 
        signal=signal, 
        participant_number=participant_number
    )

    # === Set up model & pints configuration === 
    print("Initializing HPA Axis model...")
    dde_model = Model(
        fixed_params=fixed_params, 
        suggested_params=parameters, 
        signal=signal, 
        num_days=num_days, 
        reject=reject
    )

    # Set up model boundaries
    parameter_boundaries = dde_model.get_and_create_boundaries()

    # === Set up observed data structure ===
    print("Loading observed data...")
    acth_data = pd.read_csv(ACTH_DATA_FILE)
    cort_data = pd.read_csv(CORT_DATA_FILE)
    values = []
    values.append(acth_data.loc[acth_data['test'] == participant_number_int, 'y'].tolist())
    values.append(cort_data.loc[cort_data['test'] == participant_number_int, 'y'].tolist())
    values = np.array(values)
    values = values.T

    times = np.arange(0, 1440 * num_days, 1.0)

    # === Assign error function based on phase_align_flag ===
    if phase_align_flag == 1:
        print("Using PhaseAlignedMSE error function (phase alignment ON)")
        error_function = ef.PhaseAlignedMSE(dde_model, times, values)
    else:
        print("Using ErrorMeasureTwoSignal error function (phase alignment OFF)")
        error_function = ef.ErrorMeasureTwoSignal(dde_model, times, values)

    # === Set up and run optimizer ===
    print(f"\nStarting optimization (max iterations: {MAX_ITERATIONS})...\n")
    initial_parameters = dde_model.suggested_parameters()
    logger = CustomLogger(
        error_function, 
        num_days, 
        signal, 
        output_dir=output_path, 
        plot_on_improvement=True
    )
    
    optimiser = pints.OptimisationController(
        error_function,
        initial_parameters,
        method=pints.CMAES,
        boundaries=parameter_boundaries
    )
    optimiser.set_max_iterations(MAX_ITERATIONS)
    optimiser.set_log_to_screen(True)
    optimiser.set_callback(logger)
    
    estimated_parameters, found_value = optimiser.run()
    
    print("\n" + "="*60)
    print("Optimization complete!")
    print("="*60)
    print("Estimated parameters:", estimated_parameters)
    print(f"Final error value: {found_value}")
    print("="*60 + "\n")

    # === Animation generation (optional) ===
    if config.get('animation', 0) == 1:
        print("Generating animation...")
        image_directory = output_path
        images = []
        for file_name in sorted(os.listdir(image_directory)):
            if file_name.endswith('.png') or file_name.endswith('.jpg'):
                file_path = os.path.join(image_directory, file_name)
                images.append(imageio.imread(file_path))
        output_path_gif = os.path.join(image_directory, 'output_animation.gif')
        imageio.mimsave(output_path_gif, images, duration=0.05)
        print(f"✓ Animation saved to {output_path_gif}\n")

    print("All tasks completed successfully!")

if __name__ == '__main__':
    import multiprocessing as mp
    mp.freeze_support()
    main()
