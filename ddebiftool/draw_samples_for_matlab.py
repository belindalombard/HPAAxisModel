"""
draw_samples_for_matlab.py

Important to run get_csv_file_from_configs first, to ensure the CSV file exists with parameter values in ddebiftool/helpers.

Purpose:
    Computes statistics (median, mean, percentiles) for each parameter and generates JSON config files for MATLAB runs.
    Configs are saved in subdirectories of ddebiftool/run_pipeline/config/ where MATLAB expects them:
      - matlab/mean_median/ : configs with all median or mean parameters (mean_parameters.json, median_parameters.json)
      - matlab/percentiles/ : configs with one parameter at a percentile, rest at median
      - participants/ : configs for each participant (for simulation)
      - participants_bif/ : bifurcation configs for each participant
      - matlab/config.json and matlab/config_simulation.json are also saved for reference

Input:
    - ddebiftool/results/parameters.xlsx
    - ddebiftool/run_pipeline/config/config.json
    - ddebiftool/run_pipeline/config/config_simulation.json (optional)

Output:
    - ddebiftool/results/statistics_parameters.csv
    - ddebiftool/run_pipeline/config/matlab/mean_median/mean_parameters.json
    - ddebiftool/run_pipeline/config/matlab/mean_median/median_parameters.json
    - ddebiftool/run_pipeline/config/matlab/percentiles/{param}_25th_rest_median.json
    - ddebiftool/run_pipeline/config/matlab/percentiles/{param}_75th_rest_median.json
    - ddebiftool/run_pipeline/config/participants/Participant_{participant}.json
    - ddebiftool/run_pipeline/config/participants_bif/Participant_bif_{participant}.json
    - ddebiftool/run_pipeline/config/matlab/config.json
    - ddebiftool/run_pipeline/config/matlab/config_simulation.json

Usage:
    python draw_samples_for_matlab.py
"""

import pandas as pd
import numpy as np
import json
import os
import copy

base_dir = 'ddebiftool'
result_dir = f'{base_dir}/results'
config_base_dir = f'{base_dir}/run_pipeline/config'

# --- Create all output directories (matching MATLAB expectations) ---
matlab_dir = os.path.join(config_base_dir, 'matlab')
mean_median_dir = os.path.join(matlab_dir, 'mean_median')
percentile_dir = os.path.join(matlab_dir, 'percentiles')
participants_config_dir = os.path.join(config_base_dir, 'participants')
participants_bif_dir = os.path.join(config_base_dir, 'participants_bif')
os.makedirs(mean_median_dir, exist_ok=True)
os.makedirs(percentile_dir, exist_ok=True)
os.makedirs(participants_config_dir, exist_ok=True)
os.makedirs(participants_bif_dir, exist_ok=True)

# --- Load parameters ---
parameter_file = f'{result_dir}/parameters.xlsx'
parameters = pd.read_excel(parameter_file)

# --- Compute statistics ---
median_results = {}
mean_results = {}
statistics_data = []

for column in parameters.columns:
    values = parameters[column].dropna()
    median = np.median(values)
    mean = np.mean(values)
    percentile_25 = np.percentile(values, 25)
    percentile_75 = np.percentile(values, 75)
    percentile_5 = np.percentile(values, 5)
    percentile_95 = np.percentile(values, 95)
    median_results[column] = median
    mean_results[column] = mean
    statistics_data.append({
        'Parameter': column,
        'Median': median,
        'Mean': mean,
        '25th Percentile': percentile_25,
        '75th Percentile': percentile_75,
        '5th Percentile': percentile_5,
        '95th Percentile': percentile_95
    })

statistics_df = pd.DataFrame(statistics_data)
statistics_df.to_csv(f'{result_dir}/statistics_parameters.csv', index=False)

with open(os.path.join(mean_median_dir, 'median_parameters.json'), "w") as outfile:
    json.dump(median_results, outfile)
with open(os.path.join(mean_median_dir, 'mean_parameters.json'), "w") as outfile:
    json.dump(mean_results, outfile)

# --- Load base configs ---
with open(f'{config_base_dir}/config.json', 'r') as file:
    base_config = json.load(file)
# Check if config_simulation.json exists, otherwise use config.json as template
try:
    with open(f'{config_base_dir}/config_simulation.json', 'r') as file:
        base_config_sim = json.load(file)
except FileNotFoundError:
    print('config_simulation.json not found, using config.json as template')
    base_config_sim = copy.deepcopy(base_config)

# Save copies of the base configs for reference
with open(os.path.join(matlab_dir, 'config.json'), 'w') as f:
    json.dump(base_config, f, indent=4)
with open(os.path.join(matlab_dir, 'config_simulation.json'), 'w') as f:
    json.dump(base_config_sim, f, indent=4)

# --- Median/mean configs ---
required_parameters = ['k_a', 'k_c', 'alpha', 'm_a', 'm_c', 'gamma_a', 'gamma_c']

# Median config
median_config = copy.deepcopy(base_config)
for parameter in required_parameters:
    if parameter in statistics_df['Parameter'].values:
        row = statistics_df[statistics_df['Parameter'] == parameter]
        median_config['initial_parameters'][parameter] = float(row['Median'].iloc[0])
    else:
        match parameter:
            case 'm_a': median_config['initial_parameters'][parameter] = 4
            case 'm_c': median_config['initial_parameters'][parameter] = 4
            case 'gamma_a': median_config['initial_parameters'][parameter] = 0.0462
            case 'gamma_c': median_config['initial_parameters'][parameter] = 0.00906
            case _: print(f'Could not match {parameter}.')
with open(os.path.join(mean_median_dir, 'median_parameters.json'), "w") as outfile:
    json.dump(median_config, outfile, indent=4)

# Mean config
mean_config = copy.deepcopy(base_config)
for parameter in required_parameters:
    if parameter in statistics_df['Parameter'].values:
        row = statistics_df[statistics_df['Parameter'] == parameter]
        mean_config['initial_parameters'][parameter] = float(row['Mean'].iloc[0])
    else:
        match parameter:
            case 'm_a': mean_config['initial_parameters'][parameter] = 4
            case 'm_c': mean_config['initial_parameters'][parameter] = 4
            case 'gamma_a': mean_config['initial_parameters'][parameter] = 0.0462
            case 'gamma_c': mean_config['initial_parameters'][parameter] = 0.00906
            case _: print(f'Could not match {parameter}.')
with open(os.path.join(mean_median_dir, 'mean_parameters.json'), "w") as outfile:
    json.dump(mean_config, outfile, indent=4)

# --- Percentile configs (one parameter at 25th/75th, rest at median) ---
for index, row in statistics_df.iterrows():
    param_name = row['Parameter']
    if param_name in ('lambda_s', 'delay', 'lambda_a', 'Participant', 't_s', 'sigma'):
        continue
    # 25th
    config_25th = copy.deepcopy(median_config)
    config_25th['initial_parameters'][param_name] = row['25th Percentile']
    config_25th['filename'] = f"{param_name}_25th_rest_median"
    with open(os.path.join(percentile_dir, f"{param_name}_25th_rest_median.json"), "w") as outfile:
        json.dump(config_25th, outfile, indent=4)
    # 75th
    config_75th = copy.deepcopy(median_config)
    config_75th['initial_parameters'][param_name] = row['75th Percentile']
    config_75th['filename'] = f"{param_name}_75th_rest_median"
    with open(os.path.join(percentile_dir, f"{param_name}_75th_rest_median.json"), "w") as outfile:
        json.dump(config_75th, outfile, indent=4)

# --- Participant configs for simulation ---
parameters_to_include = ['k_c', 'm_a', 'k_a', 'm_c', 'alpha', 'delay', 'gamma_a', 'gamma_c']
crh_parameters = ['lambda_s', 'sigma', 't_s', 'lambda_a']

for index, participant_row in parameters.iterrows():
    participant = participant_row['Participant']
    config_file = copy.deepcopy(base_config_sim)
    for param, value in participant_row.items():
        if param in crh_parameters:
            config_file['crh_parameters'][param] = participant_row[param]
        if param in parameters_to_include:
            config_file['simulation_parameters'][param] = participant_row[param]
    config_file['filename'] = f'participant_{participant}_bif'
    config_file['participant_number'] = participant
    with open(os.path.join(participants_config_dir, f"Participant_{participant}.json"), "w") as outfile:
        json.dump(config_file, outfile, indent=4)

# --- Bifurcation configs for each participant ---
bif_parameters = ['k_c', 'm_a', 'k_a', 'm_c', 'alpha']
for index, participant_row in parameters.iterrows():
    participant = participant_row['Participant']
    config_file = copy.deepcopy(base_config)
    for param, value in participant_row.items():
        if param in bif_parameters:
            config_file['initial_parameters'][param] = participant_row[param]
    config_file['filename'] = f'participant_{participant}_bif'
    config_file['participant_number'] = participant
    with open(os.path.join(participants_bif_dir, f"Participant_bif_{participant}.json"), "w") as outfile:
        json.dump(config_file, outfile, indent=4)

# --- END OF PIPELINE ---

