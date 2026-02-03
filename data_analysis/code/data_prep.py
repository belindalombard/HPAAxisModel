"""
HPA Axis Data Preparation
=========================
Prepares raw participant data for analysis by:
1. Loading and cleaning raw CSV files
2. Resampling to consistent time intervals
3. Shifting data to a common reference time (09:00)
4. Saving processed data for model fitting and wavelet analysis

Usage:
    python data_analysis/code/data_prep.py

Input:
    - data_analysis/data/raw/participant_*.csv files

Output:
    - data/orignal_data_combined.csv
    - data/data_resampled.csv
    - data/minute_data_cort.csv
    - data/data_per_participant_pyboat.csv
    - data/data_shifted_Cortisol_1min_09_00.csv
    - data/data_shifted_Cortisol_10min_09_00.csv
    - data/data_shifted_Cortisol_1min_00_00.csv
    - data/data_shifted_Cortisol_10min_00_00.csv
    - data/data_shifted_ACTH_1min_09_00.csv
    - data/data_shifted_ACTH_10min_09_00.csv
    - data/data_shifted_ACTH_1min_00_00.csv
    - data/data_shifted_ACTH_10min_00_00.csv
"""

import os
import pandas as pd
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import seaborn as sns


def prepare_raw_data(data_directory='data_analysis/data', interpolation='1T'):
    """
    Load and prepare raw participant data.
    
    Args:
        data_directory: Directory containing raw data
        interpolation: Resampling frequency (default: '1T' for 1 minute)
    
    Returns:
        tuple: (combined_data, combined_data_resampled)
    """
    raw_data_dir = os.path.join(data_directory, 'raw')
    data_frames_resampled = []
    data_frames = []

    print("Processing raw participant data...")
    
    for filename in os.listdir(raw_data_dir):
        if filename.endswith(".csv") and filename.startswith("participant_"):
            participant_id = filename.split('participant_')[1].split('.csv')[0]
            
            print(f"  Processing participant {participant_id}...")
            file_path = os.path.join(raw_data_dir, filename)
            df = pd.read_csv(file_path)
            
            selected_columns = ['Date', 'Time', 'Cortisol', 'ACTH']
            df_selected = df[selected_columns].copy() 

            # Convert ACTH units for participants > 7
            if int(participant_id) > 7: 
                conversion = 4.54113
                df_selected['ACTH'] = df_selected['ACTH'] / conversion

            df_selected['Time'] = pd.to_datetime(df_selected['Time'], format='%H:%M:%S').dt.time
            df_selected['Date'] = pd.to_datetime(df_selected['Date'], format='%d/%m/%Y')
            df_selected['Time'] = pd.to_datetime(df_selected['Date'].astype(str) + ' ' + df_selected['Time'].astype(str))

            df_selected.loc[:, 'participant_id'] = participant_id
            
            # Drop missing values
            df_selected = df_selected.dropna(subset=['Time', 'Cortisol', 'ACTH'])

            # Ensure numeric types
            df_selected['Cortisol'] = pd.to_numeric(df_selected['Cortisol'], errors='coerce')
            df_selected['ACTH'] = pd.to_numeric(df_selected['ACTH'], errors='coerce')

            df_selected = df_selected.sort_values(by='Time')
            df_selected.set_index('Time', inplace=True)

            # Resample and interpolate
            df_selected_resampled = df_selected[['Cortisol', 'ACTH']].resample(interpolation).mean().interpolate(method='cubic')
            df_selected_resampled.reset_index(inplace=True)
            df_selected.reset_index(inplace=True)

            df_selected_resampled.loc[:, 'participant_id'] = participant_id

            data_frames_resampled.append(df_selected_resampled)
            data_frames.append(df_selected)

    # Combine all participant data
    combined_data_resampled = pd.concat(data_frames_resampled, ignore_index=True)
    combined_data_resampled['participant_id'] = combined_data_resampled['participant_id'].astype(int)
    combined_data_resampled = combined_data_resampled.sort_values(by=['participant_id','Time'])
    combined_data_resampled = combined_data_resampled.reset_index(drop=True)

    combined_data = pd.concat(data_frames, ignore_index=True)
    combined_data['participant_id'] = combined_data['participant_id'].astype(int)
    combined_data = combined_data.sort_values(by='participant_id')
    combined_data = combined_data.reset_index(drop=True)

    return combined_data, combined_data_resampled


def shift_test_data(ref_time_str, test_data, min_times, signal='Cortisol', interpolation='10min', data_directory='data'):
    """
    Shift participant data to a common reference time.
    
    Args:
        ref_time_str: Reference time as string (e.g., "09:00")
        test_data: DataFrame with columns ['x', 'y', 'test']
        min_times: DataFrame with participant minimum times
        signal: Name of the signal being processed
        interpolation: Resampling frequency
        data_directory: Output directory
    
    Returns:
        DataFrame: Shifted and aligned data
    """
    ref_time = datetime.strptime(ref_time_str, "%H:%M").time()
    shifted_data = []  
    
    print(f"Shifting {signal} data to reference time {ref_time_str}...")
    
    for _, row in min_times.iterrows():
        participant_id = row['participant']
        
        participant_data = test_data[test_data['test'] == participant_id].copy()
        participant_data.reset_index(inplace=True)

        participant_data['date_interpolate'] = pd.Timestamp('2021-01-01') + pd.to_timedelta(
            participant_data['x'].dt.time.apply(lambda t: f"{t.hour:02d}:{t.minute:02d}:{t.second:02d}")
        )
        
        participant_data = participant_data.sort_values(by='date_interpolate', ascending=True).reset_index(drop=True)
        participant_data.set_index('date_interpolate', inplace=True)

        if isinstance(participant_data.index, pd.MultiIndex):
            participant_data.reset_index(level=0, inplace=True)  
        
        # Resample and interpolate
        resampled = participant_data.resample(interpolation).mean()
        for col in resampled.select_dtypes(include=[np.number]).columns:
            resampled[col] = resampled[col].interpolate(method='cubic')
        participant_data = resampled

        participant_data.reset_index(inplace=True)
        participant_data.rename(columns={'date_interpolate': 'x'}, inplace=True)

        # Remove duplicate columns
        participant_data = participant_data.loc[:,~participant_data.columns.duplicated()]
        participant_data['x'] = pd.to_datetime(participant_data['x'])

        def calculate_timesteps(time):
            time_dt = datetime.combine(datetime.today(), time)
            ref_time_dt = datetime.combine(datetime.today(), ref_time)
            if time >= ref_time:
                delta = (time_dt - ref_time_dt).total_seconds()
            else:
                delta = ((time_dt + timedelta(days=1)) - ref_time_dt).total_seconds()
            return int(delta // 60)  

        participant_data['timesteps_after_reftime'] = participant_data['x'].dt.time.apply(calculate_timesteps)
        shifted_data.append(participant_data)

    return_data = pd.concat(shifted_data)

    # Keep relevant columns
    keep_cols = [col for col in ['x', 'y', 'test', 'timesteps_after_reftime'] if col in return_data.columns]
    return_data = return_data[keep_cols] 
    if 'test' in return_data.columns:
        return_data['test'] = return_data['test'].round(0)

    return_data = return_data.sort_values(by=['timesteps_after_reftime', 'test'])

    # Create visualization
    g = sns.FacetGrid(return_data, col="test", col_wrap=3, height=4, sharey=True)
    g.map(sns.lineplot, "timesteps_after_reftime", "y")
    plt.savefig(f"{data_directory}/shifted_{signal}_plot.png", dpi=150, bbox_inches='tight')
    plt.close()

    # Save to file
    output_file = f"{data_directory}/data_shifted_{signal}_{interpolation}_{ref_time_str.replace(':', '_')}.csv"
    return_data.to_csv(output_file, index=False)
    print(f"  Saved: {output_file}")

    return return_data


def main():
    """Main data preparation pipeline."""
    data_directory = 'data_analysis/data'
    
    print("="*60)
    print("HPA Axis Data Preparation")
    print("="*60)
    print()
    
    # Step 1: Load and resample raw data
    combined_data, combined_data_resampled = prepare_raw_data(data_directory, interpolation='10T')
    combined_data, combined_data_resampled = prepare_raw_data(data_directory)

    # Save combined data
    print("\nSaving combined data files...")
    combined_data.to_csv(f'{data_directory}/orignal_data_combined.csv', index=False)
    print(f"  Saved: {data_directory}/orignal_data_combined.csv")
    
    combined_data_resampled.to_csv(f'{data_directory}/data_resampled.csv', index=False)
    print(f"  Saved: {data_directory}/data_resampled.csv")
    
    # Create pivot table for wavelet analysis
    extracted_data = combined_data_resampled.copy()
    extracted_data['Time'] = pd.to_datetime(extracted_data['Time'])
    extracted_data['x'] = ((extracted_data['Time'] - extracted_data.groupby('participant_id')['Time'].transform('min')).dt.total_seconds() / 60).astype(int)
    extracted_data = extracted_data[['participant_id', 'x', 'Cortisol', 'ACTH']]
    pivot_table = extracted_data.pivot_table(index='x', columns='participant_id', values=['Cortisol', 'ACTH'])
    pivot_table.columns = ['_'.join(map(str, col)).strip() for col in pivot_table.columns.values]
    pivot_table.reset_index(inplace=True)
    pivot_table.to_csv(f'{data_directory}/minute_data_cort.csv', index=False)
    pivot_table.to_csv(f'{data_directory}/data_per_participant_pyboat.csv', index=False)
    print(f"  Saved: {data_directory}/minute_data_cort.csv")
    print(f"  Saved: {data_directory}/data_per_participant_pyboat.csv")
    
    # Step 2: Shift data for Cortisol and ACTH for both 09:00 and 00:00 reference times, 1min and 10min interpolations
    test_data_cort = pd.read_csv(f'{data_directory}/data_resampled.csv')
    test_data_cort['Time'] = pd.to_datetime(test_data_cort['Time'])
    test_data_cort = test_data_cort.rename(columns={'Cortisol': "y", "Time": "x", 'participant_id': 'test'})[['x', 'y', 'test']]
    test_data_cort = test_data_cort.sort_values(['test', 'x'], ascending=[True, True])
    min_dates_cort = test_data_cort.groupby('test')['x'].min().reset_index()
    min_dates_cort.columns = ['participant', 'min_value']

    test_data_acth = pd.read_csv(f'{data_directory}/data_resampled.csv')
    test_data_acth['Time'] = pd.to_datetime(test_data_acth['Time'])
    test_data_acth = test_data_acth.rename(columns={'ACTH': "y", "Time": "x", 'participant_id': 'test'})[['x', 'y', 'test']]
    test_data_acth = test_data_acth.sort_values(['test', 'x'], ascending=[True, True])
    min_dates_acth = test_data_acth.groupby('test')['x'].min().reset_index()
    min_dates_acth.columns = ['participant', 'min_value']

    for ref_time in ['09:00', '00:00']:
        print(f"\nProcessing Cortisol data (1min interpolation, ref {ref_time})...")
        shift_test_data(ref_time, test_data_cort, signal='Cortisol', min_times=min_dates_cort, interpolation='1min', data_directory=data_directory)
        print(f"\nProcessing Cortisol data (10min interpolation, ref {ref_time})...")
        shift_test_data(ref_time, test_data_cort, signal='Cortisol', min_times=min_dates_cort, interpolation='10min', data_directory=data_directory)
        print(f"\nProcessing ACTH data (1min interpolation, ref {ref_time})...")
        shift_test_data(ref_time, test_data_acth, signal='ACTH', min_times=min_dates_acth, interpolation='1min', data_directory=data_directory)
        print(f"\nProcessing ACTH data (10min interpolation, ref {ref_time})...")
        shift_test_data(ref_time, test_data_acth, signal='ACTH', min_times=min_dates_acth, interpolation='10min', data_directory=data_directory)

    print("\n" + "="*60)
    print("Data preparation complete!")
    print("="*60)
    print(f"\nProcessed {len(combined_data_resampled['participant_id'].unique())} participants")
    print(f"Total data points: {len(combined_data_resampled)}")


if __name__ == "__main__":
    main()
