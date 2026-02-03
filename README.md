# HPA Axis Model

A comprehensive mathematical modeling & analysis framework for the Hypothalamic-Pituitary-Adrenal (HPA) axis using delay differential equations (DDEs). This repository contains everything needed to simulate, fit, analyze, and perform bifurcation analysis on the cortisol-ACTH feedback system.

## Overview

This repository provides a complete pipeline for HPA axis modeling:

1. **Data Preparation** - Process and prepare raw hormone time series data
2. **Wavelet Analysis** - Identify ultradian rhythms and periodic patterns
3. **Model Simulation** - Run the base HPA axis DDE model
4. **Parameter Fitting** - Fit model parameters to experimental data using PINTS
5. **Bifurcation Analysis** - Analyze system dynamics using DDE-BIFTOOL (MATLAB)
6. **Stressor Perturbations** - Simulate responses to external stress (acute stress, surgery, etc.)

The model describes the feedback interactions between ACTH and Cortisol, driven by CRH (external function).

## Getting Started

### Prerequisites

**Python Requirements:**
```bash
python >= 3.8
numpy
pandas
matplotlib
scipy
seaborn
pints
ddeint
pyboat  # For wavelet analysis
imageio  # Optional, for animations
```

**MATLAB Requirements** (for bifurcation analysis):
- MATLAB R2019b or later
- DDE-BIFTOOL package

### Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd HPAAxisModel
```

2. Install Python dependencies:
```bash
pip install numpy pandas matplotlib scipy seaborn pints ddeint pyboat imageio
```

3. For MATLAB bifurcation analysis, download and install [DDE-BIFTOOL](https://sourceforge.net/projects/ddebiftool/)

## Usage Guide

### 1. Data preparation

Process raw participant data and prepare for analysis:

```bash
cd data_analysis/code
python data_prep.py
```

**Input:** Raw CSV files in `data_analysis/data/raw/participant_*.csv`

**Output:**
- `data/data_resampled.csv` - Resampled time series
- `data/data_shifted_Cortisol_1min_09_00.csv` - Time-aligned cortisol data
- `data/data_shifted_ACTH_1min_09_00.csv` - Time-aligned ACTH data

Note: The timestamp & interpolation method can be adjusted as necessary in the python sciript.

### 2. Wavelet analysis

Identify ultradian rhythms and periodic patterns using a package called PyBOAT:

```bash
cd data_analysis/code
python wavelet_analysis.py --hormone Cortisol --dt 10
python wavelet_analysis.py --hormone ACTH --dt 10
```

**Arguments:**
- `--hormone`: Hormone to analyze (Cortisol or ACTH)
- `--dt`: Sampling interval in minutes (default: 10)
- `--data_folder`: Data directory (default: data)
- `--results_folder`: Output directory (default: output)

**Output:**
- Wavelet spectra for each participant
- Power distribution plots
- Statistical comparison matrices
- `output/stacked_wavelets_{hormone}.npy` - Saved wavelet data

### 3. Model Simulation

Run the HPA axis model with specified parameters:

```bash
cd model/code
python model/run_model.py --config config/parameters.json --plots 1,2 --output model/output/sim1
```

**Arguments:**
- `--config`: Path to parameter configuration file
- `--plots`: Plot types to generate (1=Separate, 2=Combined)
- `--output`: Output directory for saving plots
- `--no_show`: Skip displaying plots (only save)

**Example config file** (`config/parameters.json`):
```json
{
  "parameters": {
    "p1": 2.5,
    "p2": 1.8,
    "p3": 3.2
  },
  "fixed_params": {
    "T": 1440
  },
  "num_days": 3, 
  "signal": "Both"
}
```

### 4. Parameter Fitting

Fit model parameters to experimental data using PINTS optimization:

```bash
python pints_fitting/run_pints.py --config parameters.json 
```

**Arguments:**
- `--config`: Configuration file with initial parameter guesses
- `--phase_align`: Align model to data by phase shifting (1 or 0). In this case, do not fit t_s as well. 
- `--animation`: Generate optimization animation (1 or 0)

**Output:**
- `output/{participant}.{run}/` - Fitted parameters and plots
- `output/{participant}.{run}/parameters.txt` - Fitted parameter values

**Fitting Process:**
1. Load observed ACTH and Cortisol data
2. Initialize model with parameter guesses from config
3. Run CMA-ES optimization to minimize error
4. Save best-fit parameters and generate comparison plots

Note: This was run on an HPC cluster

### 5. Bifurcation Analysis

Analyze system dynamics and stability using DDE-BIFTOOL in MATLAB:

#### Step 1: Generate Configuration Files

```bash
python ddebiftool/get_csv_file_from_configs.py
python ddebiftool/draw_samples_for_matlab.py
```

This creates configuration files for MATLAB bifurcation analysis.

#### Step 2: Run MATLAB Analysis

```matlab
cd ddebiftool/run_pipeline
run_all_configs  % Run bifurcation analysis for all configurations
```

**Main MATLAB Scripts:**
- `run_bif_participants.m` - Individual participant bifurcation analysis
- `run_all_configs.m` - Batch processing
- `run_circ.m` - Circadian rhythm analysis

**Output:**
- Bifurcation diagrams
- Stability boundaries
- Period doubling cascades
- Limit cycle analysis

### 6. Stressor Perturbations

Simulate how the HPA axis responds to external stressors (e.g., acute stress, heart surgery). This allows you to explore how stressor timing, magnitude, and duration affect cortisol and ACTH dynamics.

**How it works:** The model runs normally for 1 day to establish baseline dynamics, then applies a stressor (CRH pulse) at a specified time. You can vary the stressor's magnitude (intensity), duration (how long it lasts), circadian phase (what time of day), and ultradian phase (where in the natural cortisol rhythm it occurs).

#### Run Stressor Simulations

```bash
cd stressors
python run_model_stressor.py --scenario ACUTE --magnitude 50 --duration 60 --start_time 14580
```

**Arguments:**
- `--scenario`: Type of stressor (ACUTE, HS1, HS2, HS3, HS4)
- `--magnitude`: Stressor intensity (CRH units)
- `--duration`: How long stressor lasts (minutes)
- `--start_time`: When to apply stressor (minutes after 09:00)
- `--config`: Model parameters file (optional)

**What each scenario means:**
- **ACUTE**: Single acute stressor (e.g., public speaking, exam)
- **HS1-HS4**: Heart surgery scenarios with different stressor profiles

**Output:** Each simulation saves to `stressors/output/{scenario}/mag{X}_dur{Y}_start{Z}/`:
- `simulated_values_original.npy` - Normal HPA axis dynamics (no stressor)
- `simulated_values_stressor.npy` - HPA dynamics with stressor
- `crh_stressor.npy` - CRH stressor profile over time
- `metrics.csv` - Response metrics (peak increase, time to peak, area under curve)
- `combined_plot.png` - Visualization showing ACTH, cortisol, and CRH

#### Analyze Stressor Results

After running multiple simulations with different parameters, analyze the results:

**Generate heatmaps** showing how stressor effects vary by timing and magnitude:
```bash
cd stressors/stressor_analysis
python heatmaps_stressors.py --scenario ACUTE
python heatmaps_stressors_heart_surgery.py --scenario HS1
```

**Create animations** showing dynamic response over time:
```bash
python make_animation_of_stressors.py --scenario ACUTE
```

**Find peaks and through timings**:
```bash
python get_stressor_times.py --scenario ACUTE --directory ../output/ACUTE/mag50_dur60_start14580/
```

All analysis scripts accept `--scenario` to target specific subdirectories (HS1, HS2, ACUTE, etc.). See `stressors/README.md` for full documentation.

## Data Format

### Input data format

Raw participant data files should be CSV with columns:
```
Date, Time, Cortisol, ACTH
```

Example:
```csv
Date,Time,Cortisol,ACTH
01/01/2021,09:00:00,15.2,45.8
01/01/2021,09:20:00,16.1,48.3
...
```

### Processed Data Format

After running `data_prep.py`, data is shifted to a common reference time:
```
x, y, test, timesteps_after_reftime
```

Where:
- `x`: Timestamp
- `y`: Hormone concentration
- `test`: Participant ID
- `timesteps_after_reftime`: Minutes after reference time (09:00)


## Visualisation Outputs

### Model Simulation Plots
- **Plot Type 1**: Separate ACTH and Cortisol time series
- **Plot Type 4**: Combined ACTCH + CORT

### Wavelet Analysis Plots
- Individual participant wavelet spectrograms
- Total power vs. period
- Mean/median power heatmaps

### Fitting Plots
- Observed vs. fitted time series
- Error evolution during optimization
- Parameter convergence plots
- Residual distributions

### Bifurcation Diagrams
- Equilibrium branches
- Limit cycles
- Two-parameter continuation

## Customization

### Adding New Parameters

1. Edit `model/config/parameters.json`
2. Update `model/code/classes/model_class.py` if new parameters require model changes

### Custom Error Functions

Create new error functions in `model/code/classes/error_function.py`:

```python
class CustomError(pints.ErrorMeasure):
    def __call__(self, parameters):
        pass
```

### Wavelet Analysis Options

Customize wavelet parameters in `data_analysis/options/participant_options.csv`:
```csv
participant_id,min_period,max_period,power_thresh,smoothing_wsize,sliding_window
1,60,360,0.0,20,100
...
```

## References

1. **Original Model**: Walker et al. (2010) - "Hypothalamic-pituitary-adrenal axis models"
2. **PINTS**: https://github.com/pints-team/pints
3. **DDE-BIFTOOL**: Engelborghs et al. (2007)
4. **PyBOAT**: Mönke et al. (2020) - "Oscillation detection in biological time series"

```

## Support

For questions or issues:
- Check existing issues on GitHub
- Contact: [bxl388@student.bham.ac.uk]

---

**Last Updated:** February 2026
# HPAAxisModel
