# HPA Axis Mathematical Model

A comprehensive mathematical modeling framework for the Hypothalamic-Pituitary-Adrenal (HPA) axis using delay differential equations (DDEs). This repository contains everything needed to simulate, fit, analyze, and perform bifurcation analysis on the cortisol-ACTH feedback system.

## 📋 Overview

This repository provides a complete pipeline for HPA axis modeling:

1. **Data Preparation** - Process and prepare raw hormone time series data
2. **Wavelet Analysis** - Identify ultradian rhythms and periodic patterns
3. **Model Simulation** - Run the base HPA axis DDE model
4. **Parameter Fitting** - Fit model parameters to experimental data using PINTS
5. **Bifurcation Analysis** - Analyze system dynamics using DDE-BIFTOOL (MATLAB)

The model describes the feedback interactions between:
- **CRH** (Corticotropin-Releasing Hormone)
- **ACTH** (Adrenocorticotropic Hormone)  
- **Cortisol** (the primary stress hormone)

## Repository Structure

```
HPAAxisModel/
├── model/                      # Core HPA axis model
│   ├── code/
│   │   ├── run_model.py       # Main simulation script
│   │   ├── classes/           # Model class definitions
│   │   └── additional_functions/
│   ├── config/                # Parameter configuration files
│   ├── data/                  # Model input data
│   └── output/                # Simulation outputs
│
├── pints_fitting/             # Parameter optimization
│   └── run_pints.py          # PINTS fitting script
│
├── data_analysis/             # Data processing and analysis
│   ├── code/
│   │   ├── data_prep.py      # Data preparation pipeline
│   │   └── wavelet_analysis.py  # PyBOAT wavelet analysis
│   ├── data/                  # Raw and processed data
│   └── output/                # Analysis outputs
│
├── ddebiftool/                # Bifurcation analysis (MATLAB)
│   ├── run_pipeline/          # MATLAB scripts
│   ├── functions/             # Helper functions
│   ├── 3_draw_samples_for_matlab.py  # Config generation
│   └── 4_parameter_ranges.py
```

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
pyboat-dtk  # For wavelet analysis
imageio     # Optional, for animations
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
pip install numpy pandas matplotlib scipy seaborn pints ddeint pyboat-dtk imageio
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
    "k1": 2.5,
    "k2": 1.8,
    "k3": 3.2,
    "tau1": 15,
    "tau2": 45,
    "tau3": 180
  },
  "fixed_params": {
    "T": 1440,
    "num_days": 3
  },
  "signal": "Both"
}
```

### 4. Parameter Fitting

Fit model parameters to experimental data using PINTS optimization:

```bash
cd pints_fitting
python run_pints.py --config parameters.json --sensitivities 1
```

**Arguments:**
- `--config`: Configuration file with initial parameter guesses
- `--sensitivities`: Run sensitivity analysis (1 or 0)
- `--phase_align`: Align model to data by phase shifting (1 or 0)
- `--animation`: Generate optimization animation (1 or 0)

**Output:**
- `output/{participant}.{run}/` - Fitted parameters and plots
- `output/{participant}.{run}/parameters.txt` - Fitted parameter values
- `output/{participant}.{run}/sensitivities.npy` - Sensitivity matrices (if enabled)

**Fitting Process:**
1. Load observed ACTH and Cortisol data
2. Initialize model with parameter guesses from config
3. Run CMA-ES optimization to minimize error
4. Save best-fit parameters and generate comparison plots

### 5. Bifurcation Analysis

Analyze system dynamics and stability using DDE-BIFTOOL in MATLAB:

#### Step 1: Generate Configuration Files

```bash
cd ddebiftool
python 3_draw_samples_for_matlab.py
python 4_parameter_ranges.py
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

## 📐 Model Equations

The HPA axis model is described by a system of delay differential equations:

```
dCRH/dt = f1(CRH(t), Cortisol(t-τ1))
dACTH/dt = f2(ACTH(t), CRH(t-τ2), Cortisol(t))
dCortisol/dt = f3(Cortisol(t), ACTH(t-τ3))
```

Where:
- τ1, τ2, τ3 are time delays representing signal transmission times
- f1, f2, f3 are nonlinear functions describing hormone production and clearance
- Feedback inhibition is modeled through cortisol's effect on CRH and ACTH

Key features:
- **Negative feedback**: Cortisol inhibits CRH and ACTH production
- **Time delays**: Account for synthesis and transport times
- **Ultradian rhythms**: Emerge naturally from the delayed feedback structure
- **Nonlinear dynamics**: Hill functions model saturation effects

## 📈 Data Format

### Input Data Format

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

## 🔬 Analysis Workflows

### Complete Pipeline Example

```bash
# 1. Prepare data
cd data_analysis/code
python data_prep.py

# 2. Wavelet analysis
python wavelet_analysis.py --hormone Cortisol

# 3. Fit parameters to participant 9
cd ../../pints_fitting
python run_pints.py --config config_participant_9.json

# 4. Run simulation with fitted parameters
cd ../model/code
python run_model.py --config ../../pints_fitting/output/9.1/best_fit_params.json --plots 1,3,4

# 5. Bifurcation analysis in MATLAB
# (Open MATLAB and navigate to ddebiftool/run_pipeline)
```

### Sensitivity Analysis

```bash
cd pints_fitting
python run_pints.py --config parameters.json --sensitivities 1
```

This computes partial derivatives of model outputs with respect to parameters, helping identify which parameters most influence the dynamics.

## 📊 Visualization Outputs

### Model Simulation Plots
- **Plot Type 1**: Separate ACTH and Cortisol time series
- **Plot Type 3**: Combined ACTH and Cortisol on same axes
- **Plot Type 4**: Combined with confidence bounds
- **Plot Type 5**: Includes CRH dynamics

### Wavelet Analysis Plots
- Individual participant wavelet spectrograms
- Total power vs. period
- Mean/median power heatmaps
- Statistical comparison matrices
- Boxplots comparing participants

### Fitting Plots
- Observed vs. fitted time series
- Error evolution during optimization
- Parameter convergence plots
- Residual distributions

### Bifurcation Diagrams
- Equilibrium branches
- Limit cycles
- Period doubling routes to chaos
- Two-parameter continuation

## 🛠️ Customization

### Adding New Parameters

1. Edit `model/config/parameters.json`
2. Update `model/code/classes/model_class.py` if new parameters require model changes
3. Adjust parameter bounds in `model/code/classes/custom_boundaries.py`

### Custom Error Functions

Create new error functions in `model/code/classes/error_function.py`:

```python
class CustomError(pints.ErrorMeasure):
    def __call__(self, parameters):
        # Your custom error calculation
        pass
```

### Wavelet Analysis Options

Customize wavelet parameters in `data_analysis/options/participant_options.csv`:
```csv
participant_id,min_period,max_period,power_thresh,smoothing_wsize,sliding_window
1,60,360,0.0,20,100
...
```

## 📚 References

1. **Original Model**: Walker et al. (2010) - "Hypothalamic-pituitary-adrenal axis models"
2. **PINTS**: https://github.com/pints-team/pints
3. **DDE-BIFTOOL**: Engelborghs et al. (2007)
4. **PyBOAT**: Mönke et al. (2020) - "Oscillation detection in biological time series"

## 🤝 Contributing

This is a research codebase. For questions or collaboration:
- Open an issue for bugs or feature requests
- Submit pull requests with clear descriptions
- Contact: [Your contact information]

## 📄 License

[Specify license - typically MIT or GPL for academic code]

## 🙏 Acknowledgments

- Laboratory/Institution name
- Funding sources
- Collaborators

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@software{hpa_axis_model,
  author = {[Your Name]},
  title = {HPA Axis Mathematical Model},
  year = {2026},
  url = {[Repository URL]}
}
```

## 📞 Support

For questions or issues:
- Check existing issues on GitHub
- Consult the documentation in `docs/`
- Contact: [your.email@university.edu]

---

**Last Updated:** February 2026
**Version:** 1.0.0
# HPAAxisModel
