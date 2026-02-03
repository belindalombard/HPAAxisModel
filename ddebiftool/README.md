# DDE-BIFTOOL Pipeline

Bifurcation analysis tools for the HPA axis delay differential equation (DDE) model using DDE-BIFTOOL (MATLAB).

## Prerequisites

- MATLAB with DDE-BIFTOOL installed
- Python 3.x with pandas, numpy, openpyxl
- Fitted model parameters in `model/config/parameters_*.json`

## Python Scripts

### 1. `get_csv_file_from_configs.py`
Extracts parameter values from participant config JSON files and creates a consolidated CSV/Excel file.

**Usage:**
```bash
python ddebiftool/get_csv_file_from_configs.py
```

### 2. `draw_samples_for_matlab.py`
Generates MATLAB config files for bifurcation analysis from parameter statistics (mean, median, percentiles, per-participant).

**Usage:**
```bash
python ddebiftool/draw_samples_for_matlab.py
```

**Outputs:**
- `run_pipeline/config/matlab/mean_median/` - Mean and median parameter configs
- `run_pipeline/config/matlab/percentiles/` - Percentile-based configs
- `run_pipeline/config/participants/` - Individual participant simulation configs
- `run_pipeline/config/participants_bif/` - Individual participant bifurcation configs

## MATLAB Scripts (run_pipeline/)

### Main Pipelines

- **`run_steps.m`** - Run bifurcation analysis for a single config file
- **`run_steps_two_delays.m`** - Bifurcation analysis with two delay parameters
- **`run_circ.m`** - Circadian rhythm bifurcation analysis
- **`run_all_configs.m`** - Batch process all configs in `config/matlab/`
- **`simulate_all_participants.m`** - Simulate all participant models
- **`run_bif_participants.m`** - Run bifurcation analysis for all participants

### Configuration

Edit `run_pipeline/config/config.json` to set:
- Initial parameter values
- Bifurcation parameters
- Branch bounds and step sizes
- Initial values for ACTH/CORT

## Typical Workflow

1. Fit model parameters using `pints_fitting/run_pints.py` → saves to `model/config/parameters_*.json`
2. Extract parameters: `python ddebiftool/get_csv_file_from_configs.py`
3. Generate MATLAB configs: `python ddebiftool/draw_samples_for_matlab.py`
4. Run bifurcation analysis in MATLAB: `run_all_configs` or `run_bif_participants`
5. Results are saved in `run_pipeline/output/`
