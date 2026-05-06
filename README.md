# HPA Axis Model

A mathematical modelling and analysis framework for the Hypothalamic-Pituitary-Adrenal (HPA) axis, built on delay differential equations (DDEs). The repository covers the full pipeline: data preparation, wavelet and cosinor analysis, model simulation, parameter fitting via CMA-ES, leave-one-out cross-validation, bifurcation analysis (DDE-BIFTOOL), stressor perturbations, and oscillatory CRH investigations.

## Overview

The model describes the ACTH–Cortisol feedback loop driven by a CRH input function. The pipeline is organised into the following components:

1. **Data Preparation** — resample and time-align raw hormone CSVs
2. **Data Visualisation** — raw hormone time series and cross-correlation plots (manuscript Fig. 2A–2B)
3. **Wavelet Analysis** — detect ultradian rhythms with PyBOAT; wavelet spectrogram (manuscript Fig. 2C); sensitivity checks on window size and interpolation
4. **Cosinor Analysis** — test for statistically significant ultradian components using CosinorPy
5. **Model Simulation** — integrate the DDE model and generate publication plots (manuscript Fig. 3A, 3B, 4A)
6. **Parameter Fitting** — fit participant parameters with CMA-ES via PINTS
7. **Leave-One-Out Cross-Validation** — assess model generalisation across the cohort
8. **Bifurcation Analysis** — trace equilibria and limit cycles with DDE-BIFTOOL (MATLAB)
9. **Stressor Perturbations** — simulate responses to acute stress or surgical stressors (manuscript Fig. 5B, 5C, 5D)
10. **CRH Oscillatory Analysis** — explore how an oscillating CRH input alters model dynamics

## Getting Started

### Prerequisites

**Python Requirements:**
```
python >= 3.8
numpy
pandas
matplotlib
scipy
seaborn
pints
ddeint
pyboat
CosinorPy
imageio   # optional, for animations
```

**MATLAB Requirements** (for bifurcation analysis only):
- MATLAB R2019b or later (The current set up is for 2026)
- [DDE-BIFTOOL](https://sourceforge.net/projects/ddebiftool/)

### Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd HPAAxisModel
```

2. Install Python dependencies:
```bash
pip install numpy pandas matplotlib scipy seaborn pints ddeint pyboat CosinorPy imageio
```

3. For MATLAB bifurcation analysis, download and place DDE-BIFTOOL in a known location (see §7 below).

## Usage Guide

### 1. Data preparation

Process raw participant data and prepare for analysis:

```bash
python data_analysis/code/data_prep.py
```

**Input:** Raw CSV files in `data_analysis/data/raw/participant_*.csv`

**Output:**
- `data_analysis/data/data_resampled.csv` — resampled time series
- `data_analysis/data/data_shifted_Cortisol_1min_09_00.csv` — time-aligned cortisol data
- `data_analysis/data/data_shifted_ACTH_1min_09_00.csv` — time-aligned ACTH data

The timestamp reference and interpolation method can be adjusted directly in the script.

---

### 2. Data visualisation (manuscript Figure 2A & 2B)

Generate raw hormone time series plots and cross-correlation analysis:

```bash
python data_analysis/code/data_plots.py
```

**Input:** Raw CSV files in `data_analysis/data/raw/participant_*.csv`

**Manuscript figures:**
| Figure | Description | Output file |
|--------|-------------|-------------|
| **2A** | Raw ACTH & Cortisol time series (all participants, shared x-axis) | `data_analysis/output/data_plots/data_acth_cort_raw_unsmoothed.pdf` |
| **2B** | ACTH–Cortisol cross-correlation (time-lag) | `data_analysis/output/data_plots/timelag.pdf` |

Additional standalone ACTH panel saved to `data_analysis/output/data_plots/data_acth_raw_unsmoothed.pdf`.

---

### 3. Wavelet Analysis

#### 3a. Main wavelet spectral diagrams (manuscript Figure 2C)

Identify ultradian rhythms using PyBOAT:

```bash
python data_analysis/code/wavelet_analysis/wavelet_analysis.py --hormone Cortisol --dt 10
python data_analysis/code/wavelet_analysis/wavelet_analysis.py --hormone ACTH --dt 10
```

**Arguments:**
- `--hormone`: `Cortisol` or `ACTH`
- `--dt`: sampling interval in minutes (default: 10)
- `--data_folder`: data directory (default: `data_analysis/data`)
- `--results_folder`: output directory (default: `data_analysis/output/{hormone}/{dt}min`)

**Manuscript figure:**
| Figure | Description | Output file |
|--------|-------------|-------------|
| **2C** | Wavelet spectrogram for Participant 9 (Cortisol, dt=10 min) | `data_analysis/output/Cortisol/10min/wavelet_Participant_9_Cortisol.pdf` |

Per-participant options (window size, period range, power threshold) can be set in `data_analysis/options/participant_options.csv`.

#### 3b. Window-size sensitivity

Tests five different window sizes (480–960 min) plus a no-normalisation reference, to confirm that the chosen window does not drive the ultradian signal:

```bash
python data_analysis/code/wavelet_analysis/wavelet_window_sensitivity.py
```

Output saved to `data_analysis/code/wavelet_analysis/output/window_sensitivity/`.

#### 3c. Interpolation-resolution sensitivity

Tests dt = 1, 5, and 10 min to confirm results are not artefacts of resampling:

```bash
python data_analysis/code/wavelet_analysis/wavelet_interpolation_sensitivity.py
```

Output saved to `data_analysis/code/wavelet_analysis/output/interpolation_sensitivity/`.

---

### 4. Cosinor Analysis

Uses CosinorPy to test for statistically significant ultradian periodicities. The scripts sinc-detrend the data before fitting to remove circadian bias.

#### 4a. Ultradian periodogram

Computes F-statistic periodograms over the ultradian range (60–240 min) for each participant:

```bash
python data_analysis/code/wavelet_analysis/cosinor_ultradian.py
```

Output saved to `data_analysis/code/wavelet_analysis/output/cosinor_ultradian/`.

#### 4b. Best-fit model selection

Compares single-component (nc=1) vs. two-component (nc=2) cosinor fits via nested F-test (α = 1 × 10⁻⁵):

```bash
python data_analysis/code/wavelet_analysis/cosinor_best_models.py
```

Output saved to `data_analysis/code/wavelet_analysis/output/cosinor_best_models/`.

---

### 5. Model Simulation

Run the HPA axis DDE model with a specified parameter set:

```bash
python model/code/run_model.py --config model/config/base/mean_parameters.json --plots 1,2
```

**Arguments:**
- `--config`: path to parameter JSON (e.g. `model/config/base/mean_parameters.json`)
- `--plots`: comma-separated list of plot types (1 = separate signals, 2 = combined with CRH)
- `--output`: directory for saved figures
- `--no_show`: suppress interactive display (save only)

**Manuscript figure:**
| Figure | Description | Command |
|--------|-------------|---------|
| **4A** | Cohort mean simulation — ACTH, CORT & CRH drive (dual-axis) | `python model/code/run_model.py --config model/config/base/mean_parameters.json --plots 2 --output model/output/figure4` |

Output saved to `model/output/figure4/plot_2.png`. The `mean_parameters.json` config has `plot_crh: true` and `plot_option: 2`, which selects the combined dual-axis CRH plot.

To generate the full set of simulation comparison figures across all fitted participants:

```bash
# Step 1: run simulations for all participants
python model/analysis/generate_simulation_plots.py --simulate

# Step 2: generate plots
python model/analysis/generate_simulation_plots.py --plot

# Or both at once:
python model/analysis/generate_simulation_plots.py --simulate --plot

# Figure 3B: violin plots of fitted parameters
python model/analysis/unified_violin_plots.py --config_dir model/config/base/
```

**Manuscript figures:**
| Figure | Description | Output file |
|--------|-------------|-------------|
| **3A** | Simulated vs. observed ACTH & Cortisol for all 10 participants | `model/analysis/output/figures/figure_simulation_grid.pdf` |
| **3B** | Violin plots of fitted model parameters across participants | `model/analysis/output/violin_plots/` |

---

### 6. Parameter Fitting

Fit participant-level parameters using CMA-ES optimisation via PINTS:

```bash
python pints_fitting/run_pints.py --config model/config/base/mean_parameters.json --participant 1
```

**Arguments:**
- `--config`: initial parameter file
- `--participant`: participant ID
- `--phase_align`: align by phase shift instead of fitting `t_s` (1 or 0)
- `--animation`: generate optimisation animation (1 or 0)

**Output:** `pints_fitting/output/{participant}.{run}/` — fitted parameters, convergence plots, residuals.

Note: Intended for use on an HPC cluster when fitting all participants.

---

### 7. Leave-One-Out Cross-Validation

Assesses how well the cohort-mean parameters generalise to held-out participants. The pipeline runs in three steps:

```bash
# Step 1: fit parameters leaving each participant out in turn
python model/analysis/cross_validation/loo_cv.py

# Step 2: aggregate individual result files
python model/analysis/cross_validation/aggregate_loo_cv_results.py

# Step 3: generate summary grid figure
python model/analysis/cross_validation/plot_loo_cv_grid.py

# Step 4: generate supplementary figure (violin + per-participant bars)
python model/analysis/cross_validation/plot_supplement_figure.py
```

Output is written to `model/analysis/cross_validation/output/`.

---

### 8. Bifurcation Analysis

Trace equilibria, Hopf bifurcations, and limit cycles using DDE-BIFTOOL in MATLAB.

#### Manuscript figure

| Figure | Description | Script | Output |
|--------|-------------|--------|--------|
| **4C** | Bifurcation diagram — circadian period continuation (Hopf boundaries, limit cycles) | `run_circ` in MATLAB | `ddebiftool/run_pipeline/` output directory |

#### Step 1: Prepare configuration files

```bash
python ddebiftool/get_csv_file_from_configs.py
python ddebiftool/draw_samples_for_matlab.py
```

#### Step 2: Run MATLAB pipeline

```matlab
cd ddebiftool/run_pipeline
run_all_configs         % batch analysis for all configurations
run_bif_participants    % individual participant bifurcation analysis
run_circ                % circadian-period continuation
```

The pipeline is split into numbered `part_*.m` scripts in `ddebiftool/run_pipeline/pipeline_files/` and can be run step by step or via the batch scripts above. DDE-BIFTOOL must be on the MATLAB path; the setup script (`part_1_setup.m`) handles this via an environment variable on both macOS and Windows.

**Output:** bifurcation diagrams, stability boundaries, period-vs-parameter curves.

---

### 9. Stressor Perturbations

Simulate how the HPA axis responds to external stressors (acute stress, surgical stressor profiles). The model runs for one day to build up baseline dynamics, then applies a CRH pulse at the specified time.

#### Manuscript figures (Figure 5B)

Three heart-surgery CRH profiles, each showing CORT vs. baseline with the perturbed CRH drive overlaid:

| Figure | Scenario | Description | Surgery start | CRH profile |
|--------|----------|-------------|---------------|-------------|
| **5B (i)**  | HS1 | Rapid onset + plateau                    | 10:45 | Rise 20 min, plateau 160 min |
| **5B (ii)** | HS3 | Gradual rise, no plateau                 | 10:45 | Rise 180 min, no plateau     |
| **5B (iii)**| HS4 | Very slow rise, capped surgery duration  | 10:45 | Rise 360 min, surgery cap 180 min |

```bash
# Figure 5B — generate all three heart-surgery plots
python stressors/run_model_stressor.py --scenario HS1 --save
python stressors/run_model_stressor.py --scenario HS3 --save
python stressors/run_model_stressor.py --scenario HS4 --save
```

Outputs saved to `stressors/output/HS{N}/<param_string>/final_style.png`.

#### Manuscript figures (Figure 5C & 5D)

Phase-response heatmaps across a sweep of stressor timing and duration, with amplitude fixed at 150 AU. The ultradian phase axis runs 0 (nadir) → π/2 (rising) → π (peak) → 3π/2 (falling).

| Figure | Metric | Output file |
|--------|--------|-------------|
| **5C** | CORT peak surge (peak_diff) | `stressors/output/ACUTE/heatmaps/fixed_magnitude_150_peakdiff.pdf` |
| **5D** | CORT AUC change (difference_in_auc) | `stressors/output/ACUTE/heatmaps/fixed_magnitude_150_auc.pdf` |

**Prerequisites:** a `results.csv` sweep file must exist at `stressors/output/ACUTE/results.csv`. This is built by running `stressor_analysis.py` over a grid of magnitude/duration/timing combinations (intended for HPC cluster). The original sweep used magnitude=150 AU with a range of durations and ultradian-phase-aligned timing values. The file is not committed to the repository as it is large and HPC-generated. 

```bash
# Figure 5C & 5D — generate heatmaps from the pre-computed sweep results
cd stressors/stressor_analysis
python heatmaps_stressors.py --scenario ACUTE
```

This produces one PDF/PNG per unique magnitude and duration value in `results.csv`. The Figure 5C and 5D panels specifically correspond to the `fixed_magnitude_150_*` outputs.

All stressor parameters are pre-configured in the `SCENARIOS` dict in `run_model_stressor.py` and can be overridden per-run with `--magnitude`, `--duration`, `--plateau`, `--time_in_scope` etc.

#### General usage

```bash
python stressors/run_model_stressor.py \
    --scenario ACUTE \
    --magnitude 50 \
    --duration 60 \
    --start_time 14580
```

**Arguments:**
- `--scenario`: `ACUTE`, `HS1`, `HS2`, `HS3`, or `HS4`
- `--magnitude`: CRH pulse amplitude
- `--duration`: pulse duration (minutes)
- `--start_time`: onset time (minutes after 09:00)
- `--config`: parameter file (optional, defaults to mean parameters)

**Output:** `stressors/output/{scenario}/mag{X}_dur{Y}_start{Z}/`
- `simulated_values_original.npy` — baseline dynamics
- `simulated_values_stressor.npy` — stressor dynamics
- `crh_stressor.npy` — CRH profile
- `metrics.csv` — peak increase, time-to-peak, AUC
- `combined_plot.png` — ACTH/Cortisol/CRH comparison

#### Stressor analysis scripts

```bash
# Heatmaps over timing × magnitude
python stressors/stressor_analysis/heatmaps_stressors.py --scenario ACUTE
python stressors/stressor_analysis/heatmaps_stressors_heart_surgery.py --scenario HS1

# Time-of-peak analysis
python stressors/stressor_analysis/get_stressor_times.py --scenario ACUTE

# Animation of dynamic response
python stressors/stressor_analysis/make_animation_of_stressors.py --scenario ACUTE
```

---

### 10. CRH Oscillatory Analysis

Investigates how replacing the constant CRH baseline with a sinusoidally oscillating input changes model dynamics.

#### Manuscript figures

| Figure | Description | Script | Output |
|--------|-------------|--------|--------|
| 4B | Low / Medium / High constant CRH — last 24 h of 40-day simulation | `simulate_crh.py constant` → `plot_crh.py constant` | `model/analysis/crh_oscillatory/output/figures/constant_crh_conditions.pdf` |

#### Usage

```bash
# Step 1 — Run simulations (saves .npy arrays to output/simulations/)
python model/analysis/crh_oscillatory/simulate_crh.py constant
python model/analysis/crh_oscillatory/simulate_crh.py oscillating 90 0.3

# Step 2 — Plot results
# Figure 4B: constant CRH (High / Medium / Low)
python model/analysis/crh_oscillatory/plot_crh.py constant

# Single oscillating simulation
python model/analysis/crh_oscillatory/plot_crh.py single

# Side-by-side constant vs oscillating comparison
python model/analysis/crh_oscillatory/plot_crh.py compare 60 0.5

# Parametric grid: 3 periods × 3 amplitudes
python model/analysis/crh_oscillatory/grid_crh_analysis.py

# Standalone CRH drive curve (for figures)
python model/analysis/crh_oscillatory/crh_plot_only.py
```

Output is written to `model/analysis/crh_oscillatory/output/`.

---

## Data Format

### Input data

Raw participant data files should be CSV with columns:

```
Date, Time, Cortisol, ACTH
```

### Processed data

After running `data_prep.py`, data is shifted to a common reference time and stored as:

```
x, y, test, timesteps_after_reftime
```

- `x` — timestamp
- `y` — hormone concentration
- `test` — participant ID
- `timesteps_after_reftime` — minutes after 09:00 reference

---

## Wavelet Analysis Options

Per-participant wavelet settings are stored in `data_analysis/options/participant_options.csv`:

```csv
participant_id,min_period,max_period,power_thresh,smoothing_wsize,sliding_window
1,60,360,0.0,20,100
...
```

---

## References

1. Walker et al. (2010) — "Hypothalamic-pituitary-adrenal axis models"
2. Pints — https://github.com/pints-team/pints
3. Engelborghs et al. (2007) — DDE-BIFTOOL
4. Mönke et al. (2020) — PyBOAT: oscillation detection in biological time series
5. CosinorPy — https://github.com/mmoskon/CosinorPy

---

**Last updated:** May 2026

