# Stressor Analysis

Simulates the HPA axis model response to acute stressors and heart surgery scenarios, and analyzes the resulting dynamics.

## Prerequisites

- Python 3.x with numpy, pandas, matplotlib, scipy, seaborn
- Model fitted parameters in `model/config/parameters_*.json`
- `ddeint` package for delay differential equation solving

## Main Script

### `run_model_stressor.py`
Runs HPA axis model simulations with stressor perturbations and saves results.

**Usage:**
```bash
# Acute stressor
python stressors/run_model_stressor.py --config parameters_9.json --scenario ACUTE --magnitude 100 --duration 120 --time_in_scope 1440

# Heart surgery scenario 1
python stressors/run_model_stressor.py --config parameters_9.json --scenario HS1

# Heart surgery scenario 2 (earlier surgery time)
python stressors/run_model_stressor.py --config parameters_9.json --scenario HS2
```

**Key Arguments:**
- `--config` - Parameter config file from model/config/
- `--scenario` - Stressor scenario: ACUTE, HS1, HS2, HS3, HS4
- `--magnitude` - Stressor magnitude (override scenario default)
- `--duration` - Stressor rise duration in minutes
- `--plateau` - Plateau duration after rise
- `--time_in_scope` - Time of stressor onset (minutes from start)
- `--save` - Save plots (default: True)
- `--plot_days` - Number of days to plot (default: 2)

**Scenarios:**
- **ACUTE**: Short acute stressor (duration: 120 min, magnitude: 100)
- **HS1**: Heart surgery at 10:45 (duration: 20 min rise + 160 min plateau, magnitude: 1000)
- **HS2**: Heart surgery at 09:45 (duration: 20 min rise + 160 min plateau, magnitude: 1000)
- **HS3**: Long surgery (duration: 180 min, magnitude: 1000)
- **HS4**: Extended surgery with custom decay (6 hour rise, magnitude: 1000)

**Outputs:**
Results are saved in `stressors/output/<scenario>/<parameter_string>/`:
- `simulated_values_original.npy` - Baseline simulation (no stressor)
- `simulated_values_stressor.npy` - Simulation with stressor
- `crh.npy`, `crh_original.npy` - CRH values with/without stressor
- `crh_stressor_only.npy` - Isolated stressor contribution
- `metrics.csv` - Response metrics (peak, AUC, time to peak, etc.)
- `combined_plot.png` - Combined ACTH/CORT/CRH visualization
- `stressor_only.png` - Isolated stressor profile

## Analysis Scripts (stressor_analysis/)

### `stressor_analysis.py`
Core analysis functions for loading and analyzing stressor simulation results.

**Key Functions:**
- `get_stressor_directories()` - Find simulation directories by magnitude/duration/time
- `parse_directory_name()` - Extract parameters from directory names
- `get_signals_from_directory()` - Load `.npy` simulation arrays
- `plot_time_series()` - Plot CORT baseline vs. perturbed with AUC and peak annotations
- `calculate_auc_after_perturbation()` - Compute AUC difference post-stressor (Figure S6B/C)
- `calculate_peak_diff_after_perturbation()` - Compute first-peak difference post-stressor (Figure S6B/C)

**Two-step workflow to generate aggregated `results.csv`:**

1. Un-comment the first `'''…'''` block (lines ~409–494) and run the script to generate
   per-directory `metrics.json` files (each contains `difference_in_auc` and `peak_diff`
   computed over the 1440 minutes following the stressor onset, plus phase labels).

2. The active code (after the second `'''`) reads those `metrics.json` files and appends
   rows to `../output/ACUTE/results.csv`, which is then consumed by `heatmaps_stressors.py`.

> **Note:** `metrics.json` files are NOT included in this repo (too large). Run step 1
> using the HPC-generated `.npy` files to regenerate them.

### `stressor_analysis_heart_surgery.py`
Specialized analysis for heart surgery scenarios.

**Functions:**
- `get_stressor_directories()` - Find HS scenario directories
- `load_signals()` - Load HS simulation results
- Analysis functions for surgery-specific metrics

### `heatmaps_stressors.py`
Generate heatmaps showing response metrics across different stressor parameters.

**Command-line Arguments:**
- `--scenario` - Scenario subdirectory (default: ACUTE). Options: ACUTE, HS1, HS2, HS3, HS4
- `--results-file` - Path to results CSV (default: ../output/{scenario}/results.csv)
- `--output-dir` - Output directory for heatmaps (default: ../output/{scenario}/heatmaps)

**Usage:**
```bash
cd stressors/stressor_analysis
# Default (ACUTE scenario)
python heatmaps_stressors.py

# Specify scenario
python heatmaps_stressors.py --scenario HS1
python heatmaps_stressors.py --scenario HS2

# Custom paths
python heatmaps_stressors.py --results-file path/to/results.csv --output-dir output/my_heatmaps
```

### `heatmaps_stressors_heart_surgery.py`
Generate heatmaps for heart surgery scenario analysis.

**Command-line Arguments:**
- `--scenario` - Scenario subdirectory (default: HS). Options: HS, HS1, HS2, HS3, HS4
- `--results-file` - Path to results CSV (default: ../output/{scenario}/results_heart_surgery.csv)

**Usage:**
```bash
# Default (HS scenario)
python heatmaps_stressors_heart_surgery.py

# Specify scenario
python heatmaps_stressors_heart_surgery.py --scenario HS3
```

### `get_stressor_times.py`
Analyze optimal stressor timing based on circadian phase.

**Command-line Arguments:**
- `--directory` - Path to simulation directory
- `--scenario` - Scenario subdirectory (default: ACUTE)

**Usage:**
```bash
# Default directory
python get_stressor_times.py

# Specify directory
python get_stressor_times.py --directory ../output/ACUTE/mag50.00_dur20.00_plat0.00_decay0.00_beta0.00_start14616.00/
python get_stressor_times.py --scenario HS1 --directory ../output/HS1/mag100.00_dur30.00_start14580.00/
```

**Outputs:**
- `stressor_analysis/points_to_consider.json` - Ultradian timing points for each circadian cycle
  (canonical location; generated here and must be copied to the HPC before running the sweep)
- `stressors/output/stressor_fig.pdf/png` - Visualization (Figure S6A)

### `make_animation_of_stressors.py`
Create animated visualizations of stressor responses over time.

**Command-line Arguments:**
- `--scenario` - Scenario subdirectory (default: ACUTE). Options: ACUTE, HS1, HS2, HS3, HS4
- `--base-dir` - Base directory for simulations (default: ../output/{scenario})
- `--output-dir` - Output directory for animations (default: ../output/animations)

**Usage:**
```bash
cd stressors/stressor_analysis
# Default (ACUTE scenario)
python make_animation_of_stressors.py

# Specify scenario
python make_animation_of_stressors.py --scenario HS1
python make_animation_of_stressors.py --scenario HS4

# Custom directories
python make_animation_of_stressors.py --base-dir ../output/HS2 --output-dir ../output/my_animations
```

**Outputs:**
- `{output-dir}/output_animation_{scenario}.mp4`
- `{output-dir}/output_animation_{scenario}.gif`

## Typical Workflow

### Local — single stressor run (quick test or heart-surgery scenario)

```bash
# Single acute stressor
python stressors/run_model_stressor.py --config parameters_9.json --scenario ACUTE \
    --magnitude 100 --duration 120 --time_in_scope 1440

# Heart surgery scenarios
python stressors/run_model_stressor.py --config parameters_9.json --scenario HS1
```

### Full ACUTE parameter sweep (HPC)

The 961-run ACUTE sweep (magnitudes 50–250, durations 20–240 min, timings at 4
ultradian phases × 7+ circadian cycles) is too large to store in the repository
and must be run on an HPC cluster.

1. **Generate stressor timings locally:**
   ```bash
   cd stressors/stressor_analysis
   python get_stressor_times.py        # writes points_to_consider.json here
   ```

2. **Upload `points_to_consider.json` and the bash sweep script to the cluster:**
   ```bash
   scp stressors/stressor_analysis/points_to_consider.json <cluster>:ddemodel/
   scp PhD_repository/ddemodel/bash_files/run_stressors.sh  <cluster>:
   ```

3. **Run the sweep on the cluster:**
   ```bash
   # On the HPC:
   bash run_stressors.sh
   ```

4. **Copy results back:**
   ```bash
   rsync -avz <cluster>:ddemodel_old/stressor_output/ACUTE/ stressors/output/ACUTE/
   ```

5. **Generate per-run analysis and `results.csv`:**
   ```bash
   cd stressors/stressor_analysis
   # Un-comment the first triple-quoted block in stressor_analysis.py, then:
   python stressor_analysis.py        # writes per-dir metrics.json + results.csv
   ```

6. **Generate heatmaps** (Figures 5C/D and related):
   ```bash
   python heatmaps_stressors.py       # reads ../output/ACUTE/results.csv
   ```

7. **Create animations:**
   ```bash
   python make_animation_of_stressors.py --scenario ACUTE
   ```

### Heart surgery analysis

```bash
cd stressors/stressor_analysis
python heatmaps_stressors_heart_surgery.py --scenario HS1
python heatmaps_stressors_heart_surgery.py --scenario HS2
python heatmaps_stressors_heart_surgery.py --scenario HS3
```

## Output Directory Structure

The ACUTE bulk output (`mag*` subdirectories) is **gitignored** — 961 directories are too
large to push. `metrics.csv.example` is the only tracked file in the ACUTE directory; it
shows the format of the per-run concatenated HPC output.

```
stressors/output/
├── ACUTE/                              # Acute stressor results (bulk dirs gitignored)
│   ├── metrics.csv.example             # ← TRACKED: sample of raw HPC output format
│   ├── results.csv                     # ← generated locally by stressor_analysis.py
│   ├── heatmaps/                       # ← generated by heatmaps_stressors.py
│   └── mag100.00_dur120.00_.../        # ← GITIGNORED: per-run simulation data
│       ├── simulated_values_original.npy
│       ├── simulated_values_stressor.npy
│       ├── crh*.npy
│       ├── metrics.csv                 # raw per-run metrics from run_model_stressor.py
│       └── metrics.json               # computed AUC/peak diff (generated by stressor_analysis.py loop)
├── HS1/                                # Heart surgery scenario 1 (gitignored)
├── HS2/                                # Heart surgery scenario 2 (gitignored)
├── HS3/                                # Heart surgery scenario 3 (gitignored)
├── HS4/                                # Heart surgery scenario 4 (gitignored)
├── animations/                         # Animated visualizations (gitignored)
└── stressor_fig.{pdf,png}             # Figure S6A — stressor timing overview
```

## Response Metrics

The following metrics are computed for each simulation:

- **Peak response** - Maximum ACTH/CORT values
- **Time to peak** - Minutes from stressor onset to peak
- **AUC (Area Under Curve)** - Total response magnitude
- **Recovery time** - Time to return to baseline
- **Overshoot** - Peak relative to baseline
- **Circadian phase** - Phase of circadian rhythm at stressor onset
