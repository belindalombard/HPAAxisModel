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
- `load_signals()` - Load simulation results
- `calculate_metrics()` - Compute response metrics (peak, AUC, recovery time)

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
- `points_to_consider.json` - Timing analysis results
- `stressors/output/stressor_fig.pdf/png` - Visualization

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

1. **Run stressor simulations** for a parameter sweep:
   ```bash
   # Loop over magnitudes and durations
   for mag in 50 100 200; do
     for dur in 60 120 240; do
       python stressors/run_model_stressor.py --scenario ACUTE --magnitude $mag --duration $dur
     done
   done
   ```

2. **Aggregate results** into a CSV file (manual or custom script)

3. **Generate heatmaps** to visualize parameter effects:
   ```bash
   cd stressors/stressor_analysis
   python heatmaps_stressors.py --scenario ACUTE
   ```

4. **Create animations** to visualize response dynamics:
   ```bash
   python make_animation_of_stressors.py --scenario ACUTE
   ```

5. **Analyze specific scenarios**:
   ```bash
   # Analyze different HS scenarios
   python heatmaps_stressors_heart_surgery.py --scenario HS1
   python heatmaps_stressors_heart_surgery.py --scenario HS2
   python heatmaps_stressors_heart_surgery.py --scenario HS3
   ```
   ```bash
   cd stressors/stressor_analysis
   python heatmaps_stressors.py
   ```

4. **Analyze specific scenarios**:
   ```bash
   python get_stressor_times.py
   python make_animation_of_stressors.py
   ```

## Output Directory Structure

```
stressors/output/
├── ACUTE/                          # Acute stressor results
│   ├── mag100.00_dur120.00_.../    # Individual simulation
│   │   ├── simulated_values_*.npy
│   │   ├── crh*.npy
│   │   ├── metrics.csv
│   │   └── *.png
│   └── ...
├── HS1/                            # Heart surgery scenario 1
├── HS2/                            # Heart surgery scenario 2
├── HS3/                            # Heart surgery scenario 3
├── HS4/                            # Heart surgery scenario 4
├── heatmaps/                       # Heatmap visualizations
├── animations/                     # Animated visualizations
├── results.csv                     # Aggregated ACUTE results
└── results_heart_surgery.csv       # Aggregated HS results
```

## Response Metrics

The following metrics are computed for each simulation:

- **Peak response** - Maximum ACTH/CORT values
- **Time to peak** - Minutes from stressor onset to peak
- **AUC (Area Under Curve)** - Total response magnitude
- **Recovery time** - Time to return to baseline
- **Overshoot** - Peak relative to baseline
- **Circadian phase** - Phase of circadian rhythm at stressor onset
