# Data Analysis

This directory contains scripts and data for analyzing HPA axis hormone time series from healthy participants.

## Data Source

The experimental data used in this analysis comes from:

This dataset provides high-frequency ACTH and cortisol measurements from healthy individuals, enabling detailed analysis of ultradian rhythms in the HPA axis.

## Analysis Methods

The wavelet analysis methods used to identify ultradian rhythms are based on:


## Data Preparation

The `data_prep.py` script processes raw participant data:

1. Loads CSV files from `data/raw/participant_*.csv`
2. Resamples to consistent time intervals (1 min or 10 min)
3. Shifts data to common reference times (00:00 or 09:00)
4. Saves processed files for model fitting and wavelet analysis

**Input format:** Raw CSV files with columns: `Date, Time, Cortisol, ACTH`

**Output files:**
- `data_shifted_Cortisol_{dt}min_{ref_time}.csv` - Cortisol aligned to reference time
- `data_shifted_ACTH_{dt}min_{ref_time}.csv` - ACTH aligned to reference time
- `data_resampled.csv` - Combined resampled data
- `data_per_participant_pyboat.csv` - Formatted for PyBOAT analysis

**Usage:**
```bash
cd data_analysis/code
python data_prep.py
```

## Wavelet Analysis

The `wavelet_analysis.py` script uses PyBOAT to identify ultradian rhythms:

1. Applies continuous wavelet transform to hormone time series
2. Identifies significant periodicities (typically 60-360 minutes)
3. Generates power spectra and statistical comparisons
4. Saves wavelet data and visualizations

**Usage:**
```bash
cd data_analysis/code
python wavelet_analysis.py --hormone Cortisol --dt 10
python wavelet_analysis.py --hormone ACTH --dt 10
```

**Arguments:**
- `--hormone` - Hormone to analyze (Cortisol or ACTH)
- `--dt` - Sampling interval in minutes (default: 10)
- `--data_folder` - Data directory (default: ../data)
- `--results_folder` - Output directory (default: ../output)
- `--ref_time` - Reference time for data alignment (default: 09_00)

**Output:**
- `output/{hormone}/{dt}min/` - Wavelet spectrograms and power plots
- `output/stacked_wavelets_{hormone}.npy` - Saved wavelet data

## Data Visualization

The `data_plots.py` script generates publication-quality plots of the hormone data.

**Usage:**
```bash
cd data_analysis/code
python data_plots.py
```

**Output:** Plots saved to `output/data_plots/`

## Configuration

Wavelet analysis parameters for individual participants can be customized in:
```
options/participant_options.csv
```

This file specifies period ranges, power thresholds, and smoothing windows for each participant's analysis.

## References

1. Henley, D. E., Leendertz, J. A., Russell, G. M., Wood, S. A., Taheri, S., Woltersdorf, W. W., & Lightman, S. L. (2009). Development of an automated blood sampling system for use in humans. *Journal of Medical Engineering & Technology*, 33(3), 199-208.

2. Upton, D. J., & Zavala, E. (2022). High-frequency hormone pulsatility in the healthy HPA axis: A review of findings and methodological considerations. *Psychoneuroendocrinology*, 142, 105786.

3. Mönke, G., Sorgenfrei, F. A., Schmal, C., & Granada, A. E. (2020). Optimal time frequency analysis for biological data - pyBOAT. *bioRxiv*, 2020.04.29.067744.
