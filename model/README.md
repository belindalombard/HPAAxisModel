# HPA Axis Model

This directory contains the delay differential equation (DDE) model of the hypothalamic-pituitary-adrenal (HPA) axis, describing the feedback interactions between CRH, ACTH, and cortisol.

## Model Overview

The model captures the key dynamics of the HPA axis.

The model uses delay differential equations to simulate hormone concentrations over time, producing realistic ultradian pulses superimposed on a circadian rhythm.


## Running Simulations

### Basic Simulation

Run the model with default or custom parameters:

```bash
cd model/code
python run_model.py --config ../config/base/parameters.json
```

**Arguments:**
- `--config` - Path to parameter configuration JSON file
- `--plots` - Comma-separated plot types (1=Separate ACTH/CORT, 3=Combined, 4=Combined+Bounds, 5=CRH)
- `--output` - Output directory to save plots (default: plots displayed but not saved)
- `--no_show` - Skip displaying plots

**Example:**
```bash
# Generate multiple plot types and save
python run_model.py --config ../config/base/parameters.json --plots 1,3,4 --output ../output/sim1

# Just show plots without saving
python run_model.py --plots 1
```

### Configuration File Format

Parameter configuration files are JSON format:

```json
{
  "parameters": {
    "p1": 2.5,
    "p2": 1.8,
    "p3": 3.2,
    "t_s": 540,
    "lambda_a": 4.5,
    "lambda_s": 0.8,
    "sigma": 2.1
  },
  "fixed_params": {
    "T": 1440
  },
  "num_days": 3,
  "signal": "Both"
}
```

The config file for each participant is available in model/config/base 

## Model Class

The core model is in `code/classes/model_class.py` and implements the PINTS `ForwardModel` interface for parameter fitting.

## Analysis Scripts

### Generate Plots

Create imulation plots for all participants:

```bash
cd model/analysis
python generate_simulation_plots.py --simulate --plot
```

**Flags:**
- `--simulate` - Run simulations for all participants
- `--plot` - Generate plots from saved simulations

**Outputs:**
- Individual participant plots (ACTH and cortisol vs. data)
- Summary plots showing all participants
- Saved in `model/analysis/output/`

**Color scheme (manuscript standard):**
- Blue solid: Simulated ACTH
- Blue dashed: Observed ACTH
- Red solid: Simulated cortisol
- Red dashed: Observed cortisol

### Violin Plots

Generate parameter distribution plots:

```bash
python unified_violin_plots.py
```

Shows fitted parameter distributions across participants.

## Plot Types

The model can generate various plot types (specified with `--plots` argument):

- **Type 1**: Separate subplots for ACTH and cortisol
- **Type 2**: Combined plot.


## Key Functions

**model_class.py:**
- `Model.simulate(parameters)` - Run DDE simulation
- `Model.crh(t)` - CRH drive function
- `Model.dde_system()` - DDE system equations

**additional_functions/model_plotting_functions_modular.py:**
- `plot_simulation_separate()` - Separate ACTH/CORT plots
- `plot_simulation_combined()` - Combined plot
- `plot_with_bounds()` - Plot with parameter-dependent bounds

## Notes

- Default simulation length is 3 days (4320 minutes)
- CRH drive combines circadian (24h) and ultradian components
