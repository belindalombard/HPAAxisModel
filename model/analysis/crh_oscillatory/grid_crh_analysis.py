"""
Parametric study of CRH oscillation effects.
Generates a grid of plots (3 periods × 3 amplitudes) showing constant vs.
oscillating CRH drive for Low, Medium, and High baseline CRH conditions.

Requires the constant CRH .npy files pre-computed by simulate_crh.py (constant
mode).  If they are absent they are computed on the fly and cached.

Outputs (in model/analysis/crh_oscillatory/output/figures/crh_parameter_sweep/):
  ampl_{amp}_per_{period}.png  – one figure per (period, amplitude) combination

Usage:
    # full parameter sweep
    python model/analysis/crh_oscillatory/grid_crh_analysis.py

    # Figure S8 only (periods 90 + 60 min, amplitude 0.3×, single combined figure)
    python model/analysis/crh_oscillatory/grid_crh_analysis.py --figure-s8
"""

import argparse
import math
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib import gridspec
from pathlib import Path

# ── path setup ─────────────────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT  = _SCRIPT_DIR.parents[2]          # HPAAxisModel/
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from model.code.classes.model_class import Model  # noqa: E402

# ── argument parsing ───────────────────────────────────────────────────────────
_parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
_parser.add_argument('--figure-s8', action='store_true',
                     help='Generate Figure S8 only (periods 90+60 min, amplitude 0.3×baseline)')
_args = _parser.parse_args()
FIGURE_S8_MODE = _args.figure_s8

plt.rcParams.update({'font.size': 18})

# ── parameter grid ─────────────────────────────────────────────────────────────
PERIODS             = [30, 60, 90]          # minutes
AMPLITUDE_FRACTIONS = [0.1, 0.3, 0.5]      # fraction of baseline

# (label, cached-npy name, lambda_a value)  — names match simulate_crh.py output
CRH_CONDITIONS = [
    ("Low",    "low_crh.npy",    0.1),
    ("Medium", "medium_crh.npy", 10.0),
    ("High",   "high_crh.npy",   10000.0),
]
CONDITIONS       = [c[0] for c in CRH_CONDITIONS]
CONST_CRH_VALUES = [c[2] for c in CRH_CONDITIONS]
_NB_NAME         = {c[2]: c[1] for c in CRH_CONDITIONS}  # baseline → filename

# ── simulation constants ───────────────────────────────────────────────────────
NUM_DAYS   = 40
TIMESTEPS  = 1440 * NUM_DAYS
TIMES      = np.linspace(0, TIMESTEPS, TIMESTEPS)
SIGNAL     = 'both'
REJECT     = False

# ── paths ──────────────────────────────────────────────────────────────────────
_SIM_CONST = _SCRIPT_DIR / 'output' / 'simulations' / 'const'
_OUT_DIR   = _SCRIPT_DIR / 'output' / 'figures' / 'crh_parameter_sweep'
_OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── load base configuration ────────────────────────────────────────────────────
config_path = _REPO_ROOT / 'model' / 'config' / 'base' / 'mean_parameters.json'
with open(config_path, 'r') as f:
    base_config = json.load(f)

fixed_params           = base_config.get('fixed_params', {})
parameters             = base_config.get('parameters', {})
parameters['lambda_s'] = 0


# ── oscillating CRH model subclass ────────────────────────────────────────────
class ModelOscillatingCRH(Model):
    """Model subclass that replaces the circadian CRH drive with an ultradian cosine."""

    def __init__(self, baseline_crh, osc_amplitude, osc_period, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.baseline_crh    = baseline_crh
        self.osc_amplitude   = osc_amplitude
        self.omega_ultradian = 2 * math.pi / osc_period

    def crh(self, t, t_s=None, lambda_a=None, lambda_s=None, sigma=None, T_c=1440):
        crh_value = self.baseline_crh + self.osc_amplitude * math.cos(self.omega_ultradian * t)
        if self.stressor_object is not None:
            if self.stressor_type == 'HS':
                crh_value += self.stressor_object.stressor_heart_surgery(t)
            else:
                crh_value += self.stressor_object.stressor(t)
        return crh_value


# ── helper: load (or compute + cache) constant CRH simulation ─────────────────
def _load_const(baseline: float) -> np.ndarray:
    path = _SIM_CONST / _NB_NAME[baseline]
    if path.exists():
        data = np.load(str(path))
    else:
        print(f"  Simulating constant CRH (baseline={baseline})...")
        params_c = dict(parameters)
        params_c['lambda_a'] = baseline
        m = Model(fixed_params=fixed_params, suggested_params=params_c,
                  signal=SIGNAL, num_days=1, reject=REJECT, length_model=TIMESTEPS)
        data = m.simulate(m.suggested_parameters(), TIMES)
        if data.ndim == 1:
            data = data[:, np.newaxis]
        _SIM_CONST.mkdir(parents=True, exist_ok=True)
        np.save(str(path), data)
    return data


# ── shared display constants ───────────────────────────────────────────────────
HOURS_SHOW   = 24
MINUTES_SHOW = HOURS_SHOW * 60
X_AXIS       = np.arange(MINUTES_SHOW)
_T0          = datetime.strptime('09:00:00', "%H:%M:%S")
X_TICKS      = list(range(180, MINUTES_SHOW, 360))
X_LABELS     = [(_T0 + timedelta(minutes=i)).strftime("%H:%M") for i in X_TICKS]

# ── CRH axis bounds for the full grid (envelope over all amplitude fractions) ──
crh_bounds_by_condition = {}
for baseline in CONST_CRH_VALUES:
    crh_min = min(baseline - baseline * f for f in AMPLITUDE_FRACTIONS)
    crh_max = max(baseline + baseline * f for f in AMPLITUDE_FRACTIONS)
    crh_bounds_by_condition[baseline] = (crh_min * 0.95, crh_max * 1.05)
    print(f"{baseline} AU CRH bounds: [{crh_min:.2e}, {crh_max:.2e}]")

print()

# ══════════════════════════════════════════════════════════════════════════════
# FULL PARAMETER SWEEP
# ══════════════════════════════════════════════════════════════════════════════
if not FIGURE_S8_MODE:
    total   = len(PERIODS) * len(AMPLITUDE_FRACTIONS)
    current = 0

    for period in PERIODS:
        omega = 2 * math.pi / period

        for amp_frac in AMPLITUDE_FRACTIONS:
            current += 1
            print(f"  [{current}/{total}] period={period} min, amplitude={amp_frac}")

            fig  = plt.figure(figsize=(20, 12))
            grid = gridspec.GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.5)

            for row, (cond, const_crh) in enumerate(zip(CONDITIONS, CONST_CRH_VALUES)):
                amplitude        = const_crh * amp_frac
                crh_min, crh_max = crh_bounds_by_condition[const_crh]

                const_data      = _load_const(const_crh)
                const_crh_trace = np.full(MINUTES_SHOW, const_crh)

                print(f"  Simulating oscillating CRH (baseline={const_crh}, amplitude={amplitude})...")
                m_osc = ModelOscillatingCRH(
                    baseline_crh=const_crh, osc_amplitude=amplitude, osc_period=period,
                    fixed_params=fixed_params, suggested_params=dict(parameters),
                    signal=SIGNAL, num_days=1, reject=REJECT, length_model=TIMESTEPS,
                )
                osc_data = m_osc.simulate(m_osc.suggested_parameters(), TIMES)
                if osc_data.ndim == 1:
                    osc_data = osc_data[:, np.newaxis]

                time_offset   = len(osc_data) - MINUTES_SHOW
                osc_crh_trace = const_crh + amplitude * np.cos(omega * (X_AXIS + time_offset))

                # left column: constant CRH
                ax1 = fig.add_subplot(grid[row, 0])
                ax2 = ax1.twinx()
                ax3 = ax1.twinx()
                ax3.spines['right'].set_position(('outward', 90))
                ax1.plot(X_AXIS, const_data[-MINUTES_SHOW:, 0], color='blue', lw=2, zorder=3)
                ax2.plot(X_AXIS, const_data[-MINUTES_SHOW:, 1], color='red',  lw=2, zorder=3)
                ax3.plot(X_AXIS, const_crh_trace, color='lightgrey', lw=2, alpha=0.7, zorder=1)
                ax1.set_ylabel('ACTH (pmol/L)', color='blue', fontsize=16)
                ax2.set_ylabel('CORT (nmol/L)', color='red',  fontsize=16)
                ax3.set_ylabel('CRH (AU)',       color='grey', fontsize=14)
                ax1.set_ylim([0, 15]); ax2.set_ylim([1, 1100]); ax3.set_ylim([crh_min, crh_max])
                ax3.tick_params(axis='y', labelcolor='grey', labelsize=12)
                if row == 0:
                    ax1.set_title('Constant CRH', fontsize=20, fontweight='bold')
                ax1.text(0.02, 0.95, f'{cond} CRH', transform=ax1.transAxes, fontsize=16,
                         verticalalignment='top',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

                # right column: oscillating CRH
                ax1o = fig.add_subplot(grid[row, 1])
                ax2o = ax1o.twinx()
                ax3o = ax1o.twinx()
                ax3o.spines['right'].set_position(('outward', 90))
                ax1o.plot(X_AXIS, osc_data[-MINUTES_SHOW:, 0], color='blue', lw=2, zorder=3)
                ax2o.plot(X_AXIS, osc_data[-MINUTES_SHOW:, 1], color='red',  lw=2, zorder=3)
                ax3o.plot(X_AXIS, osc_crh_trace, color='lightgrey', lw=2, alpha=0.7, zorder=1)
                ax1o.set_ylabel('ACTH (pmol/L)', color='blue', fontsize=16)
                ax2o.set_ylabel('CORT (nmol/L)', color='red',  fontsize=16)
                ax3o.set_ylabel('CRH (AU)',       color='grey', fontsize=14)
                ax1o.set_ylim([0, 15]); ax2o.set_ylim([1, 1100]); ax3o.set_ylim([crh_min, crh_max])
                ax3o.tick_params(axis='y', labelcolor='grey', labelsize=12)
                if row == 0:
                    ax1o.set_title(f'Oscillating CRH ({period} min period)',
                                   fontsize=20, fontweight='bold')
                ax1o.text(0.02, 0.95, f'{cond} CRH', transform=ax1o.transAxes, fontsize=16,
                          verticalalignment='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

                # x-axis on bottom row only
                if row < 2:
                    ax1.set_xticklabels([])
                    ax1o.set_xticklabels([])
                else:
                    for ax in [ax1, ax1o]:
                        ax.set_xticks(X_TICKS)
                        ax.set_xticklabels(X_LABELS, fontsize=12)
                        ax.set_xlabel('Time', fontsize=16)

            plt.suptitle(f'Amplitude: {amp_frac}\u00d7 baseline',
                         fontsize=22, fontweight='bold', y=0.995)
            out_fname = f'ampl_{amp_frac}_per_{period}.png'
            plt.savefig(str(_OUT_DIR / out_fname), bbox_inches='tight', dpi=300)
            plt.close()
            print(f"    → {out_fname}\n")

    print("=" * 60)
    print("✓ ALL PARAMETER COMBINATIONS COMPLETE!")
    print("=" * 60)
    print(f"Generated {total} plots in: {_OUT_DIR}/")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S8  (--figure-s8 flag)
# Periods 90 min (upper group) + 60 min (lower group), amplitude 0.3×baseline.
# Single combined 6-row × 2-col figure with a spacer row between groups.
# ══════════════════════════════════════════════════════════════════════════════
else:
    S8_PERIODS  = [90, 60]
    S8_AMP_FRAC = 0.3

    # tight CRH bounds for the fixed amplitude
    s8_crh_bounds = {
        baseline: ((baseline - baseline * S8_AMP_FRAC) * 0.95,
                   (baseline + baseline * S8_AMP_FRAC) * 1.05)
        for baseline in CONST_CRH_VALUES
    }

    fig = plt.figure(figsize=(20, 26))
    main_grid = gridspec.GridSpec(
        7, 2, figure=fig,
        hspace=0.45, wspace=0.65,
        height_ratios=[1, 1, 1, 0.15, 1, 1, 1],   # spacer at row 3
    )

    for p_idx, period in enumerate(S8_PERIODS):
        omega      = 2 * math.pi / period
        gs_row_off = 0 if p_idx == 0 else 4        # skip spacer row 3

        for c_idx, (cond, const_crh) in enumerate(zip(CONDITIONS, CONST_CRH_VALUES)):
            amplitude        = const_crh * S8_AMP_FRAC
            crh_min, crh_max = s8_crh_bounds[const_crh]
            row              = gs_row_off + c_idx
            is_first         = (c_idx == 0)
            is_last          = (c_idx == len(CONDITIONS) - 1)

            const_data      = _load_const(const_crh)
            const_crh_trace = np.full(MINUTES_SHOW, const_crh)

            print(f"  Running oscillating CRH  period={period}min  {cond} ...")
            m_osc = ModelOscillatingCRH(
                baseline_crh=const_crh, osc_amplitude=amplitude, osc_period=period,
                fixed_params=fixed_params, suggested_params=dict(parameters),
                signal=SIGNAL, num_days=1, reject=REJECT, length_model=TIMESTEPS,
            )
            osc_data = m_osc.simulate(m_osc.suggested_parameters(), TIMES)
            if osc_data.ndim == 1:
                osc_data = osc_data[:, np.newaxis]

            time_offset   = len(osc_data) - MINUTES_SHOW
            osc_crh_trace = const_crh + amplitude * np.cos(omega * (X_AXIS + time_offset))

            for col, (data, crh_trace, col_title) in enumerate([
                (const_data, const_crh_trace, 'Constant CRH'),
                (osc_data,   osc_crh_trace,   f'Oscillating CRH ({period} min period)'),
            ]):
                ax1 = fig.add_subplot(main_grid[row, col])
                ax2 = ax1.twinx()
                ax3 = ax1.twinx()
                ax3.spines['right'].set_position(('outward', 90))

                ax1.plot(X_AXIS, data[-MINUTES_SHOW:, 0], color='blue', lw=2, zorder=3)
                ax2.plot(X_AXIS, data[-MINUTES_SHOW:, 1], color='red',  lw=2, zorder=3)
                ax3.plot(X_AXIS, crh_trace, color='lightgrey', lw=2, alpha=0.7, zorder=1)

                ax1.set_ylabel('ACTH (pmol/L)', color='blue', fontsize=16)
                ax2.set_ylabel('CORT (nmol/L)', color='red',  fontsize=16)
                ax3.set_ylabel('CRH (AU)',       color='grey', fontsize=14)
                ax1.set_ylim([0, 15]); ax2.set_ylim([1, 1100]); ax3.set_ylim([crh_min, crh_max])
                ax3.tick_params(axis='y', labelcolor='grey', labelsize=12)

                if is_first:
                    ax1.set_title(col_title, fontsize=20, fontweight='bold')
                ax1.text(0.02, 0.95, f'{cond} CRH', transform=ax1.transAxes, fontsize=16,
                         verticalalignment='top',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

                if is_last:
                    ax1.set_xticks(X_TICKS)
                    ax1.set_xticklabels(X_LABELS, fontsize=12)
                    ax1.set_xlabel('Time', fontsize=16)
                else:
                    ax1.set_xticklabels([])

    plt.suptitle(f'Amplitude: {S8_AMP_FRAC}\u00d7 baseline',
                 fontsize=22, fontweight='bold', y=0.995)

    # period-group labels in left margin
    fig.text(0.01, 0.76, '90 min\nperiod', fontsize=15, fontweight='bold',
             ha='center', va='center', rotation=90, color='#333333')
    fig.text(0.01, 0.28, '60 min\nperiod', fontsize=15, fontweight='bold',
             ha='center', va='center', rotation=90, color='#333333')

    out_base = _SCRIPT_DIR / 'output' / 'figures'
    for ext in ('png', 'pdf'):
        out_path = out_base / f'figure_s8.{ext}'
        fig.savefig(str(out_path), bbox_inches='tight', dpi=300)
        print(f'Saved: {out_path}')
    plt.close(fig)
