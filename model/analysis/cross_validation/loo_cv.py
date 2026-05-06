#!/usr/bin/env python3
"""
Leave-One-Out Cross-Validation (LOO-CV) for the base HPA DDE model
====================================================================

Implements the LOO-CV approach for addressing reviewer comments about
generalisation:

  For each of the 10 participants (left-out one at a time):
    1. Compute cohort-mean parameters from the remaining 9 participants
    2. Simulate the DDE model using those mean parameters
    3. Evaluate prediction error against the held-out participant's data

Error is reported in three flavours:
  (a) Raw normalised SSE   — no phase correction at all
  (b) Phase-aligned SSE    — best circular lag applied before scoring
                             (reuses FFT cross-correlation logic from
                              PhaseAlignedMSE in error_function.py)
  (c) t_s-fitted SSE       — structural parameters fixed at cohort mean;
                             only the phase-offset parameter t_s is
                             re-optimised on the test subject via 1-D
                             grid search. This is the recommended
                             primary metric: it separates "structural
                             model generalises" from "chronotype is
                             person-specific".

The in-sample errors (from the individual fitted configs) are loaded for
direct comparison.

t_s rationale
--------------
t_s shifts the entire circadian CRH driver 2π(t − t_s)/T_c.  It captures
individual chronotype (early bird vs night owl) and spans −43 to +120 min
across the cohort (σ ≈ 51 min).  Fitting only t_s on the test subject
takes ≈1441 DDE evaluations (fine grid search at 1-minute resolution over
±T_c around the cohort mean).

Usage (from HPAAxisModel/)
--------------------------
    python model/analysis/cross_validation/loo_cv.py
    python model/analysis/cross_validation/loo_cv.py --no_ts_fit
    python model/analysis/cross_validation/loo_cv.py --output_dir model/analysis/cross_validation/output
    python model/analysis/cross_validation/loo_cv.py --participant 3
"""

import sys
import os
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')           # non-interactive; safe on headless machines
import matplotlib.pyplot as plt
from pathlib import Path

# ── path setup ──────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
MODEL_DIR  = SCRIPT_DIR.parent.parent          # model/
REPO_ROOT  = MODEL_DIR.parent                  # HPAAxisModel/
CODE_DIR   = MODEL_DIR / 'code'                # model/code/

sys.path.insert(0, str(CODE_DIR))
from classes.model_class import Model          # noqa: E402

# ── constants ────────────────────────────────────────────────────────────────
ALL_P       = list(range(1, 11))               # participants 1–10
CONFIG_DIR  = MODEL_DIR / 'config' / 'base'
ACTH_CSV    = REPO_ROOT / 'data_analysis' / 'data' / 'data_shifted_ACTH_1min_09_00.csv'
CORT_CSV    = REPO_ROOT / 'data_analysis' / 'data' / 'data_shifted_Cortisol_1min_09_00.csv'

# Ordered list of all fitted parameter keys (must match config JSON keys)
FITTED_KEYS = [
    'gamma_a', 'k_c', 'm_a', 'gamma_c', 'alpha', 'm_c',
    'k_a', 'delay', 'lambda_s', 'lambda_a', 't_s', 'sigma',
]

# ── data / config helpers ────────────────────────────────────────────────────

def load_configs():
    """Return dict participant_int → config dict for all 10 participants."""
    configs = {}
    for p in ALL_P:
        path = CONFIG_DIR / f'parameters_{p}.json'
        with open(path) as f:
            configs[p] = json.load(f)
    return configs


def cohort_mean_params(configs, exclude_p):
    """
    Compute (params_dict, fixed_params_dict) where params_dict contains the
    arithmetic mean of each fitted parameter across all participants except
    exclude_p.
    """
    training_ps  = [p for p in ALL_P if p != exclude_p]
    param_arrays = {k: [] for k in FITTED_KEYS}
    fixed_params = None

    for p in training_ps:
        cfg    = configs[p]
        params = cfg['parameters']
        for k in FITTED_KEYS:
            param_arrays[k].append(float(params[k]))
        if fixed_params is None:
            fixed_params = dict(cfg.get('fixed_params', {'T_c': 1440, 'h': 1}))

    mean_params = {k: float(np.mean(param_arrays[k])) for k in FITTED_KEYS}
    return mean_params, fixed_params


def load_observed(acth_df, cort_df, participant_int):
    """Return (1440, 2) array of [ACTH, Cortisol] for the given participant."""
    acth = acth_df.loc[acth_df['test'] == float(participant_int), 'y'].values.astype(float)
    cort = cort_df.loc[cort_df['test'] == float(participant_int), 'y'].values.astype(float)
    assert len(acth) == 1440, f"Expected 1440 ACTH rows for P{participant_int}, got {len(acth)}"
    assert len(cort) == 1440, f"Expected 1440 Cortisol rows for P{participant_int}, got {len(cort)}"
    return np.column_stack([acth, cort])


# ── simulation ────────────────────────────────────────────────────────────────

def run_simulation(params_dict, fixed_params, num_days=8, days_to_keep=1):
    """
    Build a Model instance and run one simulation.

    The model runs for num_days days to spin up to steady-state and returns
    the last `days_to_keep` days.

    Returns np.ndarray of shape (1440 * days_to_keep, 2).
    """
    # Start simulation at 09:00 (540 min) to align with observed data
    times = np.arange(540, 1440 * num_days + 540, 1.0)
    model = Model(
        fixed_params=fixed_params,
        suggested_params=params_dict,
        signal='Both',
        num_days=num_days,
        reject=False,       # evaluation — do not penalise unusual peaks
        days_to_keep=days_to_keep,
    )
    init_vals = list(params_dict.values())
    return model.simulate(init_vals, times)   # (1440 * days_to_keep, 2)


# ── error metrics ─────────────────────────────────────────────────────────────

def _normalise_channel(sim_ch, obs_ch):
    """Normalise by max(max_sim, max_obs) — matches ErrorMeasureTwoSignal."""
    mv = max(float(np.max(sim_ch)), float(np.max(obs_ch)))
    if mv == 0:
        return np.zeros_like(sim_ch), np.zeros_like(obs_ch)
    return sim_ch / mv, obs_ch / mv


def normalised_sse(sim, obs):
    """
    Normalised SSE over both channels, matching the loss used during fitting:
        Σ_channels Σ_t  (sim_norm[t] − obs_norm[t])²
    sim, obs: (1440, 2)
    """
    total = 0.0
    for ch in range(sim.shape[1]):
        s_n, o_n = _normalise_channel(sim[:, ch], obs[:, ch])
        total += float(np.sum((s_n - o_n) ** 2))
    return total


def _best_circular_shift_both(sim, obs):
    """
    Find the integer circular lag k (minutes) that maximises combined
    cross-correlation across ACTH and Cortisol channels.

    Mirrors _best_circular_shift_both from PhaseAlignedMSE in error_function.py.
    """
    n = sim.shape[0]
    total_corr = None
    for ch in range(sim.shape[1]):
        s = sim[:, ch].astype(float)
        o = obs[:, ch].astype(float)
        s = s - s.mean();  s /= (np.max(np.abs(s)) + 1e-12)
        o = o - o.mean();  o /= (np.max(np.abs(o)) + 1e-12)
        corr_i = np.fft.irfft(np.fft.rfft(s) * np.conj(np.fft.rfft(o)), n=n)
        total_corr = corr_i if total_corr is None else total_corr + corr_i
    k = int(np.argmax(total_corr))
    if k > n // 2:
        k -= n    # signed lag in (-n/2, n/2]
    return k


def phase_aligned_sse(sim_2day, obs):
    """
    Circular-shift a 2-day simulation to maximally align with obs (1 day),
    then extract the middle 24-hour window to avoid edge artifacts.

    sim_2day: (2880, 2) array from 2-day simulation
    obs: (1440, 2) array of observed data

    Returns (sse, shift_minutes).
    """
    # Find best shift using middle day of 2-day window for alignment
    sim_middle_day = sim_2day[720:2160, :]

    k = _best_circular_shift_both(sim_middle_day, obs)

    # Apply shift to the full 2-day simulation
    sim_shifted_2day = np.roll(sim_2day, -k, axis=0)

    # Extract middle 24 hours (avoiding wrap-around edges)
    sim_shifted_middle = sim_shifted_2day[720:2160, :]

    return normalised_sse(sim_shifted_middle, obs), k


def optimize_ts(cohort_params, fixed_params, obs, num_days=8, log_file=None, grid_points=1441):
    """
    Fit only t_s on the test participant using cohort-mean structural parameters.
    Fine grid search with 1-minute resolution over ±T_c around cohort mean t_s.

    Using 1441 points provides 1-minute precision without the need for further
    refinement.

    Returns (sse, best_t_s, sim_at_best_t_s).
    """
    params_copy = dict(cohort_params)
    T_c = float(fixed_params.get('T_c', 1440))

    def log(msg):
        print(msg)
        if log_file:
            with open(log_file, 'a') as f:
                f.write(msg + '\n')

    def objective(ts_val):
        params_copy['t_s'] = float(ts_val)
        sim = run_simulation(params_copy, fixed_params, num_days=num_days)
        if np.any(sim >= 4999):
            return 1e9
        return normalised_sse(sim, obs)

    # Fine grid over full period with user-specified resolution
    # Cover [-T_c, +T_c] with grid_points points
    grid = np.linspace(-T_c, T_c, grid_points, endpoint=False)
    log(f"    Grid search: evaluating {len(grid)} t_s values...")
    grid_sse = []
    import time
    start_time = time.time()
    for i, ts in enumerate(grid):
        sse = objective(ts)
        grid_sse.append(sse)
        if (i + 1) % 15 == 0 or i == 0:
            progress = (i + 1) / len(grid) * 100
            current_best = min(grid_sse)
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta_seconds = (len(grid) - i - 1) / rate if rate > 0 else 0
            eta_minutes = eta_seconds / 60
            log(f"      {progress:.1f}% ({i+1}/{len(grid)}) | Best SSE: {current_best:.4f} | "
                f"Rate: {rate:.1f} sim/s | ETA: {eta_minutes:.1f} min")

    best_grid_idx = int(np.argmin(grid_sse))
    best_ts  = float(grid[best_grid_idx])
    sse_val  = grid_sse[best_grid_idx]
    log(f"    Grid search complete: best t_s = {best_ts:.1f} min  (SSE = {sse_val:.4f})")

    params_copy['t_s'] = best_ts
    sim_best = run_simulation(params_copy, fixed_params, num_days=num_days)
    return sse_val, best_ts, sim_best


# ── LOO-CV loop ───────────────────────────────────────────────────────────────

def run_single_fold(test_p, configs, acth_df, cort_df, do_ts_fit=True, num_days=8, output_dir=None, ts_grid_points=1441):
    """
    Run LOO-CV for a single test participant (for array job parallelization).
    Returns a single fold result dict.
    """
    log_file = None
    if output_dir:
        participant_dir = output_dir / f'participant_{test_p}'
        participant_dir.mkdir(parents=True, exist_ok=True)
        log_file = participant_dir / 'optimization_log.txt'
        with open(log_file, 'w') as f:
            f.write(f"LOO-CV Optimization Log for Participant {test_p}\n")
            f.write(f"{'='*60}\n\n")

    print(f"\n{'='*60}")
    print(f"  Running LOO-CV fold for participant {test_p} (TEST subject)")
    print(f"{'='*60}")

    mean_params, fixed_params = cohort_mean_params(configs, exclude_p=test_p)
    training_ts = [configs[p]['parameters']['t_s'] for p in ALL_P if p != test_p]
    print(f"  Cohort mean t_s : {mean_params['t_s']:.2f} min "
          f"(training range: {min(training_ts):.1f}–{max(training_ts):.1f} min)")
    print(f"  Individual fit t_s for P{test_p}: "
          f"{configs[test_p]['parameters']['t_s']:.2f} min")

    obs = load_observed(acth_df, cort_df, test_p)

    print(f"  Running cohort model simulation (num_days={num_days})...")
    sim_cohort_1day = run_simulation(mean_params, fixed_params, num_days=num_days, days_to_keep=1)
    sim_cohort_2day = run_simulation(mean_params, fixed_params, num_days=num_days, days_to_keep=2)

    if np.any(sim_cohort_1day >= 4999):
        print("  WARNING: cohort simulation rejected (unrealistic dynamics). "
              "Setting errors to NaN.")
        raw_sse_val = np.nan
        pa_sse_val  = np.nan
        pa_shift    = np.nan
    else:
        raw_sse_val = normalised_sse(sim_cohort_1day, obs)
        pa_sse_val, pa_shift = phase_aligned_sse(sim_cohort_2day, obs)
        print(f"  Raw SSE           : {raw_sse_val:.4f}")
        print(f"  Phase-aligned SSE : {pa_sse_val:.4f}  "
              f"(best circular shift = {pa_shift} min)")

    insample = configs[test_p].get('error_metrics', {}).get(
        'mse', configs[test_p].get('error', np.nan)
    )
    print(f"  Participant MSE (individual fit): {insample}")

    fold = {
        'participant'          : test_p,
        'cohort_mean_params'   : mean_params,
        'sim_cohort'           : sim_cohort_1day,
        'obs'                  : obs,
        'raw_sse'              : raw_sse_val,
        'phase_aligned_sse'    : pa_sse_val,
        'phase_shift_minutes'  : pa_shift,
        'insample_mse'         : float(insample) if insample is not None else np.nan,
        'test_t_s'             : float(configs[test_p]['parameters']['t_s']),
        'cohort_mean_t_s'      : float(mean_params['t_s']),
    }

    if not np.isnan(raw_sse_val):
        print(f"  Fitting t_s on test subject...")
        ts_sse_val, best_ts, sim_ts = optimize_ts(
            mean_params, fixed_params, obs, num_days=num_days,
            log_file=log_file, grid_points=ts_grid_points
        )
        fold['ts_fitted_sse']   = ts_sse_val
        fold['ts_fitted_value'] = best_ts
        fold['sim_ts_fitted']   = sim_ts
        print(f"  t_s-fitted SSE    : {ts_sse_val:.4f}  "
              f"(best t_s = {best_ts:.2f} min, "
              f"Δ from cohort mean = {best_ts - mean_params['t_s']:+.2f} min, "
              f"Δ from individual = {best_ts - configs[test_p]['parameters']['t_s']:+.2f} min)")
    else:
        fold['ts_fitted_sse']   = np.nan
        fold['ts_fitted_value'] = np.nan
        fold['sim_ts_fitted']   = None

    return fold


def run_loo_cv(configs, acth_df, cort_df, do_ts_fit=True, num_days=8, output_dir=None, no_plots=False, ts_grid_points=1441):
    """
    Run leave-one-out cross-validation over all 10 participants.
    If output_dir is provided, saves each fold incrementally.

    Returns a list of per-fold result dicts.
    """
    results = []

    for test_p in ALL_P:
        print(f"\n{'='*60}")
        print(f"  Fold {test_p}/10 — participant {test_p} is the TEST subject")
        print(f"{'='*60}")

        mean_params, fixed_params = cohort_mean_params(configs, exclude_p=test_p)
        training_ts = [configs[p]['parameters']['t_s'] for p in ALL_P if p != test_p]
        print(f"  Cohort mean t_s : {mean_params['t_s']:.2f} min "
              f"(training range: {min(training_ts):.1f}–{max(training_ts):.1f} min)")
        print(f"  Individual fit t_s for P{test_p}: "
              f"{configs[test_p]['parameters']['t_s']:.2f} min")

        obs = load_observed(acth_df, cort_df, test_p)

        print(f"  Running cohort model simulation (num_days={num_days})...")
        sim_cohort_1day = run_simulation(mean_params, fixed_params, num_days=num_days, days_to_keep=1)
        sim_cohort_2day = run_simulation(mean_params, fixed_params, num_days=num_days, days_to_keep=2)

        if np.any(sim_cohort_1day >= 4999):
            print("  WARNING: cohort simulation rejected (unrealistic dynamics). "
                  "Setting errors to NaN.")
            raw_sse_val = np.nan
            pa_sse_val  = np.nan
            pa_shift    = np.nan
        else:
            raw_sse_val = normalised_sse(sim_cohort_1day, obs)
            pa_sse_val, pa_shift = phase_aligned_sse(sim_cohort_2day, obs)
            print(f"  Raw SSE           : {raw_sse_val:.4f}")
            print(f"  Phase-aligned SSE : {pa_sse_val:.4f}  "
                  f"(best circular shift = {pa_shift} min)")

        insample = configs[test_p].get('error_metrics', {}).get(
            'mse', configs[test_p].get('error', np.nan)
        )
        print(f"  Participant MSE (individual fit): {insample}")

        fold = {
            'participant'          : test_p,
            'cohort_mean_params'   : mean_params,
            'sim_cohort'           : sim_cohort_1day,
            'obs'                  : obs,
            'raw_sse'              : raw_sse_val,
            'phase_aligned_sse'    : pa_sse_val,
            'phase_shift_minutes'  : pa_shift,
            'insample_mse'         : float(insample) if insample is not None else np.nan,
            'test_t_s'             : float(configs[test_p]['parameters']['t_s']),
            'cohort_mean_t_s'      : float(mean_params['t_s']),
        }

        if not np.isnan(raw_sse_val):
            print(f"  Fitting t_s on test subject...")
            ts_sse_val, best_ts, sim_ts = optimize_ts(
                mean_params, fixed_params, obs, num_days=num_days,
                log_file=None, grid_points=ts_grid_points
            )
            fold['ts_fitted_sse']   = ts_sse_val
            fold['ts_fitted_value'] = best_ts
            fold['sim_ts_fitted']   = sim_ts
            print(f"  t_s-fitted SSE    : {ts_sse_val:.4f}  "
                  f"(best t_s = {best_ts:.2f} min, "
                  f"Δ from cohort mean = {best_ts - mean_params['t_s']:+.2f} min, "
                  f"Δ from individual = {best_ts - configs[test_p]['parameters']['t_s']:+.2f} min)")
        else:
            fold['ts_fitted_sse']   = np.nan
            fold['ts_fitted_value'] = np.nan
            fold['sim_ts_fitted']   = None

        results.append(fold)

        # Save incrementally after each fold
        if output_dir is not None:
            save_individual_errors([fold], output_dir)
            if not no_plots:
                plot_individual_folds([fold], output_dir, do_ts_fit=do_ts_fit)

    return results


def save_individual_errors(results, output_dir):
    """Save individual participant error metrics and alignment parameters to their subdirectories."""
    for fold in results:
        p = fold['participant']
        participant_dir = output_dir / f'participant_{p}'
        participant_dir.mkdir(parents=True, exist_ok=True)

        error_data = {
            'participant'          : p,
            'raw_sse'              : fold['raw_sse'],
            'phase_aligned_sse'    : fold['phase_aligned_sse'],
            'phase_shift_minutes'  : fold['phase_shift_minutes'],
            'insample_mse'         : fold['insample_mse'],
            'test_t_s'             : fold['test_t_s'],
            'cohort_mean_t_s'      : fold['cohort_mean_t_s'],
        }
        if 'ts_fitted_sse' in fold:
            error_data['ts_fitted_sse']   = fold['ts_fitted_sse']
            error_data['ts_fitted_value'] = fold['ts_fitted_value']

        error_df   = pd.DataFrame([error_data])
        error_path = participant_dir / 'errors.csv'
        error_df.to_csv(error_path, index=False)
        print(f"  Saved error metrics: {error_path}")

        alignment_data = {
            'participant': p,
            'cohort_mean_parameters': {k: float(v) for k, v in fold['cohort_mean_params'].items()},
            'alignment_results': {
                'phase_shift_minutes' : float(fold['phase_shift_minutes']) if not np.isnan(fold['phase_shift_minutes']) else None,
                'ts_fitted_value'     : float(fold.get('ts_fitted_value', np.nan)) if not np.isnan(fold.get('ts_fitted_value', np.nan)) else None,
            },
            'error_metrics': {
                'raw_sse'           : float(fold['raw_sse']) if not np.isnan(fold['raw_sse']) else None,
                'phase_aligned_sse' : float(fold['phase_aligned_sse']) if not np.isnan(fold['phase_aligned_sse']) else None,
                'ts_fitted_sse'     : float(fold.get('ts_fitted_sse', np.nan)) if not np.isnan(fold.get('ts_fitted_sse', np.nan)) else None,
                'insample_mse'      : float(fold['insample_mse']) if not np.isnan(fold['insample_mse']) else None,
            },
            'recommended_method'              : 'ts_fitted' if 'ts_fitted_sse' in fold and not np.isnan(fold.get('ts_fitted_sse', np.nan)) else 'phase_aligned',
            'test_participant_original_t_s'   : float(fold['test_t_s']),
            'cohort_mean_t_s'                 : float(fold['cohort_mean_t_s']),
        }

        alignment_path = participant_dir / 'alignment_params.json'
        with open(alignment_path, 'w') as f:
            json.dump(alignment_data, f, indent=2)
        print(f"  Saved alignment params: {alignment_path}")


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_individual_folds(results, output_dir, do_ts_fit=True):
    """
    For each fold: 2-row × 1-col grid showing observed vs t_s-fitted cohort prediction.
    Row 0 = ACTH, row 1 = Cortisol.
    """
    t_min = np.arange(1440)
    t_hrs = t_min / 60.0

    for fold in results:
        p        = fold['participant']
        obs      = fold['obs']
        sim_ts   = fold.get('sim_ts_fitted')
        ts_val   = fold.get('ts_fitted_value', np.nan)
        ts_sse   = fold.get('ts_fitted_sse', np.nan)
        insample = fold['insample_mse']

        if sim_ts is None:
            print(f"  Skipping plot for P{p} (no t_s-fitted simulation available)")
            continue

        fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
        fig.suptitle(
            f'Leave-One-Out Cross-Validation — Participant {p} (Test Subject)\n'
            f'Participant MSE (individual fit) = {insample:.2f}  |  '
            f'Cohort model with optimized t_s: SSE = {ts_sse:.4f}',
            fontsize=11
        )

        ch_labels = ['ACTH (pg/mL)', 'Cortisol (nmol/L)']

        for row, ch in enumerate(range(2)):
            ax = axes[row]
            ax.plot(t_hrs, obs[:, ch], 'k-', lw=2, label='Observed', alpha=0.8)
            ax.plot(t_hrs, sim_ts[:, ch], color='forestgreen', ls='--', lw=2,
                    label=f'Cohort model (optimized t_s = {ts_val:.1f} min)')
            ax.set_ylabel(ch_labels[ch], fontsize=10)
            ax.legend(fontsize=9, loc='best')
            ax.grid(alpha=0.3)
            if row == 0:
                ax.set_title('Cohort Model with Optimized Delay Parameter (t_s)', fontsize=10)
            if row == 1:
                ax.set_xlabel('Time (hours)', fontsize=10)
                ax.set_xticks(np.arange(0, 25, 4))
                ax.set_xticklabels([f'{int(h):02d}:00' for h in np.arange(0, 25, 4)])

        plt.tight_layout()
        participant_dir = output_dir / f'participant_{p}'
        participant_dir.mkdir(parents=True, exist_ok=True)
        out_path = participant_dir / f'fold_P{p:02d}.png'
        plt.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {out_path}")


def plot_summary(results, output_dir, do_ts_fit=True):
    """
    Two-panel figure:
      Left  — per-participant grouped bar chart
      Right — box plot comparing error distributions
    """
    ps       = [r['participant']       for r in results]
    insample = [r['insample_mse']      for r in results]
    raw_sse  = [r['raw_sse']           for r in results]
    pa_sse   = [r['phase_aligned_sse'] for r in results]

    x     = np.arange(len(ps))
    width = 0.2 if do_ts_fit else 0.25

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle('LOO-CV: In-sample vs Out-of-sample Prediction Error', fontsize=13)

    ax1.bar(x - 1.5 * width, insample, width, label='In-sample (individual fit)',
            color='steelblue', alpha=0.85)
    ax1.bar(x - 0.5 * width, raw_sse,  width, label='LOO-CV raw SSE',
            color='tomato', alpha=0.85)
    ax1.bar(x + 0.5 * width, pa_sse,   width, label='LOO-CV phase-aligned SSE',
            color='seagreen', alpha=0.85)

    data_sets  = [insample, raw_sse, pa_sse]
    box_labels = ['In-sample\n(individual fit)', 'LOO-CV\n(raw)', 'LOO-CV\n(phase-aligned)']
    box_colors = ['steelblue', 'tomato', 'seagreen']

    if do_ts_fit:
        ts_sse = [r.get('ts_fitted_sse', np.nan) for r in results]
        ax1.bar(x + 1.5 * width, ts_sse, width, label='LOO-CV t_s-fitted SSE',
                color='darkorchid', alpha=0.85)
        data_sets.append(ts_sse)
        box_labels.append('LOO-CV\n(t_s-fitted)')
        box_colors.append('darkorchid')

    ax1.set_xticks(x)
    ax1.set_xticklabels([f'P{p}' for p in ps])
    ax1.set_xlabel('Participant (test subject in LOO-CV fold)')
    ax1.set_ylabel('Normalised SSE')
    ax1.set_title('Error per participant')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.grid(axis='y', alpha=0.3)

    bp = ax2.boxplot(data_sets, labels=box_labels, patch_artist=True, notch=False,
                     medianprops={'color': 'black', 'linewidth': 2})
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    for i, d in enumerate(data_sets):
        valid = [v for v in d if not np.isnan(v)]
        if valid:
            ax2.text(i + 1, max(valid) * 1.03,
                     f'μ={np.mean(valid):.3f}', ha='center', fontsize=8)

    ax2.set_ylabel('Normalised SSE')
    ax2.set_title('Distribution across all 10 participants')
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    out_path = output_dir / 'loo_cv_summary.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path.name}")


def plot_phase_analysis(results, output_dir, do_ts_fit=True):
    """
    Two-panel figure exploring the phase dimension:
      Left  — scatter: individual fitted t_s vs cohort-mean t_s
      Right — bar: circular shift found by phase-aligned SSE (and t_s-fitted
              value if available), showing how much phase correction is needed
              per test subject
    """
    ps        = [r['participant']     for r in results]
    test_ts   = [r['test_t_s']        for r in results]
    cohort_ts = [r['cohort_mean_t_s'] for r in results]
    pa_shifts = [r['phase_shift_minutes'] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle('Phase offset (t_s) analysis', fontsize=12)

    ax1.scatter(cohort_ts, test_ts, c='steelblue', s=90, zorder=3)
    for i, p in enumerate(ps):
        ax1.annotate(f'P{p}', (cohort_ts[i], test_ts[i]),
                     textcoords='offset points', xytext=(5, 4), fontsize=8)
    lim_lo = min(min(cohort_ts), min(test_ts)) - 10
    lim_hi = max(max(cohort_ts), max(test_ts)) + 10
    ax1.plot([lim_lo, lim_hi], [lim_lo, lim_hi], 'k--', alpha=0.4, label='y = x (perfect)')
    ax1.set_xlabel('Cohort-mean t_s (min)')
    ax1.set_ylabel('Individual fitted t_s (min)')
    ax1.set_title('Individual vs cohort-mean phase parameter t_s')
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.3)

    x = np.arange(len(ps))
    w = 0.35 if do_ts_fit else 0.6
    ax2.bar(x - w / 2 if do_ts_fit else x, pa_shifts, w,
            label='Circular shift (phase-aligned)', color='coral', alpha=0.85)

    if do_ts_fit:
        ts_fitted = [r.get('ts_fitted_value', np.nan) for r in results]
        ts_delta  = [tf - cm for tf, cm in zip(ts_fitted, cohort_ts)]
        ax2.bar(x + w / 2, ts_delta, w,
                label='t_s correction (t_s-fitted − cohort mean)', color='darkorchid', alpha=0.85)

    ax2.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'P{p}' for p in ps])
    ax2.set_xlabel('Participant (test subject)')
    ax2.set_ylabel('Phase correction (minutes)')
    ax2.set_title('Phase correction required to fit each test subject')
    ax2.legend(fontsize=8)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    out_path = output_dir / 'phase_analysis.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path.name}")


def plot_generalisation_ratio(results, output_dir, do_ts_fit=True):
    """
    Bar chart of generalisation ratio = LOO-CV error / in-sample error.
    Values close to 1.0 indicate good generalisation (no overfitting).
    """
    ps       = [r['participant']       for r in results]
    insample = np.array([r['insample_mse']      for r in results], dtype=float)
    raw_sse  = np.array([r['raw_sse']           for r in results], dtype=float)
    pa_sse   = np.array([r['phase_aligned_sse'] for r in results], dtype=float)

    ratio_raw = raw_sse / insample
    ratio_pa  = pa_sse  / insample

    x     = np.arange(len(ps))
    width = 0.25 if do_ts_fit else 0.35

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.bar(x - width, ratio_raw, width, label='Raw SSE / in-sample', color='tomato', alpha=0.85)
    ax.bar(x,         ratio_pa,  width, label='Phase-aligned SSE / in-sample',
           color='seagreen', alpha=0.85)

    if do_ts_fit:
        ts_sse   = np.array([r.get('ts_fitted_sse', np.nan) for r in results], dtype=float)
        ratio_ts = ts_sse / insample
        ax.bar(x + width, ratio_ts, width, label='t_s-fitted SSE / in-sample',
               color='darkorchid', alpha=0.85)

    ax.axhline(1.0, color='black', lw=1.2, ls='--', label='Ratio = 1 (no overfitting)')
    ax.set_xticks(x)
    ax.set_xticklabels([f'P{p}' for p in ps])
    ax.set_xlabel('Participant (test subject in LOO-CV fold)')
    ax.set_ylabel('LOO-CV error / in-sample error')
    ax.set_title('Generalisation ratio  (lower overfitting → ratio closer to 1.0)')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    out_path = output_dir / 'generalisation_ratio.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_path.name}")


# ── results I/O ───────────────────────────────────────────────────────────────

def save_csv(results, output_dir, do_ts_fit=True):
    """Save per-fold numeric results to a CSV file."""
    rows = []
    for r in results:
        row = {
            'participant'          : r['participant'],
            'insample_mse'         : r['insample_mse'],
            'loo_raw_sse'          : r['raw_sse'],
            'loo_phase_aligned_sse': r['phase_aligned_sse'],
            'phase_shift_minutes'  : r['phase_shift_minutes'],
            'test_t_s'             : r['test_t_s'],
            'cohort_mean_t_s'      : r['cohort_mean_t_s'],
        }
        if do_ts_fit:
            row['loo_ts_fitted_sse'] = r.get('ts_fitted_sse', np.nan)
            row['loo_ts_fitted_t_s'] = r.get('ts_fitted_value', np.nan)
        rows.append(row)
    df = pd.DataFrame(rows)
    out_path = output_dir / 'loo_cv_results.csv'
    df.to_csv(out_path, index=False, float_format='%.6f')
    print(f"Saved: {out_path.name}")
    return df


def print_summary_table(df, do_ts_fit=True):
    """Print formatted per-fold + aggregate summary to stdout."""
    print(f"\n{'='*72}")
    print("LOO-CV RESULTS SUMMARY")
    print(f"{'='*72}")

    display_cols = ['participant', 'insample_mse', 'loo_raw_sse',
                    'loo_phase_aligned_sse']
    if do_ts_fit:
        display_cols += ['loo_ts_fitted_sse', 'loo_ts_fitted_t_s']
    display_cols += ['phase_shift_minutes', 'test_t_s', 'cohort_mean_t_s']

    print(df[display_cols].to_string(index=False, float_format='{:.4f}'.format))
    print(f"\n{'-'*72}")
    print("AGGREGATED ACROSS ALL 10 FOLDS (mean ± std):")

    metric_cols = ['insample_mse', 'loo_raw_sse', 'loo_phase_aligned_sse']
    if do_ts_fit:
        metric_cols.append('loo_ts_fitted_sse')
    for col in metric_cols:
        vals = df[col].dropna().values
        if len(vals):
            print(f"  {col:<32}: {np.mean(vals):.4f} ± {np.std(vals):.4f}")

    print("\nGeneralisation ratios (LOO-CV error / in-sample MSE):")
    insample = df['insample_mse'].values
    for col, label in [
        ('loo_raw_sse',            'raw'),
        ('loo_phase_aligned_sse',  'phase-aligned'),
        ('loo_ts_fitted_sse',      't_s-fitted'),
    ]:
        if col not in df.columns:
            continue
        ratio = df[col].values / insample
        print(f"  {label:<20}: mean = {np.nanmean(ratio):.3f}  "
              f"(1.0 = perfect generalisation, "
              f"values >1 indicate increased error on unseen data)")
    print(f"{'='*72}")


# ── CLI entry point ───────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description='Leave-One-Out Cross-Validation for the base HPA DDE model'
    )
    parser.add_argument(
        '--no_ts_fit', action='store_true',
        help='Skip t_s optimisation on test subjects (faster; only raw + phase-aligned)'
    )
    parser.add_argument(
        '--output_dir', type=str, default=None,
        help='Directory to save results (default: cross_validation/output/)'
    )
    parser.add_argument(
        '--no_plots', action='store_true',
        help='Skip plot generation — only produce the CSV results file'
    )
    parser.add_argument(
        '--num_days', type=int, default=8,
        help='Spin-up days for DDE steady-state (default: 8, matching original fits)'
    )
    parser.add_argument(
        '--participant', type=int, default=None,
        help='Run LOO-CV for a single test participant (for array job parallelization)'
    )
    parser.add_argument(
        '--force_ts', type=float, default=None,
        help='Override t_s for all participants (for testing)'
    )
    parser.add_argument(
        '--ts_grid_points', type=int, default=1441,
        help='Number of grid search points for t_s optimization (default: 1441)'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    do_ts_fit          = not args.no_ts_fit
    no_plots           = args.no_plots
    num_days           = args.num_days
    single_participant = args.participant
    output_dir = (
        Path(args.output_dir) if args.output_dir
        else SCRIPT_DIR / 'output'
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory : {output_dir.resolve()}")
    print(f"t_s optimisation : {'enabled' if do_ts_fit else 'disabled (--no_ts_fit)'}")
    print(f"Spin-up days     : {num_days}")
    if single_participant:
        print(f"Single fold mode : participant {single_participant} (array job)")

    print("\nLoading participant configs...")
    configs = load_configs()
    if args.force_ts is not None:
        print(f"Overriding t_s for all participants to {args.force_ts} min")
        for p in configs:
            configs[p]['parameters']['t_s'] = args.force_ts
    print(f"Loaded {len(configs)} configs: participants {list(configs.keys())}")

    print("Loading observed data...")
    acth_df = pd.read_csv(ACTH_CSV)
    cort_df = pd.read_csv(CORT_CSV)
    print(f"  ACTH data : {len(acth_df)} rows  |  Cortisol data : {len(cort_df)} rows")

    if single_participant:
        fold = run_single_fold(
            single_participant, configs, acth_df, cort_df,
            do_ts_fit=do_ts_fit, num_days=num_days, output_dir=output_dir,
            ts_grid_points=args.ts_grid_points
        )
        results = [fold]

        save_individual_errors(results, output_dir)

        import pickle
        fold_file = output_dir / f'fold_{single_participant}.pkl'
        with open(fold_file, 'wb') as f:
            pickle.dump(fold, f)
        print(f"\nSaved fold data to {fold_file}")

        if not no_plots:
            print("\nGenerating fold plot...")
            plot_individual_folds(results, output_dir, do_ts_fit=do_ts_fit)
    else:
        results = run_loo_cv(
            configs, acth_df, cort_df, do_ts_fit=do_ts_fit, num_days=num_days,
            output_dir=output_dir, no_plots=no_plots, ts_grid_points=args.ts_grid_points
        )

        df = save_csv(results, output_dir, do_ts_fit=do_ts_fit)
        print_summary_table(df, do_ts_fit=do_ts_fit)

        if not no_plots:
            print("\nGenerating summary plots...")
            plot_summary(results, output_dir, do_ts_fit=do_ts_fit)
            plot_phase_analysis(results, output_dir, do_ts_fit=do_ts_fit)
            plot_generalisation_ratio(results, output_dir, do_ts_fit=do_ts_fit)
            print(f"All plots saved to {output_dir.resolve()}")

    print("\nDone.")


if __name__ == '__main__':
    main()
