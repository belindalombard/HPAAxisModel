#!/usr/bin/env python3
"""
LOO-CV grid figure: cohort predictions vs observed data (all 10 participants) — Figure S13.

Representative LOO-CV predictions showing observed (solid) and predicted (dashed)
ACTH and cortisol profiles using cohort parameters.
    python model/analysis/cross_validation/plot_loo_cv_grid.py

Pre-requisites: run loo_cv.py first so that fold_*.pkl files exist in
    model/analysis/cross_validation/output/

Layout: 5 rows x 2 columns, each cell showing ACTH (top) and Cortisol (bottom).
Prediction is the t_s-fitted LOO simulation if available, otherwise the raw
cohort-mean simulation.

Usage (from HPAAxisModel/):
    python model/analysis/cross_validation/plot_loo_cv_grid.py
    python model/analysis/cross_validation/plot_loo_cv_grid.py --output_dir PATH
"""

import argparse
import pickle
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

plt.rcParams.update({'font.size': 24, 'xtick.labelsize': 18, 'ytick.labelsize': 18})

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _parse_time_of_day(df, participant_number):
    """Return (time_strings_shifted, shift_index) so that 09:00 is at index 0."""
    sub = df[df['test'] == float(participant_number)].copy()
    if 'timesteps_after_reftime' in sub.columns:
        sub = sub.sort_values('timesteps_after_reftime')
    times = sub['x'].values[-1440:]
    tod = [datetime.strptime(t, '%Y-%m-%d %H:%M:%S').time() for t in times]
    tod_str = [f"{t.hour:02d}:{t.minute:02d}" for t in tod]
    try:
        idx_9 = tod_str.index('09:00')
    except ValueError:
        idx_9 = next((i for i, s in enumerate(tod_str) if s >= '09:00'), 0)
    shifted = tod_str[idx_9:] + tod_str[:idx_9]
    return shifted, idx_9


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def make_grid_figure(output_dir, data_dir):
    output_dir = Path(output_dir)
    data_dir   = Path(data_dir)

    acth_df = pd.read_csv(data_dir / 'data_shifted_ACTH_1min_00_00.csv')
    cort_df = pd.read_csv(data_dir / 'data_shifted_Cortisol_1min_00_00.csv')

    fig = plt.figure(figsize=(13, 16))
    outer = gridspec.GridSpec(5, 2, figure=fig, hspace=0.3, wspace=0.4)

    handles, labels = [], []

    for pn in range(1, 11):
        fold_path = output_dir / f'fold_{pn}.pkl'
        if not fold_path.exists():
            print(f"  Warning: fold_{pn}.pkl not found — skipping participant {pn}")
            continue

        with open(fold_path, 'rb') as fh:
            fold = pickle.load(fh)

        # Prefer t_s-fitted simulation; fall back to raw cohort
        if fold.get('sim_ts_fitted') is not None:
            sim = fold['sim_ts_fitted']
        elif fold.get('sim_cohort') is not None:
            sim = fold['sim_cohort']
        else:
            print(f"  Warning: no simulation in fold_{pn}.pkl — skipping")
            continue

        obs = fold['obs']

        # Clip simulation to last 1440 time steps (24 hours)
        if sim.shape[0] > 1440:
            sim = sim[-1440:, :]

        # Align so 09:00 is at x=0
        tod_str, idx_9 = _parse_time_of_day(acth_df, pn)
        obs_sh = np.concatenate([obs[idx_9:], obs[:idx_9]], axis=0)
        sim_sh = np.concatenate([sim[idx_9:], sim[:idx_9]], axis=0)

        row = (pn - 1) % 5
        col = 0 if pn <= 5 else 1

        inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                    subplot_spec=outer[row, col], hspace=0.3)

        # ACTH
        ax_a = fig.add_subplot(inner[0])
        l1, = ax_a.plot(tod_str, sim_sh[:, 0], '-', color='blue',
                        linewidth=0.7, label='LOO-CV')
        l2, = ax_a.plot(tod_str, obs_sh[:, 0], ':', color='blue',
                        linewidth=1.0, label='Observed')
        ax_a.set_ylabel('ACTH (pmol/L)', fontsize=10)
        ax_a.yaxis.set_label_coords(-0.12, 0.5)
        ax_a.tick_params(axis='y', labelsize=10)
        ax_a.set_title(f'Participant {pn}', fontsize=12, pad=10)

        # Cortisol
        ax_c = fig.add_subplot(inner[1])
        l3, = ax_c.plot(tod_str, sim_sh[:, 1], '-', color='red',
                        linewidth=0.7, label='LOO-CV')
        l4, = ax_c.plot(tod_str, obs_sh[:, 1], ':', color='red',
                        linewidth=1.0, label='Observed')
        ax_c.set_ylabel('CORT (nmol/L)', fontsize=10)
        ax_c.yaxis.set_label_coords(-0.12, 0.5)
        ax_c.tick_params(axis='y', labelsize=10)

        # x ticks every 3 hours
        tick_idx = list(range(0, len(tod_str), 180))
        ax_a.set_xticks(tick_idx)
        ax_a.tick_params(axis='x', labelbottom=False)
        ax_c.set_xticks(tick_idx)
        if row == 4:
            ax_c.set_xticklabels([tod_str[i] for i in tick_idx], fontsize=10)
            ax_c.set_xlabel('Time (hour)', fontsize=10)
        else:
            ax_c.set_xticklabels([])

        if pn == 10:
            handles = [l1, l2, l3, l4]
            labels  = ['LOO-CV ACTH', 'Observed ACTH',
                       'LOO-CV CORT', 'Observed CORT']

    fig.legend(handles, labels, loc='upper center', ncol=4,
               fontsize=12, frameon=False, bbox_to_anchor=(0.5, 0.99))
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])

    for suffix in ('pdf', 'png'):
        out = output_dir / f'loo_cv_predictions_grid.{suffix}'
        fig.savefig(out)
        print(f"Saved: {out}")

    plt.close()
    print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description='Grid figure: LOO-CV predictions vs observed (all participants)'
    )
    parser.add_argument('--output_dir', default=None,
                        help='Directory containing fold_*.pkl files')
    parser.add_argument('--data_dir', default=None,
                        help='Directory containing shifted data CSVs')
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    repo_root  = script_dir.parents[2]          # HPAAxisModel/

    output_dir = Path(args.output_dir) if args.output_dir \
                 else script_dir / 'output'
    data_dir   = Path(args.data_dir) if args.data_dir \
                 else repo_root / 'data_analysis' / 'data'

    make_grid_figure(output_dir, data_dir)


if __name__ == '__main__':
    main()
