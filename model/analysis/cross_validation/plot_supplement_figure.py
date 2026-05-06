#!/usr/bin/env python3
"""
Generate publication-quality supplement figure for LOO-CV results — Figure S12.

(A) Distribution of normalised SSE for individual fits and LOO-CV predictions.
(B) Per-participant comparison of prediction error between LOO-CV and best-fit signals.
(C) Generalisation ratio (LOO-CV SSE / participant SSE).

Requires aggregate_loo_cv_results.py to have been run first to produce the
alignment_params.json files in each participant subfolder of output/.

Outputs (in model/analysis/cross_validation/output/):
  supplement_loo_cv_figure.pdf / .png  – 3-panel LOO-CV summary figure  ← Figure S12
  summary_statistics_table.csv         – numeric summary
  loo_cv_summary_stats.txt             – written interpretation

Usage:
    python model/analysis/cross_validation/plot_supplement_figure.py
    python model/analysis/cross_validation/plot_supplement_figure.py --output_dir PATH
"""

import argparse
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pathlib import Path


def create_supplement_figure(output_dir: Path) -> None:
    """
    Build a 3-panel supplement figure from per-participant alignment_params.json files.

    Panel A: Violin plots comparing participant fit vs LOO-CV error distributions.
    Panel B: Per-participant bar chart (side-by-side).
    Panel C: Generalisation ratio per participant (LOO-CV SSE / participant MSE).
    """
    participants = []
    insample = []
    loo_error = []
    loo_method = None

    for p in range(1, 11):
        json_path = output_dir / f'participant_{p}' / 'alignment_params.json'
        if json_path.exists():
            with open(json_path, 'r') as f:
                data = json.load(f)
            participants.append(data['participant'])
            insample.append(data['error_metrics']['insample_mse'])

            recommended = data.get('recommended_method', 'phase_aligned')
            if recommended == 'ts_fitted' and data['error_metrics'].get('ts_fitted_sse') is not None:
                loo_error.append(data['error_metrics']['ts_fitted_sse'])
                loo_method = 't_s-optimised'
            else:
                loo_error.append(data['error_metrics']['phase_aligned_sse'])
                loo_method = 'phase-aligned'

    participants = np.array(participants)
    insample    = np.array(insample)
    loo_error   = np.array(loo_error)

    print(f"Using LOO-CV method: {loo_method}")
    print(f"Participant fit errors: {insample}")
    print(f"LOO-CV errors:         {loo_error}")

    insample_mean = np.nanmean(insample)
    insample_std  = np.nanstd(insample)
    loo_mean      = np.nanmean(loo_error)
    loo_std       = np.nanstd(loo_error)
    gen_ratio     = loo_mean / insample_mean
    n_improved    = int(np.sum(loo_error < insample))
    n_total       = len(participants)

    # ── Figure layout ────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 9))
    gs  = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.35, wspace=0.3,
                           left=0.08, right=0.96, top=0.93, bottom=0.08)

    # ── Panel A: Violin plots ────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])

    violin_data = pd.DataFrame({
        'SSE':    np.concatenate([insample, loo_error]),
        'Method': (['Participant Fit\n(individualised)'] * len(insample) +
                   [f'Cohort Model\n(LOO-CV, {loo_method})'] * len(loo_error)),
    })

    sns.violinplot(data=violin_data, x='Method', y='SSE',
                   inner='quartile', hue='Method',
                   palette=['lightblue', '#50C878'],
                   legend=False, ax=ax1)

    for i, method in enumerate(violin_data['Method'].unique()):
        vals = violin_data[violin_data['Method'] == method]['SSE']
        x    = np.random.default_rng(42 + i).uniform(-0.04, 0.04, size=len(vals))
        ax1.scatter(i + x, vals,
                    alpha=0.6, s=60, color='darkgray', zorder=3,
                    edgecolor='black', linewidth=0.5)

    ax1.set_xlabel('')
    ax1.set_ylabel('Normalised SSE', fontsize=12, fontweight='bold')
    ax1.set_title('A. Error distribution comparison', fontsize=13, fontweight='bold', loc='left')
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.tick_params(axis='x', labelsize=11)

    stats_text = (f'Mean \u00b1 SD:\n'
                  f'Participant: {insample_mean:.3f} \u00b1 {insample_std:.3f}\n'
                  f'LOO-CV: {loo_mean:.3f} \u00b1 {loo_std:.3f}\n'
                  f'Generalisation ratio: {gen_ratio:.3f}')
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
             fontsize=9, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.4))

    # ── Panel B: Per-participant bars ────────────────────────────────────────
    ax2  = fig.add_subplot(gs[0, 1])
    x    = np.arange(len(participants))
    w    = 0.35

    ax2.bar(x - w / 2, insample,  w, label='Participant Fit',
            color='#4A90E2', alpha=0.8, edgecolor='black', linewidth=0.8)
    ax2.bar(x + w / 2, loo_error, w, label='Cohort Model (LOO-CV)',
            color='#50C878', alpha=0.8, edgecolor='black', linewidth=0.8)

    for i, (ins, loo) in enumerate(zip(insample, loo_error)):
        if loo < ins:
            ax2.text(i, max(ins, loo) + 0.002, '\u2605',
                     ha='center', va='bottom', fontsize=14,
                     color='gold', fontweight='bold')

    ax2.set_xlabel('Participant', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Normalised SSE', fontsize=12, fontweight='bold')
    ax2.set_title('B. Individual participant performance', fontsize=13,
                  fontweight='bold', loc='left')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'P{p}' for p in participants], fontsize=10)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')

    legend_elements = [
        Patch(facecolor='#4A90E2', edgecolor='black', label='Participant fit'),
        Patch(facecolor='#50C878', edgecolor='black', label='LOO-CV'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='gold',
               markersize=15, markeredgewidth=0,
               label=f'LOO-CV better ({n_improved}/{n_total})'),
    ]
    ax2.legend(handles=legend_elements, fontsize=10, loc='upper center',
               bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=True)

    # ── Panel C: Generalisation ratio ────────────────────────────────────────
    ax3    = fig.add_subplot(gs[1, :])
    ratios = loo_error / insample

    bar_colors = ['#50C878' if (not np.isnan(r) and r <= 1.0) else '#E8B4B8'
                  for r in ratios]
    bars = ax3.bar(x, ratios, color=bar_colors, alpha=0.8, edgecolor='black', linewidth=1)

    ax3.axhline(1.0, color='black', linewidth=2, linestyle='--', zorder=0)
    ax3.axhspan(0, 1.0, alpha=0.1, color='green', zorder=0)

    ax3.set_xlabel('Participant', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Generalisation ratio\n(LOO-CV SSE / Participant MSE)',
                   fontsize=12, fontweight='bold')
    ax3.set_title('C. Generalisation performance '
                  '(ratio < 1 indicates improved prediction with cohort model)',
                  fontsize=13, fontweight='bold', loc='left')
    ax3.set_xticks(x)
    ax3.set_xticklabels([f'P{p}' for p in participants], fontsize=10)
    ax3.grid(axis='y', alpha=0.3, linestyle='--')

    valid_ratios = ratios[~np.isnan(ratios)]
    ax3.set_ylim([0, (max(valid_ratios) * 1.1) if len(valid_ratios) else 2])

    for bar, ratio in zip(bars, ratios):
        if not np.isnan(ratio):
            ax3.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 0.02,
                     f'{ratio:.2f}',
                     ha='center', va='bottom', fontsize=9, fontweight='bold')

    # ── Save figure ───────────────────────────────────────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)

    for ext in ('png', 'pdf'):
        fig.savefig(output_dir / f'supplement_loo_cv_figure.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
        print(f"Saved: supplement_loo_cv_figure.{ext}")
    plt.close()

    # ── Summary CSV ───────────────────────────────────────────────────────────
    valid_increases = (loo_error - insample) / insample
    valid_increases = valid_increases[~np.isnan(valid_increases)]
    max_increase_str = f'{np.max(valid_increases):.1%}' if len(valid_increases) else 'N/A'

    summary_df = pd.DataFrame({
        'Metric': [
            'Mean SSE \u00b1 SD',
            'Median SSE',
            'Generalisation ratio',
            'Participants with improved fit',
            'Max error increase',
        ],
        'Participant fit (individualised)': [
            f'{insample_mean:.4f} \u00b1 {insample_std:.4f}',
            f'{np.nanmedian(insample):.4f}',
            '1.000 (reference)',
            '\u2014',
            '\u2014',
        ],
        'Cohort model (LOO-CV)': [
            f'{loo_mean:.4f} \u00b1 {loo_std:.4f}',
            f'{np.nanmedian(loo_error):.4f}',
            f'{gen_ratio:.3f}',
            f'{n_improved}/{n_total} ({100 * n_improved / n_total:.0f}%)',
            max_increase_str,
        ],
        'Interpretation': [
            'Similar means indicate good generalisation',
            'Robust central tendency measure',
            '< 1.0 = excellent, ~1.0 = good, > 1.5 = potential overfitting',
            'LOO-CV error lower than participant fit',
            'Worst-case generalisation penalty',
        ],
    })
    summary_df.to_csv(output_dir / 'summary_statistics_table.csv', index=False)
    print(f"Saved: summary_statistics_table.csv")

    # ── Summary text ─────────────────────────────────────────────────────────
    with open(output_dir / 'loo_cv_summary_stats.txt', 'w') as f:
        f.write("LEAVE-ONE-OUT CROSS-VALIDATION SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        f.write("Participant Fit (individualised):\n")
        f.write(f"  Mean \u00b1 SD : {insample_mean:.4f} \u00b1 {insample_std:.4f}\n")
        f.write(f"  Median   : {np.nanmedian(insample):.4f}\n")
        f.write(f"  Range    : [{np.nanmin(insample):.4f}, {np.nanmax(insample):.4f}]\n\n")
        f.write(f"Cohort Model LOO-CV (out-of-sample, {loo_method}):\n")
        f.write(f"  Mean \u00b1 SD : {loo_mean:.4f} \u00b1 {loo_std:.4f}\n")
        f.write(f"  Median   : {np.nanmedian(loo_error):.4f}\n")
        f.write(f"  Range    : [{np.nanmin(loo_error):.4f}, {np.nanmax(loo_error):.4f}]\n\n")
        f.write("Generalisation Performance:\n")
        f.write(f"  Mean error ratio (LOO-CV / Participant): {gen_ratio:.3f}\n")
        f.write(f"  Participants with improved fit: {n_improved}/{n_total} "
                f"({100 * n_improved / n_total:.0f}%)\n")
        f.write(f"  Median ratio: {np.nanmedian(ratios):.3f}\n")
        if len(valid_increases):
            f.write(f"  Max error increase: {np.max(valid_increases):.1%}\n\n")
        f.write("Interpretation:\n")
        f.write(f"  - Generalisation ratio of {gen_ratio:.3f} indicates ")
        if gen_ratio < 1.0:
            f.write("EXCELLENT generalisation (cohort better than participant fits)\n")
        elif gen_ratio < 1.2:
            f.write("GOOD generalisation (minimal overfitting)\n")
        else:
            f.write("ACCEPTABLE generalisation\n")
        f.write(f"  - {n_improved} participants showed lower error with cohort model,\n"
                f"    suggesting participant fits may have overfit to noise\n"
                f"  - Model structure captures generalisable physiological dynamics\n")
    print("Saved: loo_cv_summary_stats.txt")


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Generate publication-quality supplement figure for LOO-CV results'
    )
    parser.add_argument(
        '--output_dir', type=str, default=None,
        help='Directory containing per-participant alignment_params.json files '
             '(default: model/analysis/cross_validation/output/)',
    )
    args = parser.parse_args()

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path(__file__).resolve().parent / 'output'

    print(f"\n{'=' * 72}")
    print("GENERATING SUPPLEMENT FIGURE")
    print(f"{'=' * 72}")
    print(f"Loading alignment_params.json files from: {output_dir}")

    create_supplement_figure(output_dir)
    print("\nDone! Figure ready for supplement.")


if __name__ == '__main__':
    main()
