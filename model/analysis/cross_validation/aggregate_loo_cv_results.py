#!/usr/bin/env python3
"""
Aggregate LOO-CV results from individual fold pickle files and generate summary plots.

Run this script after all array jobs complete to:
  1. Load all fold_*.pkl files from the output directory
  2. Generate loo_cv_results.csv
  3. Generate summary plots (boxplots, phase analysis, generalisation ratios)

Usage (from HPAAxisModel/):
    python model/analysis/cross_validation/aggregate_loo_cv_results.py
    python model/analysis/cross_validation/aggregate_loo_cv_results.py --output_dir PATH
    python model/analysis/cross_validation/aggregate_loo_cv_results.py --no_plots
"""

import argparse
import pickle
from pathlib import Path
import numpy as np

# Import the plotting and I/O functions from loo_cv.py
import sys
sys.path.insert(0, str(Path(__file__).parent))
from loo_cv import (
    plot_summary, plot_phase_analysis, plot_generalisation_ratio,
    save_csv, print_summary_table,
)


def load_fold_results(output_dir):
    """Load all fold_*.pkl files from the output directory."""
    output_dir  = Path(output_dir)
    fold_files  = sorted(output_dir.glob('fold_*.pkl'))

    if not fold_files:
        raise FileNotFoundError(f"No fold_*.pkl files found in {output_dir}")

    results = []
    for fold_file in fold_files:
        with open(fold_file, 'rb') as f:
            fold = pickle.load(f)
            results.append(fold)

    results.sort(key=lambda r: r['participant'])
    print(f"Loaded {len(results)} fold results from {output_dir}")
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate LOO-CV results from pickle files and generate summary plots'
    )
    parser.add_argument(
        '--output_dir', type=str, default=None,
        help='Directory containing fold_*.pkl files (default: cross_validation/output/)'
    )
    parser.add_argument(
        '--no_plots', action='store_true',
        help='Skip plot generation — only produce the CSV results file'
    )
    args = parser.parse_args()

    output_dir = (
        Path(args.output_dir) if args.output_dir
        else Path(__file__).parent / 'output'
    )

    if not output_dir.exists():
        raise FileNotFoundError(f"Output directory not found: {output_dir}")

    print(f"\n{'='*72}")
    print("AGGREGATING LOO-CV RESULTS")
    print(f"{'='*72}")
    print(f"Output directory: {output_dir.resolve()}\n")

    results    = load_fold_results(output_dir)
    do_ts_fit  = any('ts_fitted_sse' in r for r in results)

    print("\nSaving consolidated results...")
    df = save_csv(results, output_dir, do_ts_fit=do_ts_fit)
    print_summary_table(df, do_ts_fit=do_ts_fit)

    if not args.no_plots:
        print("\nGenerating summary plots...")
        plot_summary(results, output_dir, do_ts_fit=do_ts_fit)
        plot_phase_analysis(results, output_dir, do_ts_fit=do_ts_fit)
        plot_generalisation_ratio(results, output_dir, do_ts_fit=do_ts_fit)
        print(f"\nAll plots saved to {output_dir.resolve()}")

    print("\nDone.")


if __name__ == '__main__':
    main()
