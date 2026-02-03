"""
HPA Axis Model Runner
====================
Simple interface to run cortisol-ACTH model simulations with plotting.
This is a wrapper around run_model_multi_plot for efficiency.

Usage:
    python run_model.py --config config/base/parameters.json
    python run_model.py --plots 1,2  # Generate multiple plot types
    python run_model.py --config config/base/parameters.json --output model/output/simulation_results
"""

import argparse
import os
import sys
from pathlib import Path

# Import the efficient multi-plot runner
from run_model_multi_plot import run_simulation_and_plot

# Get the project root
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))

# Default config file
params_config = os.path.join(PROJECT_ROOT, 'config', 'base', 'parameters.json')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run HPA Axis DDE model simulation')
    parser.add_argument('--config', default=params_config, help='Path to parameter config file')
    parser.add_argument('--plots', type=str, default='1', 
                        help='Comma-separated plot types to generate (1=Separate, 3=Combined, 4=Combined+Bounds, 5=CRH). Example: "1,3,4"')
    parser.add_argument('--no_show', action='store_true',
                        help='Skip plt.show() - plots will only be saved if --output is specified')
    parser.add_argument('--output', type=str, default=None,
                        help='Output directory to save plots. If not provided, plots will be shown but not saved.')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    
    # Parse plot options
    plot_options = [int(x.strip()) for x in args.plots.split(',')]
    print(f'Plot options: {plot_options}')
    # Determine config file path
    config_file = args.config
    
    # Handle relative paths
    if not os.path.isabs(config_file):
        config_file = os.path.join(PROJECT_ROOT, config_file)
    
    # Validate config file exists
    if not os.path.exists(config_file):
        print(f"ERROR: Config file not found: {config_file}")
        sys.exit(1)
    
    print(f"{'='*60}")
    print(f"HPA Axis Model Runner")
    print(f"{'='*60}")
    print(f"Config: {config_file}")
    print(f"Plot types: {plot_options}")
    if args.output:
        print(f"Output: {args.output}")
    print(f"{'='*60}\n")
    
    # Run the simulation with efficient multi-plot runner
    try:
        result = run_simulation_and_plot(
            config_file=config_file,
            plot_options=plot_options,
            output_dir=args.output
        )
        
        if result['success']:
            print(f"\n✓ SUCCESS!")
            print(f"{'='*60}")
            print(f"Generated {result['num_plots']} plot(s)")
            
            if result['plot_files']:
                print(f"\nPlot files:")
                for pf in result['plot_files']:
                    print(f"  • {pf}")
            
            # Show plots if not in no_show mode and no output specified
            if not args.no_show and not args.output:
                print(f"\nTIP: Use --output <directory> to save plots, or --no_show to suppress display")
                import matplotlib.pyplot as plt
                print(f"WARNING: Plots generated but not displayed (use --output to save them)")
            
            return 0
        else:
            print(f"\n✗ Simulation failed!")
            return 1
            
    except Exception as e:
        print(f"\n✗ ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
