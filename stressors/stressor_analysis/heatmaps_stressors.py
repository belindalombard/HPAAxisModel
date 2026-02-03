import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate heatmaps for stressor analysis results')
parser.add_argument('--scenario', type=str, default='ACUTE',
                    help='Scenario subdirectory (e.g., ACUTE, HS1, HS2, HS3, HS4). Default: ACUTE')
parser.add_argument('--results-file', type=str, default=None,
                    help='Path to results CSV file. Default: ../output/{scenario}/results.csv')
parser.add_argument('--output-dir', type=str, default=None,
                    help='Output directory for heatmaps. Default: ../output/{scenario}/heatmaps')
args = parser.parse_args()

# Set paths based on scenario
results_file = args.results_file if args.results_file else f'../output/{args.scenario}/results.csv'
out_dir = args.output_dir if args.output_dir else f'../output/{args.scenario}/heatmaps'

results = pd.read_csv(results_file)
# Ensure output directory for heatmaps exists
os.makedirs(out_dir, exist_ok=True)
for i, row in enumerate(results['circadian_phase_describe']): 
    if row == 'rising': 
        results["circadian_phase_describe"][i] = 'rise'
    if row == 'falling': 
        results["circadian_phase_describe"][i] = 'fall'

plt.rcParams.update({'font.size': 16})
def plot_facet_grid(df, facet_col, fixed_metric, fixed_value, x_axis, y_axis, heatmap_value, aggfunc=np.mean,
                    show_numbers=False, color_range=None):
    filtered_df = df[df[fixed_metric].astype(int) == int(fixed_value)]

    if facet_col == "circadian_phase_describe":
        phase_order = ["rise", "peak", "fall", "nadir"]
        filtered_df[facet_col] = pd.Categorical(filtered_df[facet_col], categories=phase_order, ordered=True)

    if color_range is None:
        global_min = filtered_df[heatmap_value].min()
        global_max = filtered_df[heatmap_value].max()
        color_range = (global_min, global_max)

    g = sns.FacetGrid(
        filtered_df, 
        col=facet_col, 
        margin_titles=True, 
        col_wrap=2, 
        sharex=False, 
        sharey=False
    )

    def heatmap(data, **kwargs):
        pivot_table = data.pivot_table(index=y_axis, columns=x_axis, values=heatmap_value, aggfunc=aggfunc)
        pivot_table = pivot_table.sort_index(ascending=False)

        if x_axis == "ultradian_phase_radians":
            def format_ultradian_label(val: float) -> str:
                step = np.pi / 2
                # Only map the simple set: 0, π/2, π, 3π/2
                if np.isfinite(val) and np.isclose(val % step, 0.0, atol=1e-9):
                    k = int(round(val / step))
                    mapping = {
                        0: "0",
                        1: "π/2",
                        2: "π",
                        3: "3π/2",
                    }
                    if k in mapping:
                        return mapping[k]
                # For anything else, keep numeric but concise
                return f"~{val:.2f}"

            pivot_table.columns = [format_ultradian_label(col) for col in pivot_table.columns]

        sns.heatmap(
            pivot_table, 
            cmap="coolwarm", 
            annot=show_numbers, 
            fmt=".2f" if show_numbers else "",
            linewidths=0.5, 
            vmin=color_range[0], 
            vmax=color_range[1], 
            cbar=False,  
            **kwargs
        )

    g.map_dataframe(heatmap)

    norm = plt.Normalize(color_range[0], color_range[1])
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    sm.set_array([])

    # Add axes for the color bar farther to the right
    cbar_ax = g.fig.add_axes([1, 0.3, 0.02, 0.4])  # [left, bottom, width, height]
    cbar = g.fig.colorbar(sm, cax=cbar_ax)
    if heatmap_value == 'peak_diff': 
        cbar.set_label('Peak diff (nM)', rotation=270, labelpad=15)
    else: 
        # Match style with peak diff and include units for AUC
        cbar.set_label('Difference in AUC (nM·min)', rotation=270, labelpad=15)

        #cbar.set_label(heatmap_value.replace("_", " ").capitalize(), rotation=270, labelpad=15)

    # Set labels and titles
    if y_axis == 'duration':
        g.set_axis_labels("Ultradian phase", 'Duration (min)')
    else: 
        g.set_axis_labels("Ultradian phase", y_axis.capitalize())

    g.set_titles(col_template="Circadian {col_name}")
    g.fig.subplots_adjust(top=0.9)  # Adjust space for the title
    #g.fig.suptitle(f'Fixed {fixed_metric} at {fixed_value}: {heatmap_value.replace("_", " ").capitalize()}')

    plt.tight_layout()
    return g



for value in results['duration'].unique(): 
    print(value)
    fig = plot_facet_grid(results, facet_col="circadian_phase_describe", fixed_metric="duration", fixed_value=int(value), 
                 x_axis="ultradian_phase_radians", y_axis="magnitude", heatmap_value="difference_in_auc")
    try:
        fig.savefig(f'{out_dir}/fixed_duration_{value}_auc.pdf')
    except Exception:
        pass
    try:
        fig.savefig(f'{out_dir}/fixed_duration_{value}_auc.png')
    except Exception:
        pass
    plt.close(fig.fig)
    
    fig = plot_facet_grid(results, facet_col="circadian_phase_describe", fixed_metric="duration", fixed_value=int(value), 
                 x_axis="ultradian_phase_radians", y_axis="magnitude", heatmap_value="peak_diff")
    try:
        fig.savefig(f'{out_dir}/fixed_duration_{value}_peakdiff.pdf')
    except Exception:
        pass
    try:
        fig.savefig(f'{out_dir}/fixed_duration_{value}_peakdiff.png')
    except Exception:
        pass
    plt.close(fig.fig)


for value in results['magnitude'].unique(): 
    print(value)
    fig = plot_facet_grid(results, facet_col="circadian_phase_describe", fixed_metric="magnitude", fixed_value=int(value), 
                 x_axis="ultradian_phase_radians", y_axis="duration", heatmap_value="difference_in_auc")
    try:
        fig.savefig(f'{out_dir}/fixed_magnitude_{value}_auc.pdf')
    except Exception:
        pass
    try:
        fig.savefig(f'{out_dir}/fixed_magnitude_{value}_auc.png')
    except Exception:
        pass
    plt.close(fig.fig)
    
    fig = plot_facet_grid(results, facet_col="circadian_phase_describe", fixed_metric="magnitude", fixed_value=int(value), 
                 x_axis="ultradian_phase_radians", y_axis="duration", heatmap_value="peak_diff")
    plt.tight_layout()
    try:
        fig.savefig(f'{out_dir}/fixed_magnitude_{value}_peakdiff.pdf')
    except Exception:
        pass
    try:
        fig.savefig(f'{out_dir}/fixed_magnitude_{value}_peakdiff.png')
    except Exception:
        pass
    plt.close(fig.fig)
