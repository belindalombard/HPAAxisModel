import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import seaborn as sns
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Generate heatmaps for heart surgery stressor analysis')
parser.add_argument('--scenario', type=str, default='HS',
                    help='Scenario subdirectory (e.g., HS, HS1, HS2, HS3, HS4). Default: HS')
parser.add_argument('--results-file', type=str, default=None,
                    help='Path to results CSV file. Default: ../output/{scenario}/results_heart_surgery.csv')
args = parser.parse_args()

# Set paths based on scenario
results_file = args.results_file if args.results_file else f'../output/{args.scenario}/results_heart_surgery.csv'

results = pd.read_csv(results_file)

def plot_facet_grid(df, facet_col, x_axis, y_axis, heatmap_value, aggfunc=np.mean,
                    show_numbers=False, color_range=None):
    if color_range is None:
        global_min = df[heatmap_value].min()
        global_max = df[heatmap_value].max()
        color_range = (global_min, global_max)

    # Create the FacetGrid
    g = sns.FacetGrid(
        df, 
        col=facet_col, 
        margin_titles=True, 
        col_wrap=2, 
        sharex=False, 
        sharey=False
    )

    def heatmap(data, **kwargs):
        pivot_table = data.pivot_table(index=y_axis, columns=x_axis, values=heatmap_value, aggfunc=aggfunc)
        pivot_table = pivot_table.sort_index(ascending=False)

        sns.heatmap(
            pivot_table, 
            cmap="coolwarm", 
            annot=show_numbers, 
            fmt=".2f" if show_numbers else "",
            linewidths=0.5, 
            vmin=color_range[0], 
            vmax=color_range[1], 
            cbar=False,  # Suppress individual color bars
            **kwargs
        )

    g.map_dataframe(heatmap)

    # Create and position a single color bar
    norm = plt.Normalize(color_range[0], color_range[1])
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    sm.set_array([])

    # Add axes for the color bar farther to the right
# Add a single color bar across all heatmaps
    cbar = g.fig.colorbar(sm, ax=g.axes.ravel().tolist(), orientation="vertical", fraction=0.02, pad=0.02)
    cbar.set_label(heatmap_value.replace("_", " ").capitalize(), rotation=270, labelpad=15)

    #cbar.set_label(heatmap_value.replace("_", " ").capitalize(), rotation=270, labelpad=15)

    # Set labels and titles
    g.set_axis_labels(x_axis.replace("_", " ").capitalize(), y_axis.replace("_", " ").capitalize())
    g.set_titles(col_template="Duration: {col_name} minutes")
    g.fig.subplots_adjust(top=0.9)  # Adjust space for the title

    plt.tight_layout()
    return g


# Generate the heatmap plots
fig = plot_facet_grid(
    results, 
    facet_col="duration", 
    x_axis="time_in_scope", 
    y_axis="magnitude", 
    heatmap_value="num_peaks" 
)
fig.savefig('../output/heatmaps/heatmap_heart_surgery.png')

crh_values = np.load("../output/HS1/250.00_120.00_1350.00/crh.npy")
cort_values = np.load("../output/HS1/250.00_120.00_1350.00/simulated_values_original.npy")



cort_values = cort_values[:, 1]

start_idx, end_idx = 0, 2880
time_range = np.arange(start_idx, end_idx)

unique_time_in_scope = results['time_in_scope'].unique()

time_in_scope_filtered = [t for t in unique_time_in_scope if start_idx <= t < end_idx]

fig, ax1 = plt.subplots(figsize=(10, 5))

ax1.plot(time_range, crh_values[start_idx:end_idx], label="CRH Levels", color="blue")
ax1.scatter(time_in_scope_filtered, crh_values[time_in_scope_filtered], color='red', label="Time in scope (CRH)", zorder=3)

ax1.set_xlabel("Time (minutes)")
ax1.set_ylabel("CRH Levels", color="blue")
ax1.tick_params(axis='y', labelcolor="blue")
ax1.grid(True)

ax2 = ax1.twinx()
ax2.plot(time_range, cort_values[start_idx:end_idx], label="CORT Levels", color="green", linestyle="dashed")
#ax2.scatter(time_in_scope_filtered, cort_values[time_in_scope_filtered], color='orange', label="Time in Scope (CORT)", zorder=3)

ax2.set_ylabel("CORT Levels", color="green")
ax2.tick_params(axis='y', labelcolor="green")

for t in time_in_scope_filtered:
    ax1.axvline(x=t, color="red", linestyle="dotted", alpha=0.8, linewidth=1.2)

max = end_idx
start_time='09:00:00'
start_time = datetime.strptime(start_time, "%H:%M:%S")
ax1.set_xticks([i for i in range(0, max + 1, 360)])
ax1.set_xticklabels([(start_time + timedelta(minutes=i)).strftime("%H:%M") for i in range(0, max+1, 360)])

#plt.title("CRH and CORT Levels with Time in Scope Markers")
plt.savefig('../output/heart_surgery_timepoints.png')