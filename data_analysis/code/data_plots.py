import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.stats import zscore
from scipy.signal import correlate
import numpy as np
import matplotlib.dates as mdates
import seaborn as sns

# Where is the data? 
data_dir = 'model/data'
acth_data_file = f'data_shifted_ACTH_1min_09_00.csv'
cort_data_file = f'data_shifted_Cortisol_1min_09_00.csv'

plt.rcParams.update({"font.size": 24})

# Get the raw data
acth_data = pd.read_csv(f'{data_dir}/{acth_data_file}')
cort_data = pd.read_csv(f'{data_dir}/{cort_data_file}')
acth_data["x"] = pd.to_datetime(acth_data["x"], errors="coerce")
cort_data["x"] = pd.to_datetime(cort_data["x"], errors="coerce")
acth_data["test"] = acth_data["test"].astype(int)
cort_data["test"] = cort_data["test"].astype(int)
cort_data = cort_data.sort_values(by=["test", "x"])
acth_data = acth_data.sort_values(by=["test", "x"])

# Merge data on time and participant
merged_data = pd.merge(acth_data, cort_data, on=["x", "test"], suffixes=("_acth", "_cort"))

# ---- Combined Plot ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout = True)

# ACTH plot
for participant, group in merged_data.groupby("test"):
    ax1.plot(group["x"], group["y_acth"], label=f"P{participant}", alpha=0.7)
ax1.set_ylabel("ACTH (pmol/L)")
ax1.grid(True, linestyle="--", alpha=0.5)

# CORT plot
for participant, group in merged_data.groupby("test"):
    ax2.plot(group["x"], group["y_cort"], label=f"P{participant}", alpha=0.7)
ax2.set_ylabel("CORT (nmol/L)")
ax2.set_xlabel("Clock time (hrs)")
ax2.grid(True, linestyle="--", alpha=0.5)

# Format shared x-axis as hours

fig.align_ylabels([ax1, ax2])
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax2.xaxis.set_major_locator(mdates.HourLocator(interval=4))


# Optional: Only one legend for clarity
# fig.legend(loc="upper right", fontsize=12)

plt.tight_layout()
plt.subplots_adjust(hspace=0.05)  # Reduce space between panels

# Save combined figure
output_dir = "data_analysis/output/data_plots"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
plt.savefig(f"{output_dir}/data_acth_cort_combined.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{output_dir}/data_acth_cort_combined.png", dpi=300, bbox_inches="tight")

plt.show()


plt.close()  # Close the figure to start fresh for the next plot
print("Created 2")

# 3. Plot Time-lag Correlation
plt.figure(figsize=(12, 4))
ax3 = plt.gca()
correlations = []
time_lags = []

for participant, group in merged_data.groupby("test"):
    acth_values = group["y_acth"].values
    cort_values = group["y_cort"].values
    lag_min, lag_max = -250, 250
    correlation = correlate(acth_values - np.mean(acth_values), cort_values - np.mean(cort_values), mode="full")
    correlation /= np.max(correlation)  # Normalize
    lags = np.arange(-len(acth_values) + 1, len(acth_values))
    valid_indices = (lags >= lag_min) & (lags <= lag_max)
    lags_filtered = lags[valid_indices]
    correlation_filtered = correlation[valid_indices]
    correlations.append(correlation)
    peak_lag = lags_filtered[np.argmax(correlation_filtered)]
    time_lags.append(peak_lag)


    ax3.plot(lags_filtered, correlation_filtered, label=f"Participant {participant}", alpha=0.7)

mean_lag = np.mean(time_lags)
std_lag = np.std(time_lags)

ax3.axvline(mean_lag, color="red", linestyle="--", linewidth=1.5, label=f"Mean lag = {mean_lag:.2f} ± {std_lag:.2f}")

#ax3.text(mean_lag + 5, 0.9, f"{mean_lag:.1f} ± {std_lag:.1f}", color="red", fontsize=10)
ax3.set_xlabel("Lag (min)")
ax3.set_ylabel("Correlation")
ax3.grid(True, linestyle="--", alpha=0.5)
plt.savefig(f"{output_dir}/timelag.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{output_dir}/timelag.png", dpi=300, bbox_inches="tight")
plt.show()
plt.close()

time_lags = np.array(time_lags)
mean_lag = np.mean(time_lags)
std_lag = np.std(time_lags)
min_lag = np.min(time_lags)
max_lag = np.max(time_lags)


