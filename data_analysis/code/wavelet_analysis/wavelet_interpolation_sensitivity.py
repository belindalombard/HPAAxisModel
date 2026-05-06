"""
Wavelet analysis sensitivity to interpolation (resampling) frequency.

Reviewer context: the data were originally sampled at 10 or 20-min intervals and
interpolated to 1 min before analysis.  This script tests whether the detected
dominant ultradian period depends on the interpolation frequency by subsampling
the 1-min data to coarser resolutions and repeating the wavelet analysis.

Resampling frequencies tested: 1, 5, 10 min.
Note: interpolation below the original Nyquist limit does not add information.
At dt=10 the Nyquist period is 20 min; the ultradian band (50–300 min) is safely
above this limit for all tested dt values.

Outputs (in data_analysis/code/wavelet_analysis/output/interpolation_sensitivity/):
  interp_spectra_Cortisol.pdf         – group-mean spectra for each dt
  interp_boxplot_Cortisol.pdf         – dominant period distribution per dt
  interp_individual_Cortisol.pdf      – per-participant dominant periods across dt
  interp_summary_Cortisol.csv         – participant-level dominant periods (table)
  (same for ACTH)

Usage:
    python data_analysis/code/wavelet_analysis/wavelet_interpolation_sensitivity.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyboat import WAnalyzer

# ── configuration ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_BASE_DIR   = os.path.dirname(os.path.dirname(_SCRIPT_DIR))  # data_analysis/

DATA_FILE    = os.path.join(_BASE_DIR, "data", "data_per_participant_pyboat.csv")
OUTPUT_DIR   = os.path.join(_SCRIPT_DIR, "output", "interpolation_sensitivity")
N_PARTICIPANTS = 10

BASE_DT      = 1       # resolution of the stored data (min)
CUTOFF       = 1440    # sinc-filter cutoff (min)
WINDOW_SIZE  = 720     # amplitude normalisation window (min)

MIN_PERIOD   = 50
MAX_PERIOD   = 300
N_PERIODS    = 500

ULTRADIAN_LO = 50
ULTRADIAN_HI = 300

# Resampling frequencies to test (in minutes)
DT_VALUES  = [1, 5, 10]
DT_LABELS  = ["1 min\n(as-is)", "5 min", "10 min\n(original)"]

HORMONES = ["Cortisol", "ACTH"]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── helpers ────────────────────────────────────────────────────────────────────
def interpolate_nans(signal):
    s = np.array(signal, dtype=float)
    mask = np.isnan(s)
    if mask.any():
        idx = np.arange(len(s))
        s[mask] = np.interp(idx[mask], idx[~mask], s[~mask])
    return s


def subsample(signal, dt, base_dt=BASE_DT):
    """Subsample a 1-min signal to the target resolution by taking every nth sample."""
    step = max(1, int(round(dt / base_dt)))
    return signal[::step]


def run_wavelet(signal, dt, cutoff=CUTOFF, window_size=WINDOW_SIZE):
    """Detrend, normalise, and compute wavelet power spectrum at a given dt."""
    periods = np.linspace(MIN_PERIOD, MAX_PERIOD, N_PERIODS)
    wAn = WAnalyzer(periods, dt, time_unit_label="min")
    trend      = wAn.sinc_smooth(signal, T_c=cutoff)
    detrended  = signal - trend
    norm_sig   = wAn.normalize_amplitude(detrended, window_size=window_size)
    modulus, _ = wAn.compute_spectrum(norm_sig, do_plot=False)
    return modulus, periods


def dominant_ultradian_period(modulus, periods,
                               lo=ULTRADIAN_LO, hi=ULTRADIAN_HI):
    mask  = (periods >= lo) & (periods <= hi)
    power = modulus[mask].mean(axis=1)
    return periods[mask][np.argmax(power)]


def mean_power_spectrum(modulus):
    return modulus.mean(axis=1)


# ── main analysis ──────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_FILE)

cmap   = plt.cm.plasma
colors = [cmap(v) for v in np.linspace(0.15, 0.75, len(DT_VALUES))]

for hormone in HORMONES:
    print(f"\n{'='*60}\n  {hormone}\n{'='*60}")

    results = {}

    for dt in DT_VALUES:
        dom_periods  = []
        mean_spectra = []
        for p in range(1, N_PARTICIPANTS + 1):
            raw     = interpolate_nans(df[f"{hormone}_{p}"].values)
            signal  = subsample(raw, dt)
            modulus, periods = run_wavelet(signal, dt)
            dom_periods.append(dominant_ultradian_period(modulus, periods))
            mean_spectra.append(mean_power_spectrum(modulus))
        results[dt] = dict(dominant_periods=dom_periods, mean_spectra=mean_spectra,
                           periods=periods)
        print(f"  dt={dt:2d} min  n_samples={len(subsample(interpolate_nans(df[f'{hormone}_1'].values), dt)):5d}  "
              f"median period={np.median(dom_periods):.1f} min")

    ref_periods = results[DT_VALUES[0]]["periods"]

    # ── Figure 1: Group-mean spectra ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, dt in enumerate(DT_VALUES):
        stack    = np.vstack(results[dt]["mean_spectra"])
        for_plot = np.array([np.interp(ref_periods, results[dt]["periods"], row)
                             for row in stack])
        mean_s   = for_plot.mean(axis=0)
        sem_s    = for_plot.std(axis=0) / np.sqrt(N_PARTICIPANTS)
        ax.plot(ref_periods, mean_s, color=colors[i], lw=1.8,
                label=f"dt = {dt} min")
        ax.fill_between(ref_periods, mean_s - sem_s, mean_s + sem_s,
                        color=colors[i], alpha=0.15)
    ax.axvspan(ULTRADIAN_LO, ULTRADIAN_HI, alpha=0.07, color="gray",
               label=f"Ultradian range ({ULTRADIAN_LO}–{ULTRADIAN_HI} min)")
    ax.set_xlabel("Period (min)", fontsize=13)
    ax.set_ylabel("Mean wavelet power (a.u.)", fontsize=13)
    ax.set_title(
        f"{hormone}: sensitivity of mean power spectrum\n"
        f"to interpolation frequency  (n = {N_PARTICIPANTS})",
        fontsize=13,
    )
    ax.legend(fontsize=10, framealpha=0.9)
    ax.set_xlim(MIN_PERIOD, MAX_PERIOD)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/interp_spectra_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/interp_spectra_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved group spectra figure.")

    # ── Figure 2: Boxplot of dominant periods ─────────────────────────────────
    all_data = [results[dt]["dominant_periods"] for dt in DT_VALUES]
    fig, ax_bp = plt.subplots(figsize=(7, 5))
    bp = ax_bp.boxplot(all_data, patch_artist=True,
                       medianprops=dict(color="black", lw=2))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
    rng = np.random.default_rng(42)
    for i, vals in enumerate(all_data):
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax_bp.scatter(np.full(len(vals), i + 1) + jitter, vals,
                      color="black", s=18, zorder=3, alpha=0.6)
    ax_bp.set_xticks(range(1, len(DT_VALUES) + 1))
    ax_bp.set_xticklabels(DT_LABELS, fontsize=11)
    ax_bp.set_xlabel("Interpolation frequency (dt)", fontsize=13)
    ax_bp.set_ylabel("Dominant ultradian period (min)", fontsize=13)
    ax_bp.set_title(
        f"{hormone}: dominant period vs. interpolation frequency\n"
        f"(n = {N_PARTICIPANTS}, band = {ULTRADIAN_LO}–{ULTRADIAN_HI} min)",
        fontsize=12,
    )
    ax_bp.set_ylim(ULTRADIAN_LO - 10, ULTRADIAN_HI + 10)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/interp_boxplot_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/interp_boxplot_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved boxplot figure.")

    # ── Figure 3: Individual participant lines ────────────────────────────────
    fig, ax_ind = plt.subplots(figsize=(7, 5))
    x_pos  = np.arange(1, len(DT_VALUES) + 1)
    cmap_p = plt.cm.tab10
    for p in range(N_PARTICIPANTS):
        vals = [results[dt]["dominant_periods"][p] for dt in DT_VALUES]
        ax_ind.plot(x_pos, vals,
                    color=cmap_p(p / max(N_PARTICIPANTS, 10)),
                    marker="o", markersize=5, lw=1.2, alpha=0.75,
                    label=f"P{p + 1}")
    median_vals = [np.median(results[dt]["dominant_periods"]) for dt in DT_VALUES]
    ax_ind.plot(x_pos, median_vals, color="black", marker="D",
                markersize=7, lw=2.2, zorder=5, label="Median")
    ax_ind.set_xticks(x_pos)
    ax_ind.set_xticklabels(DT_LABELS, fontsize=11)
    ax_ind.set_xlabel("Interpolation frequency (dt)", fontsize=13)
    ax_ind.set_ylabel("Dominant ultradian period (min)", fontsize=13)
    ax_ind.set_title(
        f"{hormone}: individual dominant periods vs. interpolation frequency\n"
        f"(n = {N_PARTICIPANTS})",
        fontsize=12,
    )
    ax_ind.set_ylim(ULTRADIAN_LO - 10, ULTRADIAN_HI + 10)
    ax_ind.legend(fontsize=9, ncol=2, framealpha=0.9, loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/interp_individual_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/interp_individual_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved individual figure.")

    # ── CSV ────────────────────────────────────────────────────────────────────
    rows = []
    for p in range(N_PARTICIPANTS):
        row = {"participant": p + 1}
        for dt in DT_VALUES:
            row[f"dt_{dt}min"] = round(results[dt]["dominant_periods"][p], 1)
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.loc["median"] = ["median"] + \
        [round(np.median(results[dt]["dominant_periods"]), 1) for dt in DT_VALUES]
    summary.loc["std"] = ["std"] + \
        [round(np.std(results[dt]["dominant_periods"]), 1) for dt in DT_VALUES]
    summary.to_csv(f"{OUTPUT_DIR}/interp_summary_{hormone}.csv", index=False)
    print(f"  Saved CSV.")

print("\nDone.")
