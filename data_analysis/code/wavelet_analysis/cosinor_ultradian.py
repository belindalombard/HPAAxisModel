"""
Multi-component cosinor analysis of ultradian rhythmicity.

The cosinor method provides an independent confirmation of ultradian periodicity
that does not rely on a sliding window (addressing the reviewer concern about
circularity in window selection for the wavelet analysis).

Method
------
For each candidate period T, a single-component cosinor model is fitted by OLS
using CosinorPy (Moskon, BMC Bioinformatics, 2020; doi:10.1186/s12859-020-03830-w):

    y(t) = M  +  A·cos(2π t / T)  +  B·sin(2π t / T)

An F-test (ANOVA) evaluates whether the rhythmic terms (A, B) are jointly
non-zero (H0: no rhythm at period T).  Scanning across periods produces a
cosinor periodogram analogous to a power spectrum.

The minimum period tested is 40 min, which is the Nyquist limit for the
20-min original sampling interval (the more conservative of the two sampling
rates in the dataset).  This avoids spurious rhythms that can arise from
interpolation artefacts at shorter periods.

Pre-processing matches the wavelet analysis pipeline:
  1. Linear NaN interpolation
  2. Sinc-filter detrending at T_c = 1440 min (removes circadian trend;
     identical to pyboat.core.sinc_smooth used in wavelet_window_sensitivity.py)

Significant rhythms are identified using the raw F-test p-value (α = 0.05)
and the dominant ultradian period (highest F-statistic within 100–200 min)
is reported per participant.

Outputs (in data_analysis/code/wavelet_analysis/output/cosinor_analysis/):
  cosinor_periodogram_Cortisol.pdf      – F-statistic scans for all participants
  cosinor_dominant_period_Cortisol.pdf  – boxplot of dominant periods
  cosinor_group_mean_Cortisol.pdf       – group-median F-statistic periodogram
  cosinor_summary_Cortisol.csv          – per-participant results table
  (same for ACTH)

Usage:
    python data_analysis/code/wavelet_analysis/cosinor_ultradian.py
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from CosinorPy import cosinor1
from pyboat import core as pyboat_core

# ── configuration ──────────────────────────────────────────────────────────────
_SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
_BASE_DIR      = os.path.dirname(os.path.dirname(_SCRIPT_DIR))  # data_analysis/

DATA_FILE      = os.path.join(_BASE_DIR, "data", "data_per_participant_pyboat.csv")
OUTPUT_DIR     = os.path.join(_SCRIPT_DIR, "output", "cosinor_analysis")
N_PARTICIPANTS = 10

DT             = 1          # data resolution (1 min)
MIN_PERIOD     = 50         # minimum period (min) – Nyquist-safe
MAX_PERIOD     = 300        # maximum period (min)
PERIOD_STEP    = 1.0        # step between tested periods (min)

ULTRADIAN_LO   = 50        # ultradian band for dominant-period extraction
ULTRADIAN_HI   = 300

CUTOFF         = 1440       # sinc-filter cutoff (min) – removes circadian trend

ALPHA          = 0.05       # significance threshold (before correction)
HORMONES       = ["Cortisol", "ACTH"]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── helpers ────────────────────────────────────────────────────────────────────
def interpolate_nans(signal):
    """Linearly interpolate any NaN values."""
    s = np.array(signal, dtype=float)
    mask = np.isnan(s)
    if mask.any():
        idx = np.arange(len(s))
        s[mask] = np.interp(idx[mask], idx[~mask], s[~mask])
    return s


def sinc_detrend(signal, cutoff=CUTOFF, dt=DT):
    """Remove low-frequency trend using a sinc filter at *cutoff* period.

    Identical to pyboat WAnalyzer.sinc_smooth(signal, T_c=cutoff) used in
    wavelet_window_sensitivity.py, ensuring consistent pre-processing.
    """
    trend = pyboat_core.sinc_smooth(signal, cutoff, dt)
    return signal - trend


def fit_cosinor_period(t, y, period):
    """Fit a single-component cosinor at *period* using CosinorPy.

    Returns F-statistic, its p-value, amplitude and acrophase.
    """
    data = pd.DataFrame({"x": t, "y": y})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results, amp, acr, statistics = cosinor1.test_cosinor_single(
            data, period=period, corrected=True
        )
    return results.fvalue, statistics["F-test"], amp, acr


def cosinor_periodogram(t, signal, periods):
    """Scan over *periods* and return the F-statistic at each one."""
    f_stats = np.empty(len(periods))
    p_vals  = np.empty(len(periods))
    for i, T in enumerate(periods):
        F, p, *_ = fit_cosinor_period(t, signal, T)
        f_stats[i] = F if np.isfinite(F) else 0.0
        p_vals[i]  = p
    return f_stats, p_vals


# ── main analysis ──────────────────────────────────────────────────────────────
periods = np.arange(MIN_PERIOD, MAX_PERIOD + PERIOD_STEP, PERIOD_STEP)
df      = pd.read_csv(DATA_FILE)
t_vec   = df["x"].values.astype(float)   # minutes from reference time

for hormone in HORMONES:
    print(f"\n{'='*60}\n  {hormone}\n{'='*60}")

    all_f_stats   = []
    dominant_info = []

    for p in range(1, N_PARTICIPANTS + 1):
        signal  = interpolate_nans(df[f"{hormone}_{p}"].values)
        f_stats, p_vals = cosinor_periodogram(t_vec, signal, periods)
        all_f_stats.append(f_stats)

        ultra_mask = (periods >= ULTRADIAN_LO) & (periods <= ULTRADIAN_HI)
        f_ultra    = f_stats[ultra_mask]
        best_idx   = np.argmax(f_ultra)
        best_T     = periods[ultra_mask][best_idx]
        best_F     = f_ultra[best_idx]
        p_val      = p_vals[ultra_mask][best_idx]

        _, _, amplitude, acrophase = fit_cosinor_period(t_vec, signal, best_T)

        dominant_info.append(dict(
            participant   = p,
            dominant_T    = round(best_T, 1),
            F_stat        = round(best_F, 2),
            p_value       = round(p_val, 4),
            significant   = p_val < ALPHA,
            amplitude     = round(float(amplitude), 3),
            acrophase_rad = round(float(acrophase), 3),
            acrophase_min = round(float(acrophase) / (2 * np.pi) * best_T, 1),
        ))
        status = "✓ sig" if p_val < ALPHA else "  n.s."
        print(f"  P{p:2d}: dominant T = {best_T:.1f} min  "
              f"F = {best_F:.2f}  p = {p_val:.4f}  {status}")

    summary_df = pd.DataFrame(dominant_info)
    summary_df.to_csv(f"{OUTPUT_DIR}/cosinor_summary_{hormone}.csv", index=False)

    # ── Figure 1 : Periodograms for all participants ──────────────────────────
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.35)
    axes = [fig.add_subplot(gs[i // 4, i % 4]) for i in range(N_PARTICIPANTS)]

    all_f_flat = np.concatenate(all_f_stats)
    ymax = np.percentile(all_f_flat, 99) * 1.15
    cmap_p = plt.cm.tab10

    for i, (ax, f_stats) in enumerate(zip(axes, all_f_stats)):
        ax.plot(periods, f_stats, lw=1.2, color=cmap_p(i / 10))
        ax.axvline(dominant_info[i]["dominant_T"], color="black", lw=1.0, ls=":")
        ax.set_title(f"P{i + 1}", fontsize=10)
        ax.set_xlim(MIN_PERIOD, MAX_PERIOD)
        ax.set_ylim(0, ymax)
        ax.set_xticks([100, 200, 300])
        if i % 4 == 0:
            ax.set_ylabel("F-statistic", fontsize=9)
        if i >= 8:
            ax.set_xlabel("Period (min)", fontsize=9)
        ax.tick_params(labelsize=8)

    ax_legend = fig.add_subplot(gs[2, 3])
    ax_legend.axis("off")
    ax_legend.legend(handles=[
        plt.Line2D([0], [0], color="black", ls=":", lw=1.2, label="Dominant period"),
    ], loc="center", fontsize=9, framealpha=0.9)

    fig.suptitle(
        f"{hormone}: cosinor periodogram (F-statistic vs. period)\n"
        f"Vertical dashed line = dominant ultradian period  "
        f"(sinc-detrended at {CUTOFF} min)",
        fontsize=12, y=1.01,
    )
    plt.savefig(f"{OUTPUT_DIR}/cosinor_periodogram_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/cosinor_periodogram_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved periodogram figure.")

    # ── Figure 2 : Dominant period summary ────────────────────────────────────
    dom_periods = [d["dominant_T"] for d in dominant_info]
    sig_flags   = [d["significant"] for d in dominant_info]

    fig, ax = plt.subplots(figsize=(5, 6))
    bp = ax.boxplot([dom_periods], patch_artist=True, widths=0.4,
                    medianprops=dict(color="black", lw=2.5))
    bp["boxes"][0].set_facecolor("steelblue")
    bp["boxes"][0].set_alpha(0.55)

    rng = np.random.default_rng(0)
    for j, (val, sig) in enumerate(zip(dom_periods, sig_flags)):
        jitter = rng.uniform(-0.12, 0.12)
        ax.scatter(1 + jitter, val,
                   color="black" if sig else "#888",
                   s=40, zorder=3,
                   marker="o" if sig else "x",
                   label="significant" if (sig and j == 0) else
                         ("n.s." if (not sig and j == 0) else ""))

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if by_label:
        ax.legend(by_label.values(), by_label.keys(), fontsize=10)

    ax.set_xticks([1])
    ax.set_xticklabels([hormone], fontsize=12)
    ax.set_ylabel("Dominant ultradian period (min)", fontsize=12)
    ax.set_title(
        f"{hormone}: dominant ultradian period\n"
        f"median = {np.median(dom_periods):.1f} min  "
        f"(n = {N_PARTICIPANTS})",
        fontsize=11,
    )
    ax.set_ylim(ULTRADIAN_LO - 10, ULTRADIAN_HI + 10)
    ax.axhspan(ULTRADIAN_LO, ULTRADIAN_HI, alpha=0.06, color="steelblue")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/cosinor_dominant_period_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/cosinor_dominant_period_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved dominant period figure.")

    # ── Figure 3 : Group-median F-statistic periodogram ───────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))
    stack    = np.vstack(all_f_stats)
    median_f = np.median(stack, axis=0)
    q25_f    = np.percentile(stack, 25, axis=0)
    q75_f    = np.percentile(stack, 75, axis=0)

    ax.plot(periods, median_f, color="steelblue", lw=2, label="Median F-statistic")
    ax.fill_between(periods, q25_f, q75_f,
                    color="steelblue", alpha=0.25, label="IQR (25th–75th percentile)")
    ax.set_xlabel("Period (min)", fontsize=13)
    ax.set_ylabel("F-statistic", fontsize=13)
    ax.set_title(
        f"{hormone}: group-median cosinor periodogram  (n = {N_PARTICIPANTS})\n"
        "Shaded region = IQR across participants",
        fontsize=12,
    )
    ax.set_xlim(MIN_PERIOD, MAX_PERIOD)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/cosinor_group_mean_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/cosinor_group_mean_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved group-mean F-statistic figure.")

    print(f"\n  Summary ({hormone}):")
    print(f"  Median dominant period : {np.median(dom_periods):.1f} min")
    print(f"  Mean   dominant period : {np.mean(dom_periods):.1f} min")
    print(f"  Std                    : {np.std(dom_periods):.1f} min")
    n_sig = sum(d["significant"] for d in dominant_info)
    print(f"  Participants with significant ultradian rhythm "
          f"(α={ALPHA}): {n_sig}/{N_PARTICIPANTS}")

print("\nDone.")
