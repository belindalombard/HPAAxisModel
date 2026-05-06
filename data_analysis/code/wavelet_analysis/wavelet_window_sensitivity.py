"""
Wavelet analysis sensitivity to sliding window size — Figure S9.

Reviewer concern: the optimised sliding window was chosen to highlight ultradian
periodicity, which is circular when the wavelet is the method used to *detect*
ultradian rhythms.

This script re-runs the wavelet analysis using a range of fixed window sizes
(including the half-day window of 720 min suggested as a more principled choice)
and shows that the dominant ultradian period is robust to this choice.

Note on period range:  The original data were sampled at 10/20-min intervals and
interpolated to 1 min.  Interpolation does not add information below the
original Nyquist limit (≥ 40 min for 20-min sampling; ≥ 20 min for 10-min
sampling).  Period analysis is therefore restricted to 40–300 min, which still
spans the full ultradian range of interest (100–250 min).

Outputs (in data_analysis/code/wavelet_analysis/output/sensitivity_analysis/):
  sensitivity_spectra_Cortisol.pdf       – group-mean power spectra for each window size
  sensitivity_spectra_participants_Cortisol.pdf  – per-participant power spectra (all windows)
  sensitivity_boxplot_Cortisol.pdf       – dominant period distribution per window size  ← Figure S9
  sensitivity_individual_Cortisol.pdf    – per-participant dominant periods across window sizes
  sensitivity_summary_Cortisol.csv       – participant-level dominant periods (table)
  (same for ACTH)

Usage:
    python data_analysis/code/wavelet_analysis/wavelet_window_sensitivity.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pyboat import WAnalyzer

# ── configuration ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_BASE_DIR   = os.path.dirname(os.path.dirname(_SCRIPT_DIR))  # data_analysis/

DATA_FILE    = os.path.join(_BASE_DIR, "data", "data_per_participant_pyboat.csv")
OUTPUT_DIR   = os.path.join(_SCRIPT_DIR, "output", "sensitivity_analysis")
N_PARTICIPANTS = 10

DT            = 1           # data resolution (1 min, after interpolation)
CUTOFF        = 1440        # sinc-smooth cutoff to remove circadian trend (24 h)
MIN_PERIOD    = 50          # min period to analyse (min) – Nyquist-safe
MAX_PERIOD    = 250         # max period (min)
N_PERIODS     = 500

ULTRADIAN_LO  = 50         # ultradian window for dominant-period extraction
ULTRADIAN_HI  = 250

# Window sizes to test (min); 720 = half a day (principled choice)
WINDOW_SIZES  = [480, 600, 720, 840, 960]
WINDOW_LABELS = ["480", "600", "720\n(half-day)", "840", "960"]

# Reference window used to add the unnormalised case to the plots
NORM_REF_WIN  = 720

HORMONES = ["Cortisol", "ACTH"]

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


def run_wavelet(signal, window_size, periods, dt=DT, cutoff=CUTOFF, normalise=True):
    """Detrend, optionally normalise, and compute wavelet power spectrum.

    Parameters
    ----------
    normalise : bool
        If True (default), apply sliding-window amplitude normalisation before
        computing the spectrum.  If False, fit the detrended signal directly.

    Returns
    -------
    modulus : ndarray, shape (n_periods, n_times)
        Wavelet power at each period and time point.
    """
    wAn = WAnalyzer(periods, dt, time_unit_label="min")
    trend     = wAn.sinc_smooth(signal, T_c=cutoff)
    detrended = signal - trend
    if normalise:
        sig_to_fit = wAn.normalize_amplitude(detrended, window_size=window_size)
    else:
        sig_to_fit = detrended
    modulus, _ = wAn.compute_spectrum(sig_to_fit, do_plot=False)
    return modulus


def dominant_ultradian_period(modulus, periods,
                               lo=ULTRADIAN_LO, hi=ULTRADIAN_HI):
    """Mean power over time, then argmax within the ultradian band."""
    mask  = (periods >= lo) & (periods <= hi)
    power = modulus[mask].mean(axis=1)        # time-average
    return periods[mask][np.argmax(power)]


def mean_power_spectrum(modulus):
    """Return time-averaged power for each period."""
    return modulus.mean(axis=1)


# ── main analysis ──────────────────────────────────────────────────────────────
periods  = np.linspace(MIN_PERIOD, MAX_PERIOD, N_PERIODS)
df       = pd.read_csv(DATA_FILE)

for hormone in HORMONES:
    print(f"\n{'='*60}\n  {hormone}\n{'='*60}")

    # Store results: {label: {dominant_periods, mean_spectra}}
    results = {}

    # --- fixed window sizes ---
    for win in WINDOW_SIZES:
        dom_periods   = []
        mean_spectra  = []
        for p in range(1, N_PARTICIPANTS + 1):
            signal  = interpolate_nans(df[f"{hormone}_{p}"].values)
            modulus = run_wavelet(signal, win, periods)
            dom_periods.append(dominant_ultradian_period(modulus, periods))
            mean_spectra.append(mean_power_spectrum(modulus))
        results[win] = dict(dominant_periods=dom_periods, mean_spectra=mean_spectra)
        print(f"  window={win:4d} min  |  median dominant period = "
              f"{np.median(dom_periods):.1f} min  "
              f"(range {min(dom_periods):.0f}–{max(dom_periods):.0f} min)")

    # ── Unnormalised reference (NORM_REF_WIN, no amplitude normalisation) ─────
    dom_nn  = []
    spec_nn = []
    for p in range(1, N_PARTICIPANTS + 1):
        signal  = interpolate_nans(df[f"{hormone}_{p}"].values)
        modulus = run_wavelet(signal, NORM_REF_WIN, periods, normalise=False)
        dom_nn.append(dominant_ultradian_period(modulus, periods))
        spec_nn.append(mean_power_spectrum(modulus))
    results["no_norm"] = dict(dominant_periods=dom_nn, mean_spectra=spec_nn)
    print(f"  no-norm ({NORM_REF_WIN} min)    |  median dominant period = "
          f"{np.median(dom_nn):.1f} min")

    # ── Figure 1 : Mean power spectra ─────────────────────────────────────────
    cmap   = plt.cm.viridis
    colors = [cmap(v) for v in np.linspace(0.15, 0.85, len(WINDOW_SIZES))]

    fig, ax = plt.subplots(figsize=(8, 5))

    for i, win in enumerate(WINDOW_SIZES):
        stack     = np.vstack(results[win]["mean_spectra"])
        mean_spec = stack.mean(axis=0)
        sem_spec  = stack.std(axis=0) / np.sqrt(N_PARTICIPANTS)
        ax.plot(periods, mean_spec, color=colors[i], lw=1.8, label=f"{win} min")
        ax.fill_between(periods,
                        mean_spec - sem_spec, mean_spec + sem_spec,
                        color=colors[i], alpha=0.15)

    stack_nn = np.vstack(results["no_norm"]["mean_spectra"])
    mean_nn  = stack_nn.mean(axis=0)
    sem_nn   = stack_nn.std(axis=0) / np.sqrt(N_PARTICIPANTS)
    ax.plot(periods, mean_nn, color="crimson", lw=2.0, ls="--",
            label="No normalisation")
    ax.fill_between(periods, mean_nn - sem_nn, mean_nn + sem_nn,
                    color="crimson", alpha=0.12)

    ax.axvspan(ULTRADIAN_LO, ULTRADIAN_HI, alpha=0.07, color="gray",
               label=f"Ultradian range ({ULTRADIAN_LO}–{ULTRADIAN_HI} min)")
    ax.set_xlabel("Period (min)", fontsize=13)
    ax.set_ylabel("Mean wavelet power (a.u.)", fontsize=13)
    ax.set_title(
        f"{hormone}: sensitivity of mean power spectrum\n"
        f"to sliding window size  (n = {N_PARTICIPANTS})",
        fontsize=13,
    )
    ax.legend(fontsize=10, framealpha=0.9)
    ax.set_xlim(MIN_PERIOD, MAX_PERIOD)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_spectra_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_spectra_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved power spectra figure.")

    # ── Figure 2 : Per-participant power spectra grid ──────────────────────────
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.35)
    axes_pp = [fig.add_subplot(gs[i // 4, i % 4]) for i in range(N_PARTICIPANTS)]

    all_spectra_flat = np.concatenate(
        [results[win]["mean_spectra"][p]
         for win in WINDOW_SIZES for p in range(N_PARTICIPANTS)]
    )
    ymax_pp = np.percentile(all_spectra_flat, 99) * 1.15

    for p_idx, ax_pp in enumerate(axes_pp):
        for i, win in enumerate(WINDOW_SIZES):
            ax_pp.plot(periods, results[win]["mean_spectra"][p_idx],
                       color=colors[i], lw=1.2, alpha=0.85,
                       label=f"{win} min" if p_idx == 0 else None)
        ax_pp.plot(periods, results["no_norm"]["mean_spectra"][p_idx],
                   color="crimson", lw=1.2, ls="--", alpha=0.85,
                   label="No normalisation" if p_idx == 0 else None)

        ref_win = 720 if 720 in WINDOW_SIZES else WINDOW_SIZES[len(WINDOW_SIZES) // 2]
        dom_T   = results[ref_win]["dominant_periods"][p_idx]
        ax_pp.axvline(dom_T, color="black", lw=1.0, ls=":")

        ax_pp.set_title(f"P{p_idx + 1}", fontsize=10)
        ax_pp.set_xlim(MIN_PERIOD, MAX_PERIOD)
        ax_pp.set_ylim(0, ymax_pp)
        if p_idx % 4 == 0:
            ax_pp.set_ylabel("Mean power (a.u.)", fontsize=9)
        if p_idx >= 8:
            ax_pp.set_xlabel("Period (min)", fontsize=9)
        ax_pp.tick_params(labelsize=8)

    ax_leg = fig.add_subplot(gs[2, 3])
    ax_leg.axis("off")
    legend_handles = [
        plt.Line2D([0], [0], color=colors[i], lw=1.8, label=f"{win} min")
        for i, win in enumerate(WINDOW_SIZES)
    ] + [
        plt.Line2D([0], [0], color="crimson", lw=1.2, ls="--",
                   label="No normalisation"),
        plt.Line2D([0], [0], color="black", lw=1.0, ls=":",
                   label=f"Dominant\nperiod ({ref_win} min)"),
    ]
    ax_leg.legend(handles=legend_handles, loc="center", fontsize=9, framealpha=0.9)

    fig.suptitle(
        f"{hormone}: mean wavelet power spectrum per participant\n"
        f"Colours = window size; dashed line = dominant period at {ref_win} min window",
        fontsize=12, y=1.01,
    )
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_spectra_participants_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_spectra_participants_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved per-participant spectra figure.")

    # ── Figure 3 : Boxplot of dominant periods ─────────────────────────────────
    all_keys   = WINDOW_SIZES + ["no_norm"]
    all_labels = WINDOW_LABELS + ["No\nnormalisation"]
    all_data   = [results[k]["dominant_periods"] for k in all_keys]

    fig, ax = plt.subplots(figsize=(9, 5))
    bp = ax.boxplot(all_data, patch_artist=True, notch=False,
                    medianprops=dict(color="black", lw=2))

    box_colors = list(colors) + ["crimson"]
    for patch, c in zip(bp["boxes"], box_colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
    bp["boxes"][-1].set_hatch("///")

    rng = np.random.default_rng(42)
    for i, vals in enumerate(all_data):
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   color="black", s=18, zorder=3, alpha=0.6)

    ax.axvline(len(WINDOW_SIZES) + 0.5, color="gray", lw=0.8, ls=":")
    ax.set_xticks(range(1, len(all_keys) + 1))
    ax.set_xticklabels(all_labels, fontsize=11)
    ax.set_xlabel("Sliding window size (min)", fontsize=13)
    ax.set_ylabel("Dominant ultradian period (min)", fontsize=13)
    ax.set_title(
        f"{hormone}: dominant ultradian period vs. window size "
        f"(n = {N_PARTICIPANTS})\nUltradian band: {ULTRADIAN_LO}–{ULTRADIAN_HI} min",
        fontsize=12,
    )
    ax.set_ylim(ULTRADIAN_LO - 10, ULTRADIAN_HI + 10)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_boxplot_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_boxplot_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved boxplot figure.")

    # ── Figure 4 : Individual participant lines ────────────────────────────────
    all_cases_4  = WINDOW_SIZES + ["no_norm"]
    all_labels_4 = WINDOW_LABELS + ["No\nnormalisation"]
    x_pos        = np.arange(1, len(all_cases_4) + 1)
    fig, ax = plt.subplots(figsize=(9, 5))
    cmap_ind = plt.cm.tab10

    for p in range(N_PARTICIPANTS):
        vals = [results[k]["dominant_periods"][p] for k in all_cases_4]
        ax.plot(x_pos, vals, color=cmap_ind(p / max(N_PARTICIPANTS, 10)),
                marker="o", markersize=5, lw=1.2, alpha=0.75,
                label=f"P{p + 1}")

    median_vals = [np.median(results[k]["dominant_periods"]) for k in all_cases_4]
    ax.plot(x_pos, median_vals, color="black", marker="D",
            markersize=7, lw=2.2, zorder=5, label="Median")

    ax.axvline(len(WINDOW_SIZES) + 0.5, color="gray", lw=0.8, ls=":")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(all_labels_4, fontsize=11)
    ax.set_xlabel("Sliding window size (min)", fontsize=13)
    ax.set_ylabel("Dominant ultradian period (min)", fontsize=13)
    ax.set_title(
        f"{hormone}: individual dominant periods across window sizes "
        f"(n = {N_PARTICIPANTS})\nUltradian band: {ULTRADIAN_LO}–{ULTRADIAN_HI} min",
        fontsize=12,
    )
    ax.set_ylim(ULTRADIAN_LO - 10, ULTRADIAN_HI + 10)
    ax.legend(fontsize=9, ncol=2, framealpha=0.9, loc="upper right")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_individual_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/sensitivity_individual_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved individual participant figure.")

    # ── Summary table ─────────────────────────────────────────────────────────
    rows = []
    for p in range(N_PARTICIPANTS):
        row = {"participant": p + 1}
        for win in WINDOW_SIZES:
            row[f"window_{win}min"] = round(results[win]["dominant_periods"][p], 1)
        row[f"no_norm_{NORM_REF_WIN}min"] = round(results["no_norm"]["dominant_periods"][p], 1)
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.loc["median"] = ["median"] + \
        [round(np.median(results[win]["dominant_periods"]), 1) for win in WINDOW_SIZES] + \
        [round(np.median(results["no_norm"]["dominant_periods"]), 1)]
    summary.loc["std"] = ["std"] + \
        [round(np.std(results[win]["dominant_periods"]), 1) for win in WINDOW_SIZES] + \
        [round(np.std(results["no_norm"]["dominant_periods"]), 1)]

    summary.to_csv(f"{OUTPUT_DIR}/sensitivity_summary_{hormone}.csv", index=False)
    print(f"  Saved summary CSV.")

print("\nDone.")
