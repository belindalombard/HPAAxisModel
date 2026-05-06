"""
Best-model cosinor analysis: single vs. multi-component ultradian rhythms — Figure S10.

Uses CosinorPy's cosinor.fit_group to scan ultradian periods and fit models
with 1 and 2 harmonic components for each participant.

Model selection is a two-stage hybrid:
  1. cosinor.get_best_fits (R²_adj, per nc) finds the best-fitting period
     separately for each number of harmonic components.
  2. A sequential nested F-test (α = ALPHA_NESTED) compares nc=1 vs nc=2,
     only accepting the more complex model if the extra harmonic explains
     significantly more residual variance.
     This guards against the large-n over-fitting that occurs when using
     R²_adj or p-value directly to select nc (with n≈1431, α=0.05 is too
     permissive; a stricter threshold is used instead).

A model with n_components=1 is a pure sinusoid.
A model with n_components=2 adds a 2nd harmonic, capturing asymmetric
waveform shape.

Pre-processing:
  1. Linear NaN interpolation
  2. Sinc-filter detrending at T_c = 1440 min (removes circadian trend).

NOTE: CosinorPy's cosinor.fit_me uses np.round_ which was removed in NumPy 2.0.
      The monkey-patch `np.round_ = np.round` at the top of this file fixes this
      and must appear BEFORE any CosinorPy import.

Outputs (in data_analysis/code/wavelet_analysis/output/cosinor_best_models/):
  best_models_Cortisol.pdf      – 3×4 per-participant grid: signal + best fit  ← Figure S10
  periodogram_Cortisol.pdf      – 3×4 per-participant R²_adj periodograms  ← Figure S11
  best_models_summary.pdf       – summary: best n_components and period per participant
  best_models_Cortisol.csv      – per-participant best model statistics
  best_models_all.csv           – combined table across both hormones
  (same for ACTH)

Usage:
    python data_analysis/code/wavelet_analysis/cosinor_best_models.py
"""

import numpy as np
np.round_ = np.round   # NumPy 2.0 compatibility – MUST be before CosinorPy import

import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy import stats as scipy_stats
from CosinorPy import cosinor
from pyboat import core as pyboat_core

# ── configuration ──────────────────────────────────────────────────────────────
_SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
_BASE_DIR      = os.path.dirname(os.path.dirname(_SCRIPT_DIR))  # data_analysis/

DATA_FILE      = os.path.join(_BASE_DIR, "data", "data_per_participant_pyboat.csv")
OUTPUT_DIR     = os.path.join(_SCRIPT_DIR, "output", "cosinor_best_models")
N_PARTICIPANTS = 10

DT             = 1          # 1-min resolution
CUTOFF         = 1440       # sinc-filter cutoff for display detrending (min)

ULTRADIAN_LO   = 50
ULTRADIAN_HI   = 250
PERIOD_STEP    = 1          # step for period scan (min)
N_COMPONENTS   = [1, 2]    # harmonic model complexities to compare

ALPHA          = 0.05       # omnibus significance threshold for final model
ALPHA_NESTED   = 1e-5      # nested F-test threshold for nc selection
HORMONES       = ["Cortisol", "ACTH"]

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── helpers ────────────────────────────────────────────────────────────────────
def interpolate_nans(signal):
    s = np.array(signal, dtype=float)
    mask = np.isnan(s)
    if mask.any():
        idx = np.arange(len(s))
        s[mask] = np.interp(idx[mask], idx[~mask], s[~mask])
    return s


def sinc_detrend(signal, cutoff=CUTOFF, dt=DT):
    """Remove circadian trend via sinc-filter."""
    trend = pyboat_core.sinc_smooth(signal, T_c=cutoff, dt=dt)
    return signal - trend


def cosinor_predict(fitted_model, t, period, n_components):
    """Evaluate a fitted CosinorPy statsmodels model at arbitrary timepoints.

    CosinorPy's fit_me only returns Y_test over ~1 period; this rebuilds the
    design matrix (matching CosinorPy's convention) and predicts across all of t.
    """
    cols = [np.ones(len(t))]
    for k in range(1, n_components + 1):
        cols.append(np.cos(2 * np.pi * k * t / period))
        cols.append(np.sin(2 * np.pi * k * t / period))
    return fitted_model.predict(np.column_stack(cols))


# Colours: one per n_components value
COMP_COLORS = {1: "steelblue", 2: "darkorange"}
COMP_LABELS = {1: "1 component\n(simple sinusoid)",
               2: "2 components\n(+ 2nd harmonic)"}

# ── data loading ───────────────────────────────────────────────────────────────
df_data      = pd.read_csv(DATA_FILE)
scan_periods = np.arange(ULTRADIAN_LO, ULTRADIAN_HI + PERIOD_STEP, PERIOD_STEP,
                         dtype=float)

# ── per-hormone analysis ────────────────────────────────────────────────────────
summary_rows = []

for hormone in HORMONES:
    print(f"\n{'='*60}\n  {hormone}\n{'='*60}")

    signals   = {}
    detrended = {}
    t_arr     = None
    for p in range(1, N_PARTICIPANTS + 1):
        signals[p]   = interpolate_nans(df_data[f"{hormone}_{p}"].values)
        detrended[p] = sinc_detrend(signals[p])
        if t_arr is None:
            t_arr = np.arange(len(signals[p]), dtype=float)

    records = []
    for p in range(1, N_PARTICIPANTS + 1):
        for ti, yi in zip(t_arr, detrended[p]):
            records.append({"x": ti, "y": yi, "test": f"P{p}"})
    df_long = pd.DataFrame(records)

    print(f"  Scanning {len(scan_periods)} periods × {N_COMPONENTS} components …")
    all_fits = []
    for per in scan_periods:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = cosinor.fit_group(df_long, n_components=N_COMPONENTS,
                                  period=[float(per)], plot=False)
        r["period"] = per
        all_fits.append(r)
    df_fits = pd.concat(all_fits, ignore_index=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df_best_fits = cosinor.get_best_fits(
            df_fits, criterium='R2_adj', reverse=False, n_components=N_COMPONENTS
        )

    plot_info = {}
    n_obs = len(t_arr)

    for p in range(1, N_PARTICIPANTS + 1):
        p_label = f"P{p}"
        rows_p  = df_best_fits[df_best_fits["test"] == p_label].sort_values("n_components")
        best_row = rows_p[rows_p["n_components"] == N_COMPONENTS[0]].iloc[0]

        for nc in N_COMPONENTS[1:]:
            cand = rows_p[rows_p["n_components"] == nc]
            if cand.empty:
                continue
            cand_row = cand.iloc[0]
            rss_red  = float(best_row["RSS"])
            rss_full = float(cand_row["RSS"])
            df_red   = n_obs - (int(best_row["n_components"]) * 2 + 1)
            df_full  = n_obs - (nc * 2 + 1)
            if rss_red > rss_full and df_red > df_full:
                F_stat = ((rss_red - rss_full) / (df_red - df_full)) / (rss_full / df_full)
                p_F    = 1 - scipy_stats.f.cdf(F_stat, df_red - df_full, df_full)
                if p_F < ALPHA_NESTED:
                    best_row = cand_row

        best_T  = float(best_row["period"])
        best_nc = int(best_row["n_components"])
        r2      = float(best_row["R2_adj"])
        p_val   = float(best_row["p"])
        sig     = p_val < ALPHA

        plot_info[p] = dict(best_T=best_T, best_nc=best_nc, r2=r2, p_model=p_val)
        summary_rows.append(dict(
            hormone=hormone, participant=p,
            best_period=best_T, n_components=best_nc,
            R2_adj=r2, p=p_val, significant=sig,
        ))

    print(f"\n  {'P':>4}  {'T':>6}  {'nc':>4}  {'R²_adj':>8}  {'p':>10}  {'sig':>5}")
    for p in range(1, N_PARTICIPANTS + 1):
        info    = plot_info[p]
        sig_str = '✓' if info['p_model'] < ALPHA else 'n.s.'
        print(f"  P{p:2d}  {info['best_T']:6.0f}  {info['best_nc']:4d}  "
              f"{info['r2']:8.3f}  {info['p_model']:10.4f}  {sig_str:>5}")

    # ── Figure 1 : Per-participant grid (signal + best-fit) ───────────────────
    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.35)
    axes = [fig.add_subplot(gs[i // 4, i % 4]) for i in range(N_PARTICIPANTS)]

    for i in range(N_PARTICIPANTS):
        p_num   = i + 1
        ax      = axes[i]
        info    = plot_info[p_num]
        best_T  = info["best_T"]
        best_nc = info["best_nc"]
        r2      = info["r2"]
        p_model = info["p_model"]
        sig     = p_model < ALPHA
        color   = COMP_COLORS[best_nc]

        sig_full = detrended[p_num]
        ax.plot(t_arr, sig_full, color="lightgray", lw=0.6, zorder=1)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit_result = cosinor.fit_me(
                t_arr, sig_full,
                n_components=best_nc, period=best_T,
                plot=False, return_model=True,
            )
        t_dense = np.linspace(t_arr[0], t_arr[-1], 3000)
        Y_pred  = cosinor_predict(fit_result[0], t_dense, best_T, best_nc)
        ax.plot(t_dense, Y_pred, color=color, lw=1.8, zorder=3)

        ax.set_title(
            f"P{p_num}  T={best_T:.0f} min  nc={best_nc}"
            f"  R²={r2:.2f}  {'✓' if sig else 'n.s.'}",
            fontsize=9,
        )
        ax.set_xlim(t_arr[0], t_arr[-1])
        ax.tick_params(labelsize=7)
        if i % 4 == 0:
            ax.set_ylabel("Signal", fontsize=8)
        if i >= 8:
            ax.set_xlabel("Time (min)", fontsize=8)

    ax_leg = fig.add_subplot(gs[2, 3])
    ax_leg.axis("off")
    patches = [mpatches.Patch(color=COMP_COLORS[nc], label=COMP_LABELS[nc])
               for nc in N_COMPONENTS]
    ax_leg.legend(handles=patches, loc="center", fontsize=8.5, framealpha=0.9)

    fig.suptitle(
        f"{hormone}: best-fit cosinor model per participant\n"
        f"Period selected by R²_adj; nc by nested F-test (α={ALPHA_NESTED})",
        fontsize=12, y=1.01,
    )
    plt.savefig(f"{OUTPUT_DIR}/best_models_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/best_models_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved per-participant best-model figure.")

    # ── Figure 2 : Per-participant R²_adj periodograms ────────────────────────
    fig_pg  = plt.figure(figsize=(14, 10))
    gs_pg   = gridspec.GridSpec(3, 4, figure=fig_pg, hspace=0.55, wspace=0.35)
    axes_pg = [fig_pg.add_subplot(gs_pg[i // 4, i % 4]) for i in range(N_PARTICIPANTS)]

    for i in range(N_PARTICIPANTS):
        p_num  = i + 1
        ax_pg  = axes_pg[i]
        info   = plot_info[p_num]
        best_T = info["best_T"]
        best_nc = info["best_nc"]

        sub_p  = df_fits[df_fits["test"] == f"P{p_num}"]
        sub_nc = sub_p[sub_p["n_components"] == best_nc].sort_values("period")
        ax_pg.plot(sub_nc["period"].values, sub_nc["R2_adj"].values,
                   color=COMP_COLORS[best_nc], lw=1.4)
        ax_pg.axvline(best_T, color=COMP_COLORS[best_nc], lw=1.6,
                      linestyle="--", zorder=3, alpha=0.85)

        ax_pg.set_title(f"P{p_num}  T={best_T:.0f} min  nc={best_nc}", fontsize=9)
        ax_pg.set_xlim(ULTRADIAN_LO, ULTRADIAN_HI)
        ax_pg.tick_params(labelsize=7)
        if i % 4 == 0:
            ax_pg.set_ylabel("R²_adj", fontsize=8)
        if i >= 8:
            ax_pg.set_xlabel("Period (min)", fontsize=8)

    ax_leg_pg = fig_pg.add_subplot(gs_pg[2, 3])
    ax_leg_pg.axis("off")
    _pg_ncs = sorted({plot_info[p]["best_nc"] for p in range(1, N_PARTICIPANTS + 1)})
    line_handles = [mlines.Line2D([0], [0], color=COMP_COLORS[nc], lw=1.5, label=f"nc={nc}")
                    for nc in _pg_ncs]
    ax_leg_pg.legend(handles=line_handles, loc="center", fontsize=8.5, framealpha=0.9)

    fig_pg.suptitle(f"{hormone}: R²_adj periodogram", fontsize=12, y=1.01)
    plt.savefig(f"{OUTPUT_DIR}/periodogram_{hormone}.pdf",
                dpi=300, bbox_inches="tight")
    plt.savefig(f"{OUTPUT_DIR}/periodogram_{hormone}.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved periodogram figure.")

    summary_hormone = pd.DataFrame(
        [r for r in summary_rows if r["hormone"] == hormone]
    )
    summary_hormone.to_csv(f"{OUTPUT_DIR}/best_models_{hormone}.csv", index=False)
    print(f"  Saved CSV.")

# ── Cross-hormone summary figure ───────────────────────────────────────────────
summary_df = pd.DataFrame(summary_rows)

fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)

for ax, hormone in zip(axes, HORMONES):
    sub = summary_df[summary_df["hormone"] == hormone].sort_values("participant")
    participants = sub["participant"].values
    best_periods = sub["best_period"].values
    ncs          = sub["n_components"].values
    colors_bar   = [COMP_COLORS[nc] for nc in ncs]

    x_pos = np.arange(len(participants))
    bars  = ax.bar(x_pos, best_periods, color=colors_bar, alpha=0.8,
                   edgecolor="black", lw=0.7)

    for bar, nc, per in zip(bars, ncs, best_periods):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 1.5,
                f"{per:.0f} min\n(nc={nc})",
                ha="center", va="bottom", fontsize=7.5, fontweight="bold")

    for xi, (nc, per) in enumerate(zip(ncs, best_periods)):
        if nc >= 2:
            h2 = per / 2
            ax.scatter(xi, h2, marker="v", color=COMP_COLORS[nc],
                       s=70, zorder=4, edgecolors="black", linewidths=0.8)
            ax.text(xi, h2 - 2, f"{h2:.0f}", ha="center", va="top",
                    fontsize=7, color="dimgray")

    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"P{p}" for p in participants], fontsize=10)
    ax.set_xlabel("Participant", fontsize=12)
    ax.set_ylabel("Period (min)", fontsize=12)
    ax.set_title(f"{hormone}", fontsize=13)
    ax.set_ylim(35, ULTRADIAN_HI + 30)
    ax.axhspan(ULTRADIAN_LO, ULTRADIAN_HI, alpha=0.06, color="gray",
               label="Ultradian band")

_ncs_used = set(summary_df["n_components"].values)
patches = [mpatches.Patch(color=COMP_COLORS[nc], label=COMP_LABELS[nc])
           for nc in N_COMPONENTS if nc in _ncs_used]
harmonic_handles = [
    plt.scatter([], [], marker="v", color="gray", s=70, edgecolors="black",
                linewidths=0.8, label="T/2 harmonic (nc≥2)"),
]
axes[-1].legend(handles=patches + harmonic_handles, fontsize=9,
                loc="upper right", framealpha=0.9)

fig.suptitle(
    "Best-fit cosinor period and harmonics per participant\n"
    "(bar = fundamental period; markers = harmonic periods for multi-component models)",
    fontsize=12,
)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/best_models_summary.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{OUTPUT_DIR}/best_models_summary.png", dpi=150, bbox_inches="tight")
plt.close()
print("\n  Saved summary figure.")

summary_df.to_csv(f"{OUTPUT_DIR}/best_models_all.csv", index=False)
print("  Saved combined CSV.")
print("\nDone.")
