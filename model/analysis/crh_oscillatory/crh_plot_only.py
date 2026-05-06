"""
Plot the CRH drive Φ(t) = λ_a · exp(λ_s · cos(2π(t-t_s)/T_c + σ·cos(...)))
over one circadian cycle (1440 min), shaping rising phases green and falling
phases red.

Uses mean cohort parameters from model/config/base/mean_parameters.json.
Output: model/analysis/crh_oscillatory/output/figures/crh_drive.png

Usage:
    python model/analysis/crh_oscillatory/crh_plot_only.py
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from pathlib import Path

# ── path setup ─────────────────────────────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT  = _SCRIPT_DIR.parents[2]          # HPAAxisModel/
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from model.code.classes.model_class import Model  # noqa: E402

# ── load mean parameters ───────────────────────────────────────────────────────
config_path = _REPO_ROOT / 'model' / 'config' / 'base' / 'mean_parameters.json'
with open(config_path, 'r') as f:
    config = json.load(f)

parameters   = config.get('parameters', {})
fixed_params = config.get('fixed_params', {})

# ── build model and evaluate CRH ──────────────────────────────────────────────
dde_model = Model(
    fixed_params=fixed_params,
    suggested_params=parameters,
    reject=False,
    num_days=2,
)

timespan = np.linspace(0, 1440, 1440)
crh_vals = [dde_model.crh(t=t, t_s=0, lambda_a=2, lambda_s=2, sigma=0.5, T_c=1440)
            for t in timespan]

# ── figure ────────────────────────────────────────────────────────────────────
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 40

fig, ax = plt.subplots(figsize=(13, 10))
ax.plot(timespan, crh_vals, color='black', linewidth=2)

for i in range(1, len(timespan)):
    if crh_vals[i] > crh_vals[i - 1]:
        ax.fill_between(timespan[i - 1:i + 1], crh_vals[i - 1:i + 1],
                        color='green', alpha=0.1)
    elif crh_vals[i] < crh_vals[i - 1]:
        ax.fill_between(timespan[i - 1:i + 1], crh_vals[i - 1:i + 1],
                        color='red', alpha=0.1)

t0 = datetime.strptime('09:00:00', "%H:%M:%S")
ax.set_xticks([i for i in range(0, 1441, 420)])
ax.set_xticklabels([(t0 + timedelta(minutes=i)).strftime("%H:%M")
                    for i in range(0, 1441, 420)])
ax.set_ylabel('CRH (Arbitrary)')
ax.set_xlabel('Time (hours)')

plt.tight_layout()
out_dir = _SCRIPT_DIR / 'output' / 'figures'
out_dir.mkdir(parents=True, exist_ok=True)
plt.savefig(str(out_dir / 'crh_drive.png'), dpi=300, bbox_inches='tight')
print(f"Saved: {out_dir / 'crh_drive.png'}")
