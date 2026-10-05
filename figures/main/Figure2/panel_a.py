"""Panel A: 82 accepted ligand-level molecular-score/24-h Force pairs."""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import t

from shared_style import DATA_DIR, PALETTE as P, save_panel, scale_panel_fonts, style_axes


def descriptive_linear_fit(x, y, grid):
    """OLS descriptive guide and pointwise 95% CI for the mean fitted value."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    beta = np.polyfit(x, y, 1)
    fitted = np.polyval(beta, x)
    dof = len(x) - 2
    residual_variance = np.sum((y - fitted) ** 2) / dof
    sxx = np.sum((x - x.mean()) ** 2)
    estimate = np.polyval(beta, grid)
    se_mean = np.sqrt(residual_variance * (1 / len(x) + (grid - x.mean()) ** 2 / sxx))
    margin = t.ppf(0.975, dof) * se_mean
    return estimate, estimate - margin, estimate + margin


def build_panel_a(target_width_px=None):
    data = pd.read_csv(DATA_DIR / "organoid_ligand_points.csv")
    x = data["mean_frozen_score"].to_numpy(float)
    # The accepted source stores Force as a control-normalized ratio; ×100
    # displays the reference figure's percent-of-control unit.
    y = data["median_force_24h"].to_numpy(float) * 100
    xlim = (-400, 400)
    ylim = (0, 250)

    fig, ax = plt.subplots(figsize=(13.5, 4.95), constrained_layout=False)
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.20, top=0.95)
    ax.scatter(x, y, s=54, color=P["blue"], edgecolor=P["white"], linewidth=0.55,
               alpha=0.82, zorder=3, label=f"Ligand summaries (n={len(data)})")
    grid = np.linspace(xlim[0], xlim[1], 300)
    fit, lo, hi = descriptive_linear_fit(x, y, grid)
    ax.fill_between(grid, lo, hi, color=P["blue_band"], alpha=0.22, linewidth=0, zorder=1)
    ax.plot(grid, fit, color=P["blue"], linewidth=1.8, alpha=0.82, zorder=2)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xticks([-400, -200, 0, 200, 400])
    ax.set_yticks([0, 50, 100, 150, 200, 250])
    ax.set_xlabel("Frozen molecular projection score (arbitrary units)", labelpad=6)
    ax.set_ylabel("24-h Force (% of control)", labelpad=7)
    style_axes(ax)
    ax.text(0.03, 0.08, "Spearman ρ = −0.705\n82 ligand summaries (382 paired organoids)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=9.2, color="#16243A")

    legend_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=8.0,
               markerfacecolor=P["blue"], markeredgecolor=P["white"], label="Ligand summaries"),
        Line2D([], [], color=P["blue"], linewidth=1.8, alpha=0.82,
               label="Descriptive linear fit"),
        Patch(facecolor=P["blue_band"], alpha=0.22, edgecolor="none",
              label="95% CI for mean fit"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=8.2,
              handlelength=1.8, labelspacing=0.3, borderpad=0.2)
    scale_panel_fonts(ax)
    save_panel(fig, "Figure2_panel_a", target_png_width_px=target_width_px)


if __name__ == "__main__":
    build_panel_a()
