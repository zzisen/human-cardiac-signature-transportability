"""Panel B: 42 compound-level Spearman estimates and accepted pooled summary."""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from shared_style import DATA_DIR, PALETTE as P, save_panel, scale_panel_fonts, style_axes


def build_panel_b():
    compounds = pd.read_csv(DATA_DIR / "pharmacology_compound_associations.csv")
    summary = pd.read_csv(DATA_DIR / "pharmacology_primary_summary.csv").query(
        "analysis == 'Primary frozen score'"
    ).iloc[0]
    values = compounds["spearman_rho"].to_numpy(float)
    # Deterministic display jitter only separates overlapping estimates.
    jitter = ((np.arange(len(values)) * 29) % 15 - 7) * 0.009

    compound_y = 1.35
    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=False)
    # Keep the plotting area the same physical size while reserving more room
    # for the doubled x-axis label on the right.
    fig.subplots_adjust(left=0.34, right=0.832, bottom=0.22, top=0.88)
    ax.axvline(0, color=P["mid_grey"], linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)
    ax.scatter(values, compound_y + jitter, s=36, color=P["grey"], edgecolor=P["white"],
               linewidth=0.45, alpha=0.88, zorder=3)
    estimate = float(summary["rho"])
    lo = float(summary["ci_lo"])
    hi = float(summary["ci_hi"])
    ax.errorbar(estimate, 0, xerr=[[estimate - lo], [hi - estimate]], fmt="D",
                markersize=7.2, color=P["blue"], ecolor=P["blue"],
                elinewidth=1.65, capsize=4, zorder=4)
    ax.set_xlim(-0.65, 0.65)
    ax.set_ylim(-0.60, 1.95)
    ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])
    ax.set_yticks([0, compound_y])
    ax.set_yticklabels(["Pooled within-compound\n310 dose groups /\n42 compounds",
                         "Compound-specific estimates\n42 compounds"])
    ax.set_xlabel("Within-compound\nSpearman ρ", labelpad=6)
    style_axes(ax)
    ax.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markersize=6.1,
               markerfacecolor=P["grey"], markeredgecolor=P["white"], label="Compound estimates"),
        Line2D([], [], marker="D", linestyle="-", markersize=5.8,
               markerfacecolor=P["blue"], markeredgecolor=P["blue"],
               color=P["blue"], label="Pooled estimate ± 95% CI"),
    ], loc="lower right", bbox_to_anchor=(1.0, 1.03),
       frameon=True, facecolor=P["white"], edgecolor="#CCD2D7",
       framealpha=0.97, fontsize=7.8, handlelength=1.4, borderpad=0.55, labelspacing=0.45)
    scale_panel_fonts(ax)
    save_panel(fig, "Figure2_panel_b")


if __name__ == "__main__":
    build_panel_b()
