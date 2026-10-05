"""Panel C: the two accepted engineered-tissue frozen-score contrasts."""

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from shared_style import DATA_DIR, PALETTE as P, save_panel, scale_panel_fonts, style_axes


def build_panel_c():
    data = pd.read_csv(DATA_DIR / "engineered_tissue_contrasts.csv")
    fig, ax = plt.subplots(figsize=(6.5, 5.8), constrained_layout=False)
    fig.subplots_adjust(left=0.37, right=0.98, bottom=0.31, top=0.93)
    ax.axvline(0, color=P["mid_grey"], linewidth=1.0, linestyle=(0, (3, 3)), zorder=1)
    colors = ["#6E7780", P["blue"]]
    markers = ["o", "s"]
    y = [1, 0]
    for i, row in data.iterrows():
        estimate = float(row["mean_score_difference"])
        lo = float(row["ci_lo"])
        hi = float(row["ci_hi"])
        ax.errorbar(estimate, y[i], xerr=[[estimate - lo], [hi - estimate]], fmt=markers[i],
                    markersize=8, color=colors[i], ecolor=colors[i], elinewidth=1.75,
                    capsize=4, zorder=3)
    ax.set_xlim(-800, 600)
    ax.set_ylim(-0.55, 1.55)
    ax.set_xticks([-800, -600, -400, -200, 0, 200, 400, 600])
    ax.set_yticks(y)
    ax.set_yticklabels(["Tachypaced EHT\n(day 3; n = 2/2)",
                         "Dystrophin-deficient EHT\n(n = 3/3)"])
    ax.set_xlabel("Frozen projection difference\n(worse − better mechanics)", labelpad=7)
    ax.set_ylabel("")
    style_axes(ax)
    ax.legend(handles=[
        Line2D([], [], color=colors[0], marker="o", markersize=6, linestyle="none",
               label="Tachypaced EHT"),
        Line2D([], [], color=colors[1], marker="s", markersize=6, linestyle="none",
               label="Dystrophin-deficient EHT"),
    ], loc="upper right", frameon=True, facecolor=P["white"], edgecolor="#CCD2D7",
       framealpha=0.96, fontsize=7.5, handlelength=1.5, labelspacing=0.35)
    scale_panel_fonts(ax)
    save_panel(fig, "Figure2_panel_c")


if __name__ == "__main__":
    build_panel_c()
