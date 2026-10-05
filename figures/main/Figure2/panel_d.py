"""Panel D: accepted between- and within-unit associations and 95% CIs."""

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from shared_style import DATA_DIR, PALETTE as P, save_panel, scale_panel_fonts, style_axes


def build_panel_d():
    data = pd.read_csv(DATA_DIR / "hierarchy_associations.csv")
    systems = ["Cardiopedia organoids", "hiPSC-CM pharmacology"]
    yloc = {"Cardiopedia organoids": 1, "hiPSC-CM pharmacology": 0}
    colors = {"between": P["blue"], "within": P["mid_grey"]}
    markers = {"between": "o", "within": "o"}

    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=False)
    # Preserve the axes footprint and make enough canvas room for the doubled
    # x-axis label without clipping it.
    fig.subplots_adjust(left=0.289, right=0.824, bottom=0.20, top=0.85)
    ax.axvline(0, color=P["mid_grey"], linewidth=0.9, linestyle=(0, (3, 3)), zorder=1)
    for system in systems:
        group = data.loc[data["system"] == system].set_index("hierarchy")
        y = yloc[system]
        between = group.loc["between"]
        within = group.loc["within"]
        for kind, row in [("between", between), ("within", within)]:
            rho = float(row["rho"])
            lo = float(row["ci_lo"])
            hi = float(row["ci_hi"])
            ax.errorbar(rho, y, xerr=[[rho - lo], [hi - rho]], fmt=markers[kind],
                        markersize=7.2, color=colors[kind], ecolor=colors[kind],
                        markerfacecolor=colors[kind] if kind == "between" else P["white"],
                        markeredgecolor=colors[kind],
                        elinewidth=1.5, capsize=3.5, zorder=3)
            label = f"− {abs(rho):.2f}" if rho < 0 else f"{rho:.2f}"
            ax.annotate(label, (rho, y), xytext=(0, 7),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=8.2, fontweight="bold", color=colors[kind],
                        bbox={"facecolor": P["white"], "edgecolor": "none", "pad": 0.15})

    ax.set_xlim(-0.9, 0.8)
    ax.set_ylim(-0.55, 1.55)
    ax.set_xticks([-0.8, -0.4, 0, 0.4, 0.8])
    ax.set_yticks([1, 0])
    ax.set_yticklabels(["Cardiopedia organoids\n24-h Force",
                         "hiPSC-CM pharmacology\nCaT amplitude"])
    ax.set_ylabel("")
    ax.set_xlabel("Molecular–functional association\n(Spearman ρ)", labelpad=7)
    style_axes(ax)
    ax.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", markersize=6.4,
               markerfacecolor=colors["between"], markeredgecolor=colors["between"],
               label="Between perturbations"),
        Line2D([], [], marker="o", linestyle="none", markersize=6.4,
               markerfacecolor=P["white"], markeredgecolor=colors["within"],
               label="Within perturbations"),
    ], loc="lower right", bbox_to_anchor=(1.0, 1.03),
       frameon=True, facecolor=P["white"], edgecolor="#CCD2D7",
       framealpha=0.96, fontsize=7.7, handlelength=1.0, borderpad=0.5, labelspacing=0.35)
    scale_panel_fonts(ax)
    save_panel(fig, "Figure2_panel_d")


if __name__ == "__main__":
    build_panel_d()
