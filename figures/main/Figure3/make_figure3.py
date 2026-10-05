#!/usr/bin/env python3
"""Recreate three Figure 3 panels in the visual style of Figure_Tem/Figure3.png.

The plotted observations and accepted summary estimates are read from the
project's Figure 3 source tables. This script only changes the presentation.
Requires Python, NumPy, and Matplotlib.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch


OUT_DIR = Path(__file__).resolve().parent
ROOT = next(
    (
        parent
        for parent in (OUT_DIR, *OUT_DIR.parents)
        if (parent / "06_Figures" / "COMM_BIOL_V4" / "Main" / "Figure_3" / "data").is_dir()
    ),
    None,
)
SOURCE_DATA_DIR = OUT_DIR / "source_data"
if SOURCE_DATA_DIR.is_dir():
    # The final deliverable bundles its small input tables for standalone use.
    DATA_DIR = SOURCE_DATA_DIR
    FIGURE3_SOURCE_DIR = SOURCE_DATA_DIR
elif ROOT is not None:
    # Keep the original working copy runnable from its existing project sources.
    DATA_DIR = ROOT / "06_Figures" / "COMM_BIOL_V4" / "Main" / "Figure_3" / "data"
    FIGURE3_SOURCE_DIR = ROOT / "06_Figures" / "COMM_BIOL_V3" / "Main" / "Figure_3"
else:
    raise RuntimeError("Input tables are missing and the project root could not be located.")

# Colors sampled from the supplied Figure_Tem/Figure3.png reference.
BLUE = "#1E6CC0"
TEAL = "#199891"
RED = "#D43F28"
DARK_BLUE = "#174F8C"
TEXT = "#17212B"
GUIDE = "#BAC7D7"
TRACK = "#D0D9E3"
WHITE = "#FFFFFF"

# Accepted Figure 3a statistic reported with the project source data.
SPEARMAN_RHO = 0.2385
SPEARMAN_CI = (0.1508, 0.2674)


def read_table(filename: str, data_dir: Path = DATA_DIR) -> list[dict[str, str]]:
    with (data_dir / filename).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def padded_limits(values: np.ndarray, lo: float, hi: float, pad: float = 0.06) -> tuple[float, float]:
    span = hi - lo
    margin = span * pad if span else 1.0
    return lo - margin, hi + margin


def linear_fit_ci(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """OLS line and pointwise 95% confidence interval for the fitted mean."""
    slope, intercept = np.polyfit(x, y, 1)
    fitted = intercept + slope * x
    residuals = y - fitted
    n = len(x)
    x_mean = float(np.mean(x))
    sxx = float(np.sum((x - x_mean) ** 2))
    residual_variance = float(np.sum(residuals**2) / max(n - 2, 1))
    se_mean = np.sqrt(residual_variance * (1.0 / n + (grid - x_mean) ** 2 / sxx))
    line = intercept + slope * grid
    # With 12,929 observations, the two-sided t critical value is effectively 1.96.
    return line, line - 1.96 * se_mean, line + 1.96 * se_mean


def style_axes(ax: plt.Axes, *, top: bool = False, right: bool = False) -> None:
    ax.set_facecolor(WHITE)
    ax.tick_params(axis="both", colors=TEXT, labelsize=9, width=0.9, length=4, direction="out")
    for side in ("left", "bottom", "top", "right"):
        ax.spines[side].set_color(TEXT if (side in {"left", "bottom"} or (side == "top" and top) or (side == "right" and right)) else "none")
        ax.spines[side].set_linewidth(1.0)
    ax.grid(False)


def draw_panel_a(ax: plt.Axes, genes: list[dict[str, str]]) -> None:
    x = np.asarray([float(row["frozen_projection_weight"]) for row in genes])
    y = np.asarray([float(row["adult_hf_expression_effect"]) for row in genes])

    # Main view keeps the dense central cloud legible; the inset preserves every
    # raw gene, including the long-tail expression effects.
    x_lo, x_hi = np.quantile(x, [0.025, 0.975])
    y_lo, y_hi = np.quantile(y, [0.025, 0.975])
    x_lim = padded_limits(x, float(x_lo), float(x_hi), 0.06)
    y_lim = padded_limits(y, float(y_lo), float(y_hi), 0.06)
    grid = np.linspace(x_lim[0], x_lim[1], 250)
    line, low, high = linear_fit_ci(x, y, grid)

    ax.scatter(x, y, s=7, color=BLUE, alpha=0.24, linewidths=0, rasterized=True, zorder=2)
    ax.fill_between(grid, low, high, color=BLUE, alpha=0.17, linewidth=0, zorder=3)
    ax.plot(grid, line, color=DARK_BLUE, linewidth=1.7, zorder=4)
    ax.axhline(0, color=GUIDE, linestyle=(0, (5, 3)), linewidth=1.0, zorder=1)
    ax.axvline(0, color=GUIDE, linestyle=(0, (5, 3)), linewidth=1.0, zorder=1)
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xticks([-0.3, -0.15, 0, 0.15, 0.3])
    ax.set_yticks([-1.0, -0.5, 0, 0.5])
    ax.set_xlabel("Frozen projection weight (a.u.)", fontsize=11, color=TEXT, labelpad=8)
    ax.set_ylabel("Adult-HF expression effect", fontsize=11, color=TEXT, labelpad=8)
    style_axes(ax, top=True, right=True)

    legend_items = [
        Line2D([], [], marker="o", linestyle="none", color=BLUE, markersize=4.5, alpha=0.65, label="Genes"),
        Line2D([], [], color=DARK_BLUE, linewidth=1.7, label="Linear fit"),
        Patch(facecolor=BLUE, edgecolor="none", alpha=0.17, label="OLS mean-fit 95% CI"),
    ]
    ax.legend(handles=legend_items, loc="upper right", frameon=True, framealpha=0.92,
              facecolor=WHITE, edgecolor="#B6C5D5", fontsize=7.5, borderpad=0.45,
              handlelength=1.5, labelspacing=0.35)
    ax.text(0.50, 0.96, "Main axes: central 95%", transform=ax.transAxes, ha="left", va="top",
            fontsize=7.5, color=TEXT)
    ax.text(0.98, 0.035,
            f"Spearman $\\rho$ = {SPEARMAN_RHO:.4f}\n"
            f"95% donor-bootstrap CI {SPEARMAN_CI[0]:.4f}–{SPEARMAN_CI[1]:.4f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.2, color=TEXT)

    inset = ax.inset_axes([0.065, 0.66, 0.39, 0.29])
    inset.set_facecolor(WHITE)
    inset.scatter(x, y, s=1.7, color=BLUE, alpha=0.28, linewidths=0, rasterized=True, zorder=2)
    inset_grid = np.linspace(float(np.min(x)), float(np.max(x)), 180)
    inset_line, inset_low, inset_high = linear_fit_ci(x, y, inset_grid)
    inset.fill_between(inset_grid, inset_low, inset_high, color=BLUE, alpha=0.14, linewidth=0, zorder=3)
    inset.plot(inset_grid, inset_line, color=DARK_BLUE, linewidth=1.0, zorder=4)
    inset.axhline(0, color=GUIDE, linestyle=(0, (4, 3)), linewidth=0.75, zorder=1)
    inset.axvline(0, color=GUIDE, linestyle=(0, (4, 3)), linewidth=0.75, zorder=1)
    inset.set_xlim(*padded_limits(x, float(np.min(x)), float(np.max(x)), 0.04))
    inset.set_ylim(*padded_limits(y, float(np.min(y)), float(np.max(y)), 0.04))
    inset.set_title(f"Full range (all {len(genes):,} genes)", fontsize=7.5, color=TEXT, pad=3, weight="bold")
    inset.tick_params(axis="both", labelsize=6.3, length=2.5, width=0.7, colors=TEXT, pad=1.5)
    inset.set_xticks([round(float(np.min(x)), 1), 0, round(float(np.max(x)), 1)])
    inset.set_yticks([round(float(np.min(y)), 0), 0, round(float(np.max(y)), 0)])
    for spine in inset.spines.values():
        spine.set_color("#91A9C3")
        spine.set_linewidth(0.9)


def draw_panel_b(ax: plt.Axes, sign_rows: list[dict[str, str]], genes: list[dict[str, str]]) -> None:
    metrics = {row["metric"]: row for row in sign_rows}
    denominator = sum(as_bool(row["nonzero_sign_eligible"]) for row in genes)
    concordant_count = sum(as_bool(row["sign_concordant"]) for row in genes if as_bool(row["nonzero_sign_eligible"]))
    rows = [
        ("Sign concordance", metrics["Sign concordance"], BLUE,
         f"{concordant_count:,} / {denominator:,} shared eligible genes"),
        ("Absolute-weighted sign concordance", metrics["Absolute-weighted sign concordance"], TEAL,
         "weighted by |frozen projection weight|"),
    ]
    y_values = [1.55, 0.52]

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.02, 2.58)
    ax.set_yticks(y_values, [row[0] for row in rows])
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Sign concordance (%)", fontsize=10.5, color=TEXT, labelpad=7)
    style_axes(ax)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=9, pad=9)

    for (label, row, color, note), y in zip(rows, y_values):
        estimate = float(row["estimate"]) * 100
        lo = float(row["ci_lo"]) * 100
        hi = float(row["ci_hi"]) * 100
        track = FancyBboxPatch((0, y - 0.12), 100, 0.24,
                               boxstyle="round,pad=0.012,rounding_size=0.12",
                               linewidth=0, facecolor=TRACK, transform=ax.transData, zorder=1, clip_on=True)
        fill = FancyBboxPatch((0, y - 0.12), estimate, 0.24,
                              boxstyle="round,pad=0.012,rounding_size=0.12",
                              linewidth=0, facecolor=color, transform=ax.transData, zorder=2, clip_on=True)
        ax.add_patch(track)
        ax.add_patch(fill)
        ax.text(min(estimate + 1.8, 97), y, f"{estimate:.1f}%", ha="left", va="center",
                color=color, fontsize=10.5, weight="bold", clip_on=True)
        ax.text(1, y - 0.27, f"95% donor-bootstrap CI {lo:.1f}–{hi:.1f}% · {note}",
                ha="left", va="top", color=TEXT, fontsize=7.0)

    legend_items = [Patch(facecolor=BLUE, edgecolor="none", label="Sign concordance"),
                    Patch(facecolor=TEAL, edgecolor="none", label="Absolute-weighted sign concordance")]
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.56, 1.02), ncol=2,
              frameon=False, fontsize=7.3, handlelength=0.9, handleheight=0.8,
              columnspacing=1.1, borderaxespad=0)


def draw_panel_c(ax: plt.Axes, f3c_rows: list[dict[str, str]]) -> None:
    contrast_to_group = {
        "Ischemic cardiomyopathy": "ICM",
        "Dilated cardiomyopathy": "DCM",
    }
    estimates = {
        contrast_to_group[row["contrast"]]: row
        for row in f3c_rows
        if row["contrast"] in contrast_to_group
    }
    specs = [("ICM", RED, 1.0), ("DCM", BLUE, 0.0)]

    ax.set_xlim(-0.05, 0.85)
    ax.set_ylim(-0.62, 1.62)
    ax.set_yticks([1.0, 0.0], ["Ischemic\ncardiomyopathy", "Dilated\ncardiomyopathy"])
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.axvline(0, color=GUIDE, linestyle=(0, (5, 3)), linewidth=1.0, zorder=0)
    ax.set_xlabel("Frozen-weight vs HF-effect Spearman $\\rho$\n(95% donor-bootstrap CI)",
                  fontsize=10.5, color=TEXT, labelpad=8)
    style_axes(ax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0, labelsize=8.8, pad=8)

    for group, color, y in specs:
        row = estimates[group]
        est = float(row["estimate"])
        lo = float(row["ci_lo"])
        hi = float(row["ci_hi"])
        ax.errorbar(est, y, xerr=np.array([[est - lo], [hi - est]]), fmt="o", color=color,
                    markersize=7.3, elinewidth=1.6, capsize=4.5, capthick=1.5, zorder=3)
        ax.text(0.57, y + 0.045, f"{est:.2f}", color=color, fontsize=10.5, weight="bold", va="bottom")
        ax.text(0.57, y - 0.045, f"({lo:.2f}–{hi:.2f})", color="#45556A", fontsize=8.5, va="top")

    legend_items = [
        Line2D([], [], marker="o", linestyle="none", color=RED, markersize=6.5, label="ICM"),
        Line2D([], [], marker="o", linestyle="none", color=BLUE, markersize=6.5, label="DCM"),
    ]
    ax.legend(handles=legend_items, loc="upper right", frameon=False, fontsize=8,
              handletextpad=0.35, columnspacing=0.8, ncol=2)


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=400, bbox_inches="tight", facecolor=WHITE)
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight", facecolor=WHITE)


def main() -> None:
    matplotlib.rcParams.update({
        "font.family": "Arial",
        "font.size": 9,
        "axes.labelcolor": TEXT,
        "text.color": TEXT,
        "figure.facecolor": WHITE,
        "savefig.facecolor": WHITE,
        "svg.fonttype": "none",
    })
    genes = read_table("adult_lv_gene_effects.csv")
    sign = read_table("sign_concordance.csv")
    # Reproduce the accepted Figure 3 source table; the later BC02 replay uses
    # a different bootstrap seed and therefore has slightly different limits.
    f3c = read_table("Source_F3c_Etiology_Estimates.csv", FIGURE3_SOURCE_DIR)

    fig_a, ax_a = plt.subplots(figsize=(7.5, 7.5))
    draw_panel_a(ax_a, genes)
    fig_a.subplots_adjust(left=0.18, bottom=0.15, right=0.97, top=0.97)
    save_figure(fig_a, "Figure3a_scatter")
    plt.close(fig_a)

    fig_b, ax_b = plt.subplots(figsize=(7.2, 3.65))
    draw_panel_b(ax_b, sign, genes)
    fig_b.subplots_adjust(left=0.33, bottom=0.23, right=0.97, top=0.89)
    save_figure(fig_b, "Figure3b_concordance")
    plt.close(fig_b)

    fig_c, ax_c = plt.subplots(figsize=(7.2, 4.15))
    draw_panel_c(ax_c, f3c)
    fig_c.subplots_adjust(left=0.35, bottom=0.22, right=0.97, top=0.92)
    save_figure(fig_c, "Figure3c_subtypes")
    plt.close(fig_c)

    fig = plt.figure(figsize=(14.48, 10.86), facecolor=WHITE)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.0, 0.94],
                            left=0.055, right=0.985, bottom=0.075, top=0.98,
                            wspace=0.18, hspace=0.19)
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])
    draw_panel_a(ax_a, genes)
    draw_panel_b(ax_b, sign, genes)
    draw_panel_c(ax_c, f3c)
    save_figure(fig, "Figure3_three_panels")
    plt.close(fig)


if __name__ == "__main__":
    main()
