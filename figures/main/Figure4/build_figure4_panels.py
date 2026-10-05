#!/usr/bin/env python3
"""Rebuild three Figure 4 panels from the project's accepted source data.

Inputs are the Figure 4 source tables plus the accepted stable matched-set
distributions and benchmark summaries. The script only formats existing
results; it does not recompute the biological analyses.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

NAVY = "#13243A"
TEAL = "#078F83"
ORANGE = "#F16F53"
BLUE = "#3766B2"
PALE_BLUE = "#DCE5F7"
AXIS = "#657184"
ZERO = "#AEB7C4"
CONNECTOR = "#A9AEB7"


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    return float(row[key])


def style_axis(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_color(AXIS)
        spine.set_linewidth(0.85)
    ax.tick_params(axis="both", colors=NAVY, labelsize=9, length=4, width=0.8)
    ax.xaxis.label.set_color(NAVY)
    ax.yaxis.label.set_color(NAVY)
    ax.grid(False)


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(HERE / name, dpi=300, facecolor="white")
    fig.savefig(HERE / Path(name).with_suffix(".svg").name, facecolor="white")
    plt.close(fig)


def panel_a_associations() -> None:
    """Panel A: pooled and pan-etiology component associations with Force."""
    rows = read_csv("Source_F4a_Component_Associations.csv")
    definitions = ["Pooled", "Pan-etiology"]
    y_positions = {"Pooled": 1.0, "Pan-etiology": 0.0}
    component_colors = {"Concordant": TEAL, "Discordant": ORANGE}

    fig, ax = plt.subplots(figsize=(12.2, 3.65))
    ax.axvline(0, color=ZERO, lw=1.0, ls=(0, (4, 3)), zorder=0)

    for definition in definitions:
        group = [row for row in rows if row["definition"] == definition]
        y = y_positions[definition]
        concordant = next(row for row in group if row["component"] == "Concordant")
        discordant = next(row for row in group if row["component"] == "Discordant")
        ax.hlines(
            y,
            number(concordant, "rho"),
            number(discordant, "rho"),
            color=CONNECTOR,
            lw=1.2,
            zorder=1,
        )

        for row in (concordant, discordant):
            rho = number(row, "rho")
            lo = number(row, "ci_lo")
            hi = number(row, "ci_hi")
            color = component_colors[row["component"]]
            ax.errorbar(
                rho,
                y,
                xerr=np.array([[rho - lo], [hi - rho]]),
                fmt="o",
                color=color,
                markersize=7.4,
                elinewidth=1.45,
                capsize=3.3,
                capthick=1.25,
                zorder=3,
            )
            minus = "\u2212"
            value_label = f"{rho:.2f}".replace("-", minus)
            lo_label = f"{lo:.2f}".replace("-", minus)
            hi_label = f"{hi:.2f}".replace("-", minus)
            ax.text(
                rho,
                y - 0.17,
                f"{value_label}\n[{lo_label}, {hi_label}]",
                ha="center",
                va="top",
                color=color,
                fontsize=8.5,
                fontweight="semibold",
                linespacing=1.2,
            )

    legend = [
        Line2D(
            [0], [0], marker="o", color=TEAL, lw=1.3, markersize=7,
            label="HF-concordant association",
        ),
        Line2D(
            [0], [0], marker="o", color=ORANGE, lw=1.3, markersize=7,
            label="HF-discordant association",
        ),
    ]
    ax.legend(
        handles=legend,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=True,
        fancybox=True,
        framealpha=0.98,
        edgecolor="#D0D7E1",
        fontsize=9,
        handlelength=1.4,
        columnspacing=1.8,
    )
    ax.set_xlim(-0.97, 0.69)
    ax.set_ylim(-0.55, 1.57)
    ax.set_xticks(np.arange(-0.8, 0.61, 0.2))
    ax.set_yticks([1, 0], ["Pooled", "Pan-etiology"])
    for label in ax.get_yticklabels():
        label.set_fontweight("semibold")
    ax.set_xlabel("Spearman \u03c1 with 24-hour Force", fontsize=10)
    style_axis(ax)
    fig.subplots_adjust(left=0.13, right=0.99, bottom=0.22, top=0.97)
    save(fig, "Figure4A_component_associations.png")


def panel_b_preference() -> None:
    """Panel B: accepted paired-bootstrap component preference estimates."""
    rows = read_csv("Source_F4a_Component_Preferences.csv")
    y_positions = {"Pooled": 1.0, "Pan-etiology": 0.0}

    fig, ax = plt.subplots(figsize=(5.9, 4.0))
    ax.axvline(0, color=ZERO, lw=1.0, ls=(0, (4, 3)), zorder=0)

    for row in rows:
        y = y_positions[row["definition"]]
        estimate = number(row, "preference")
        lo = number(row, "ci_lo")
        hi = number(row, "ci_hi")
        ax.fill_betweenx(
            [y - 0.13, y + 0.13],
            lo,
            hi,
            color=PALE_BLUE,
            alpha=0.9,
            linewidth=0,
            zorder=1,
        )
        ax.hlines(y, lo, hi, color=BLUE, lw=1.5, zorder=2)
        ax.vlines([lo, hi], y - 0.055, y + 0.055, color=BLUE, lw=1.25, zorder=2)
        ax.scatter([estimate], [y], s=52, color=BLUE, zorder=3)

        estimate_label = f"{estimate:.2f}"
        lo_label = f"{lo:.2f}"
        hi_label = f"{hi:.2f}"
        ax.text(
            estimate,
            y + 0.09,
            f"{estimate_label}\n[{lo_label}, {hi_label}]",
            ha="center",
            va="bottom",
            color="#234A91",
            fontsize=8.5,
            fontweight="semibold",
            linespacing=1.2,
        )

    legend = [
        Line2D(
            [0], [0], marker="o", color=BLUE, lw=0, markersize=6.8,
            label="Preference estimate",
        ),
        Patch(
            facecolor=PALE_BLUE,
            edgecolor="none",
            label="95% paired ligand-bootstrap CI",
        ),
    ]
    ax.legend(
        handles=legend,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=2,
        frameon=True,
        fancybox=True,
        framealpha=0.98,
        edgecolor="#D0D7E1",
        fontsize=8.5,
        handlelength=1.4,
        columnspacing=1.7,
    )
    ax.set_xlim(-0.2, 1.0)
    ax.set_ylim(-0.52, 1.6)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks([1, 0], ["Pooled", "Pan-etiology"])
    for label in ax.get_yticklabels():
        label.set_fontweight("semibold")
    ax.set_xlabel(
        "Preference (|concordant \u03c1| \u2212 |discordant \u03c1|)",
        fontsize=10,
    )
    style_axis(ax)
    fig.subplots_adjust(left=0.18, right=0.99, bottom=0.20, top=0.97)
    save(fig, "Figure4B_component_preference.png")


def matched_density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Small Gaussian-kernel display smoother, avoiding extra dependencies."""
    bandwidth = 1.06 * np.std(values, ddof=1) * len(values) ** (-1 / 5)
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError("Matched-set values do not support a density estimate.")
    z = (grid[:, None] - values[None, :]) / bandwidth
    return np.exp(-0.5 * z * z).mean(axis=1) / (
        bandwidth * np.sqrt(2 * np.pi)
    )


def panel_c_matched_nulls() -> None:
    """Panel C: stable HF-concordant matched-set distributions for both definitions."""
    source_rows = {
        row["definition"]: row
        for row in read_csv("Source_F4b_Matched_Gene_Benchmark.csv")
    }
    sources = [
        {
            "definition": "Pooled",
            "distribution": "BC03a_STABLE_MATCHED_SET_DISTRIBUTION.csv",
            "benchmark": "BC03a_STABLE_MATCHED_SET_BENCHMARK.csv",
            "rho_key": "stable_concordant_rho",
        },
        {
            "definition": "Pan-etiology",
            "distribution": "BC05a_MATCHED_SET_DISTRIBUTION.csv",
            "benchmark": "BC05a_MATCHED_SET_BENCHMARK.csv",
            "rho_key": "target_rho",
        },
    ]

    fig, axes = plt.subplots(1, 2, figsize=(6.3, 4.0), sharey=True)
    for i, (ax, spec) in enumerate(zip(axes, sources)):
        distribution_rows = read_csv(spec["distribution"])
        values = np.asarray(
            [float(row["B5_between_rho"]) for row in distribution_rows],
            dtype=float,
        )
        benchmark = next(
            row for row in read_csv(spec["benchmark"])
            if row["system"] == "B5" and row["hierarchy"] == "between"
        )
        observed = number(benchmark, spec["rho_key"])
        empirical_p = number(benchmark, "empirical_p_absrho")
        expected_p = number(source_rows[spec["definition"]], "empirical_p")
        if len(values) != 2000:
            raise ValueError(
                f"{spec['distribution']} contains {len(values)} draws; expected 2000."
            )
        if not np.isclose(empirical_p, expected_p, atol=5e-5):
            raise ValueError(
                f"{spec['definition']} benchmark disagrees with accepted Figure 4 source."
            )

        x_low, x_high = -0.85, 0.05
        grid = np.linspace(x_low, x_high, 901)
        density = matched_density(values, grid)
        ax.fill_between(grid, density, color=TEAL, alpha=0.18, linewidth=0, zorder=1)
        ax.plot(
            grid,
            density,
            color=TEAL,
            lw=1.65,
            zorder=2,
        )
        ax.axvline(0, color=ZERO, lw=1.0, ls=(0, (4, 3)), zorder=0)
        ax.axvline(
            observed,
            color=TEAL,
            lw=2.0,
            label="Observed association",
            zorder=3,
        )

        minus = "\u2212"
        observed_label = f"{observed:.3f}".replace("-", minus)
        p_label = f"{expected_p:.3f}"
        ax.text(
            0.23,
            0.955,
            f"{spec['definition']} stable HF-concordant",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color=NAVY,
            fontsize=8.3,
            fontweight="semibold",
        )
        ax.text(
            0.97,
            0.79,
            f"Observed \u03c1 = {observed_label}\nEmpirical P = {p_label}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            color=NAVY,
            fontsize=8.2,
            linespacing=1.3,
            bbox={
                "boxstyle": "round,pad=0.3",
                "facecolor": "white",
                "edgecolor": "#D7DDE5",
                "alpha": 0.95,
            },
        )
        ax.legend(
            handles=[
                Patch(
                    facecolor=TEAL,
                    edgecolor="none",
                    alpha=0.18,
                    label="Mathcad-set null (n = 2,000)",
                ),
                Line2D(
                    [0], [0], color=TEAL, lw=2.0,
                    label="Observed association",
                ),
            ],
            loc="lower right",
            frameon=True,
            fancybox=True,
            framealpha=0.97,
            edgecolor="#D0D7E1",
            fontsize=7.6,
            handlelength=1.4,
        )
        ax.set_xlim(x_low, x_high)
        ax.set_ylim(0, float(density.max()) * 1.28)
        ax.set_xticks([-0.8, -0.6, -0.4, -0.2, 0])
        if i == 0:
            ax.set_ylabel("Density", fontsize=9.5)
        ax.set_xlabel("Spearman \u03c1 with 24-hour Force", fontsize=9.2)
        style_axis(ax)

    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.19, top=0.98, wspace=0.08)
    save(fig, "Figure4C_matched_gene_nulls.png")


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.labelcolor": NAVY,
            "text.color": NAVY,
            "mathtext.fontset": "dejavusans",
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )
    panel_a_associations()
    panel_b_preference()
    panel_c_matched_nulls()
    print(f"Saved three Figure 4 panels under: {HERE}")


if __name__ == "__main__":
    main()
