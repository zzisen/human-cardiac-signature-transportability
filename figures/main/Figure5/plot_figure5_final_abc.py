"""Rebuild the final Figure 5 panels a-c from accepted source data.

The attached Figure_Tem/Figure5.png is used only as the visual reference.
Values are read from the COMM_BIOL_V4 Figure 5 data directory and are never
transcribed from the reference image.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"

COLORS = {
    "Unadjusted": "#104F9C",
    "Class-conditioned": "#2C8BE0",
    "Plate-adjusted": "#6543B5",
    "Generic-response-burden adjusted": "#008C89",
    "Joint burden + class": "#E8533D",
    "concordant": "#2C86D5",
    "discordant": "#F05B43",
    "pooled": "#104F9C",
    "pan": "#F05B43",
    "text": "#142846",
    "axis": "#61758E",
    "guide": "#9BAEC3",
    "white": "#FFFFFF",
}

CROSS_HALF_MODELS = [
    ("Unadjusted", "unadjusted_cross_half_partitions.csv"),
    ("Class-conditioned", "class_adjusted_cross_half_partitions.csv"),
    ("Plate-adjusted", "plate_adjusted_cross_half_partitions.csv"),
    ("Generic-response-burden adjusted", "burden_adjusted_cross_half_partitions.csv"),
    ("Joint burden + class", "joint_burden_class_cross_half_partitions.csv"),
]
PANEL_A_WIDTH_IN = 18.83
PANEL_A_HEIGHT_IN = 7.2
PANEL_B_WIDTH_IN = 8.7
PANEL_C_WIDTH_IN = 8.7
PANEL_HEIGHT_IN = 5.8
FONT_SCALE = 2.0


def read_csv(name: str) -> pd.DataFrame:
    path = DATA_DIR / name
    if not path.is_file():
        raise FileNotFoundError(f"Required Figure 5 source data not found: {path}")
    return pd.read_csv(path)


def check_sources() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    joint_components = read_csv("joint_component_effects.csv")
    temporal = read_csv("temporal_preferences.csv")
    cross_half: dict[str, pd.DataFrame] = {}

    if len(joint_components) != 4:
        raise ValueError("joint_component_effects.csv should contain four accepted component summaries")
    if len(temporal) != 6:
        raise ValueError("temporal_preferences.csv should contain two definitions across three endpoints")

    for label, filename in CROSS_HALF_MODELS:
        frame = read_csv(filename)
        if len(frame) != 1000 or not np.array_equal(frame["partition"].to_numpy(), np.arange(1000)):
            raise ValueError(f"{filename} must preserve all 1,000 ordered cross-half partitions")
        cross_half[label] = frame

    return joint_components, temporal, cross_half


def style_axis(ax: plt.Axes, *, zero_line: bool = True) -> None:
    ax.set_facecolor(COLORS["white"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(COLORS["axis"])
        ax.spines[side].set_linewidth(1.15)
    ax.tick_params(axis="both", colors=COLORS["text"], labelsize=11.5 * FONT_SCALE,
                   direction="out", length=4.5, width=1.0, pad=6)
    ax.set_axisbelow(True)
    if zero_line:
        ax.axhline(0, color=COLORS["guide"], linewidth=1.2,
                   linestyle=(0, (4, 3)), zorder=0)


def set_preference_x(ax: plt.Axes, xmax: float = 0.8) -> None:
    ax.set_xlim(-0.1, xmax)
    ticks = np.arange(0, xmax + 1e-9, 0.2)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{x:.1f}" for x in ticks])
    ax.set_xlabel(
        "Association-magnitude preference\n"
        r"$|\rho_{\mathrm{concordant}}|-|\rho_{\mathrm{discordant}}|$",
        fontsize=14 * FONT_SCALE, color=COLORS["text"], labelpad=10,
    )
    ax.axvline(0, color=COLORS["guide"], linewidth=1.8,
               linestyle=(0, (4, 3)), zorder=0)


def write_figure(fig: plt.Figure, stem: str) -> None:
    for suffix in ("png", "svg"):
        fig.savefig(HERE / f"{stem}.{suffix}", dpi=300, facecolor="white",
                    bbox_inches="tight", pad_inches=0.13)
    plt.close(fig)


def marker_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], marker="o", linestyle="none", markersize=7,
               markerfacecolor=COLORS["pooled"], markeredgecolor=COLORS["pooled"],
               label="Pooled"),
        Line2D([0], [0], marker="D", linestyle="none", markersize=6.5,
               markerfacecolor="white", markeredgecolor=COLORS["text"],
               markeredgewidth=1.2, label="Pan-etiology"),
    ]


def panel_a(cross_half: dict[str, pd.DataFrame]) -> None:
    fig, ax = plt.subplots(figsize=(PANEL_A_WIDTH_IN, PANEL_A_HEIGHT_IN))
    style_axis(ax, zero_line=False)
    set_preference_x(ax, xmax=1.0)
    offsets = {"pooled_preference": -0.14, "pan_preference": 0.14}
    markers = {"pooled_preference": "o", "pan_preference": "D"}
    definitions = (("pooled_preference", "Pooled"),
                   ("pan_preference", "Pan-etiology"))

    for i, (model, _) in enumerate(CROSS_HALF_MODELS):
        frame = cross_half[model]
        color = COLORS[model]
        for column, _label in definitions:
            values = frame[column].to_numpy(dtype=float)
            lo, median, hi = np.quantile(values, [0.025, 0.5, 0.975])
            y = i + offsets[column]
            is_pooled = column == "pooled_preference"
            ax.errorbar(median, y, xerr=[[median - lo], [hi - median]],
                        fmt=markers[column], color=color, ecolor=color,
                        markerfacecolor=color if is_pooled else "white",
                        markeredgecolor=color, markeredgewidth=1.4,
                        markersize=7.2, capsize=4.0, elinewidth=1.8,
                        linewidth=0, zorder=3)
            label_offset = (0, 20) if is_pooled else (0, -20)
            ax.annotate(f"{median:.2f}", (median, y), xytext=label_offset,
                        textcoords="offset points", ha="center", va="center",
                        fontsize=10.4 * FONT_SCALE, color=color, fontweight="semibold")

    labels = [model for model, _ in CROSS_HALF_MODELS]
    ax.set_ylim(-0.55, 4.55)
    ax.set_yticks(range(5), labels)
    ax.invert_yaxis()
    for label, (model, _) in zip(ax.get_yticklabels(), CROSS_HALF_MODELS):
        label.set_color(COLORS[model])
        label.set_fontweight("semibold")
        label.set_fontsize(10.2 * FONT_SCALE)
    ax.set_ylabel("Cross-half analysis", fontsize=14 * FONT_SCALE, color=COLORS["text"], labelpad=10)
    ax.legend(handles=marker_handles(), loc="upper left", ncol=2,
              frameon=True, facecolor="white", edgecolor="#C8D2DF",
              framealpha=0.98, fontsize=9.8 * FONT_SCALE, borderpad=0.5,
              handletextpad=0.45, columnspacing=0.9)
    write_figure(fig, "Figure5a_cross_half_distributions")


def panel_b(joint_components: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(PANEL_B_WIDTH_IN, PANEL_HEIGHT_IN))
    style_axis(ax, zero_line=False)
    ax.set_xlim(-1, 0.5)
    ax.set_xticks([-1, -0.5, 0, 0.5])
    ax.set_xticklabels(["−1.0", "−0.5", "0.0", "0.5"])
    ax.set_xlabel("Component association (ρ)", fontsize=14 * FONT_SCALE, color=COLORS["text"], labelpad=10)
    ax.axhspan(-0.5, 1.5, color=COLORS["concordant"], alpha=0.075, zorder=0)
    ax.axhspan(2.5, 4.5, color=COLORS["discordant"], alpha=0.075, zorder=0)
    ax.axvline(0, color=COLORS["guide"], linewidth=1.8,
               linestyle=(0, (4, 3)), zorder=0)
    ax.axhline(2, color="#D6DDE6", linewidth=1.0, zorder=1)

    rows_to_plot = [
        ("HF-concordant", "concordant", "pooled", 0),
        ("HF-concordant", "concordant", "pan", 1),
        ("HF-discordant", "discordant", "pooled", 3),
        ("HF-discordant", "discordant", "pan", 4),
    ]
    band_backgrounds = {"concordant": "#EFF6FC", "discordant": "#FEF3F1"}
    definition_markers = {"pooled": "o", "pan": "D"}
    for component, color_key, definition_key, y in rows_to_plot:
        color = COLORS[color_key]
        row = joint_components[(joint_components["component"] == component) &
                               (joint_components["definition"].str.lower().str.startswith(definition_key))].iloc[0]
        estimate = float(row["median"])
        lo = float(row["partition_q025"])
        hi = float(row["partition_q975"])
        is_pooled = definition_key == "pooled"
        ax.errorbar(estimate, y, xerr=[[estimate - lo], [hi - estimate]],
                    fmt=definition_markers[definition_key], color=color, ecolor=color,
                    markerfacecolor=color if is_pooled else "white",
                    markeredgecolor=color, markeredgewidth=1.4,
                    markersize=7.2, capsize=4.0, elinewidth=1.8,
                    linewidth=0, zorder=3)
        value_label = f"−\u2009{abs(estimate):.2f}" if estimate < 0 else f"{estimate:.2f}"
        label_offset = (0, 20) if is_pooled else (0, -20)
        ax.annotate(value_label, (estimate, y), xytext=label_offset,
                    textcoords="offset points", ha="center", va="center",
                    fontsize=9.6 * FONT_SCALE, color=color, fontweight="semibold",
                    bbox=dict(facecolor=band_backgrounds[color_key],
                              edgecolor="none", pad=0.35))

    ax.set_ylim(-0.5, 4.5)
    ax.set_yticks([])
    ax.invert_yaxis()
    group_label_transform = ax.get_yaxis_transform()
    ax.text(0.025, 0.5, "HF-\nconcordant", transform=group_label_transform,
            ha="left", va="center",
            color=COLORS["concordant"], fontsize=11.5 * FONT_SCALE, fontweight="semibold")
    ax.text(0.025, 3.5, "HF-\ndiscordant", transform=group_label_transform,
            ha="left", va="center",
            color=COLORS["discordant"], fontsize=11.5 * FONT_SCALE, fontweight="semibold")
    ax.legend(handles=marker_handles(), loc="lower center", bbox_to_anchor=(0.5, 1.03), ncol=2,
              frameon=True, facecolor="white", edgecolor="#C8D2DF",
              framealpha=0.98, fontsize=9.8 * FONT_SCALE, borderpad=0.48,
              handletextpad=0.45, columnspacing=1.0)
    write_figure(fig, "Figure5b_joint_component_effects")


def panel_c(temporal: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(PANEL_C_WIDTH_IN, PANEL_HEIGHT_IN))
    style_axis(ax, zero_line=False)
    endpoints = ["4 h Force", "24 h Force", "24 h − 4 h Force"]
    xbase = np.arange(len(endpoints), dtype=float)
    ax.set_xlim(-0.5, 2.5)
    ax.set_xticks(xbase, ["4 h Force", "24 h Force", "24 h − 4 h\nForce change"])
    ax.set_xlabel("Temporal Force endpoint", fontsize=14 * FONT_SCALE,
                  color=COLORS["text"], labelpad=10)
    ax.set_ylim(-0.1, 1.2)
    preference_ticks = np.arange(0, 1.2001, 0.2)
    ax.set_yticks(preference_ticks)
    ax.set_yticklabels([f"{value:.1f}" for value in preference_ticks])
    ax.axhline(0, color=COLORS["guide"], linewidth=1.8,
               linestyle=(0, (4, 3)), zorder=0)
    series = [
        ("Pooled definition", "Pooled", COLORS["pooled"], -0.14, "o"),
        ("Pan-etiology definition", "Pan-etiology", COLORS["pan"], 0.14, "o"),
    ]
    legend_handles = []

    for definition, label, color, offset, marker in series:
        rows = temporal[temporal["definition"] == definition].set_index("endpoint").loc[endpoints]
        x = xbase + offset
        estimates = rows["preference"].to_numpy(dtype=float)
        lo = rows["ci_lo"].to_numpy(dtype=float)
        hi = rows["ci_hi"].to_numpy(dtype=float)
        ax.plot(x, estimates, color=color, linewidth=2.0, zorder=2)
        ax.errorbar(x, estimates, yerr=[estimates - lo, hi - estimates], fmt=marker,
                    color=color, ecolor=color, markerfacecolor=color,
                    markeredgecolor="white", markeredgewidth=0.9, markersize=7.5,
                    capsize=4.0, elinewidth=1.8, linewidth=0, zorder=3)
        legend_handles.append(Line2D(
            [0], [0], color=color, linewidth=2.0, marker=marker,
            markerfacecolor=color, markeredgecolor="white", markersize=10,
            label=label,
        ))
        for xx, yy, lower, upper in zip(x, estimates, lo, hi):
            if label == "Pooled":
                label_anchor, label_offset = upper, (0, 20)
            else:
                label_anchor, label_offset = lower, (0, -20)
            ax.annotate(f"{yy:.2f}", (xx, label_anchor), xytext=label_offset,
                        textcoords="offset points", ha="center", va="center",
                        fontsize=9.4 * FONT_SCALE, color=color, fontweight="semibold")

    ax.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, 1.03), ncol=2,
              frameon=True, facecolor="white",
              edgecolor="#C8D2DF", framealpha=0.98, fontsize=10 * FONT_SCALE,
              borderpad=0.5, handletextpad=0.45, columnspacing=0.9)
    write_figure(fig, "Figure5c_temporal_preferences")


def main() -> None:
    if not DATA_DIR.is_dir():
        raise FileNotFoundError(f"Figure 5 source directory not found: {DATA_DIR}")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 11.5 * FONT_SCALE,
        "axes.labelsize": 14 * FONT_SCALE,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })
    joint_components, temporal, cross_half = check_sources()
    panel_a(cross_half)
    panel_b(joint_components)
    panel_c(temporal)
    print(f"Wrote final Figure 5 panels a-c to {HERE}")


if __name__ == "__main__":
    main()
