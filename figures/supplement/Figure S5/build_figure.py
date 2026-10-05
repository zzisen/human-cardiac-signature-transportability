from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox


SCRIPT_DIR = Path(__file__).resolve().parent
BUNDLE_MODE = bool(re.fullmatch(r"Figure S[1-8]", SCRIPT_DIR.name))
if BUNDLE_MODE:
    SOURCE = SCRIPT_DIR
    OUT = SCRIPT_DIR
else:
    ROOT = SCRIPT_DIR.parents[1]
    SOURCE = ROOT / "06_Figures" / "COMM_BIOL_V3" / "Supplement"
    OUT = ROOT / "06_Figures" / "COMM_BIOL_V3" / "Reference_Style_Redraw"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#0867D5"
VIOLET = "#7B42D8"
GRAY = "#87919B"
TEAL = "#00A6A6"
CORAL = "#FF625E"
NAVY = "#0A2458"
PALE = "#D9E1EB"
WHITE = "#FFFFFF"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8.2,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 8,
    "text.color": NAVY,
    "axes.labelcolor": NAVY,
    "axes.edgecolor": NAVY,
    "xtick.color": NAVY,
    "ytick.color": NAVY,
    "savefig.facecolor": WHITE,
    "figure.facecolor": WHITE,
})


def table(fig, name):
    source_file = SOURCE / name if BUNDLE_MODE else SOURCE / f"Supplementary_Figure_{fig}" / name
    return pd.read_csv(source_file)


def style_axis(ax, zero=None, xlim=None, ylim=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.75)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.tick_params(axis="both", width=0.65, length=3, pad=3)
    ax.grid(False)
    if zero is not None:
        ax.axvline(zero, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)


def section(ax, letter, label):
    ax.text(0, 1.10, f"{letter}.  {label}", transform=ax.transAxes,
            ha="left", va="bottom", color=NAVY, fontsize=9.2,
            fontweight="bold", clip_on=False)


def point_interval(ax, x, lo, hi, y, color, marker="o", ms=5.1,
                   open_marker=False, lw=1.05, cap=2.1, z=3):
    x, lo, hi, y = map(float, (x, lo, hi, y))
    ax.errorbar(x, y, xerr=[[max(0.0, x - lo)], [max(0.0, hi - x)]],
                fmt=marker, color=color, ecolor=color, markersize=ms,
                markerfacecolor=(WHITE if open_marker else color),
                markeredgecolor=color, markeredgewidth=0.8,
                elinewidth=lw, capsize=cap, capthick=lw, zorder=z)


def save(fig, number):
    stem = f"Figure S{number}" if BUNDLE_MODE else f"Supplementary_Figure_{number}"
    if number == 5:
        fig.canvas.draw()
        tight = fig.get_tightbbox(fig.canvas.get_renderer())
        # Preserve the pre-patch canvas bounds and data-region position after
        # the S5 labels were expanded; only whitespace around the plot is adjusted.
        bbox = Bbox.from_extents(
            tight.x0 - (64 / 300), tight.y0 - (26 / 300),
            tight.x1 + (24 / 300), tight.y1 + (21 / 300))
    else:
        bbox = "tight"
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches=bbox, pad_inches=0.08)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches=bbox, pad_inches=0.08)
    plt.close(fig)

def fig1():
    h = table(1, "Source_S1a_Hierarchy.csv")
    sh = table(1, "Source_S1b_Split_Half.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.25, 6.4),
                               gridspec_kw={"height_ratios": [1.15, 0.65], "hspace": 0.62})
    section(a, "a", "Hierarchy effect differs by system")
    a.set_yticks([1, 0], ["Organoid screen\n(24-h Force)", "Cardiomyocyte screen\n(CaT amplitude)"])
    a.set_ylim(-0.55, 1.55)
    a.set_xlim(-0.9, 0.82)
    a.axvline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    for system, yy in [("B5", 1), ("B6", 0)]:
        entries = h[h.system == system]
        for _, r in entries.iterrows():
            between = r.hierarchy == "between"
            point_interval(a, r.rho, r.ci_lo, r.ci_hi, yy + (0.055 if between else -0.055),
                           BLUE if between else GRAY, marker="s" if between else "o", ms=5.0)
    a.set_xlabel("Spearman ρ (total score)")
    a.legend(handles=[
        Line2D([0], [0], marker="s", color=BLUE, lw=0, markersize=5.2, label="Between perturbations"),
        Line2D([0], [0], marker="o", color=GRAY, lw=0, markersize=5.2, label="Within perturbations"),
    ], loc="upper right", frameon=True, framealpha=1, facecolor=WHITE,
       edgecolor=PALE, fontsize=7.4, borderpad=0.45, handletextpad=0.5)
    style_axis(a)
    section(b, "b", "Split-half component repeatability")
    b.set_yticks([1, 0], ["Pooled", "Pan-etiology"])
    b.set_xlim(0.65, 0.92)
    b.set_ylim(-0.55, 1.55)
    for y, definition in [(1, "Pooled"), (0, "Pan")]:
        conc = sh[sh.component.str.startswith(definition) & sh.component.str.contains("concordant")].iloc[0]
        disc = sh[sh.component.str.startswith(definition) & sh.component.str.contains("discordant")].iloc[0]
        x_conc = float(conc.median_split_half_rho)
        x_disc = float(disc.median_split_half_rho)
        b.plot([x_disc, x_conc], [y, y], color=GRAY, lw=1.35, zorder=1)
        b.plot(x_conc, y, marker="o", color=TEAL, markersize=5.8, lw=0, zorder=3)
        b.plot(x_disc, y, marker="s", color=CORAL, markersize=5.5, lw=0, zorder=3)
    b.set_xlabel("Median split-half correlation")
    b.legend(handles=[
        Line2D([0], [0], marker="o", color=TEAL, lw=0, markersize=5.5, label="HF-concordant"),
        Line2D([0], [0], marker="s", color=CORAL, lw=0, markersize=5.3, label="HF-discordant"),
    ], loc="lower right", frameon=True, framealpha=1, facecolor=WHITE,
       edgecolor=PALE, fontsize=7.2, borderpad=0.4, handletextpad=0.45)
    style_axis(b)
    fig.subplots_adjust(left=0.29, right=0.98, top=0.94, bottom=0.10)
    save(fig, 1)


def fig2():
    cf = table(2, "Source_S2a_Class_Conditioning.csv")
    adj = table(2, "Source_S2b_Adjustment_Summaries.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.35, 8.0),
                               gridspec_kw={"height_ratios": [0.95, 1.15], "hspace": 0.66})
    section(a, "a", "Class-conditioned and leave-one-class-out preference")
    a.set_xlim(0, 1.2)
    a.set_ylim(-0.45, 1.55)
    a.axvline(0, color=GRAY, lw=0.75, ls=(0, (3, 3)), zorder=0)
    a.set_yticks([1, 0], ["LOO range across 7 classes", "Class-conditioned"])
    r = adj[(adj.analysis == "Leave-one-class-out range") & (adj.definition == "Pooled")].iloc[0]
    a.hlines(1, r.estimate_or_min, r.lo_or_max, color=BLUE, lw=2.2)
    a.vlines([r.estimate_or_min, r.lo_or_max], 0.88, 1.12, color=BLUE, lw=1.35)
    for y, definition, offset in [(0.085, "Pooled", 0.085), (-0.085, "Pan-etiology", -0.085)]:
        r = cf[cf.definition == definition].iloc[0]
        point_interval(a, r.rho_or_preference, r.ci_lo, r.ci_hi, offset,
                       BLUE if definition == "Pooled" else VIOLET,
                       marker="o" if definition == "Pooled" else "D", ms=5.0)
    a.set_xlabel("Preference (|ρ concordant| − |ρ discordant|)")
    for ax in (a,):
        ax.legend(handles=[
            Line2D([0], [0], marker="o", color=BLUE, lw=0, markersize=5.2, label="Pooled"),
            Line2D([0], [0], marker="D", color=VIOLET, lw=0, markersize=5.0, label="Pan-etiology"),
        ], loc="upper right", frameon=True, framealpha=1, facecolor=WHITE,
           edgecolor=PALE, fontsize=7.3, borderpad=0.4, handletextpad=0.45)
    style_axis(a)

    section(b, "b", "Separate same-support technical adjustments")
    b.set_xlim(0, 1.2)
    b.set_ylim(-0.45, 1.55)
    b.axvline(0, color=GRAY, lw=0.75, ls=(0, (3, 3)), zorder=0)
    b.set_yticks([1, 0], ["Plate-adjusted", "Generic-response-burden adjusted"])
    for y, analysis in [(1, "Plate-adjusted"), (0, "Response-burden adjusted")]:
        for definition, offset in [("Pooled", 0.085), ("Pan-etiology", -0.085)]:
            r = adj[(adj.analysis == analysis) & (adj.definition == definition)].iloc[0]
            point_interval(b, r.estimate_or_min, r.lo_or_max, r.hi, y + offset,
                           BLUE if definition == "Pooled" else VIOLET,
                           marker="o" if definition == "Pooled" else "D", ms=4.9)
    b.set_xlabel("Preference (|ρ concordant| − |ρ discordant|)")
    b.legend(handles=[
        Line2D([0], [0], marker="o", color=BLUE, lw=0, markersize=5.2, label="Pooled"),
        Line2D([0], [0], marker="D", color=VIOLET, lw=0, markersize=5.0, label="Pan-etiology"),
    ], loc="upper right", frameon=True, framealpha=1, facecolor=WHITE,
       edgecolor=PALE, fontsize=7.3, borderpad=0.4, handletextpad=0.45)
    style_axis(b)
    fig.subplots_adjust(left=0.39, right=0.98, top=0.94, bottom=0.09)
    save(fig, 2)


def fig3():
    d = table(3, "Source_S3_CaT_Endpoint_Support.csv")
    fig, ax = plt.subplots(figsize=(6.25, 5.4))
    order = ["Primary frozen score", "Sub-cytotoxic support", "No observed MEC",
             "≤25× Cmax", "≤10× Cmax"]
    labels = ["Primary score", "Strict sub-cytotoxic support", "No observed MEC", "≤25× Cmax", "≤10× Cmax"]
    ax.set_yticks(np.arange(len(order) - 1, -1, -1), labels)
    ax.set_ylim(-0.65, 4.65)
    ax.set_xlim(-0.4, 0.6)
    ax.axvline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    for i, name in enumerate(order):
        r = d[d.analysis == name].iloc[0]
        primary = i == 0
        point_interval(ax, r.rho, r.ci_lo, r.ci_hi, len(order) - 1 - i,
                       BLUE if primary else GRAY, marker="D" if primary else "o", ms=5.4 if primary else 4.8)
    ax.set_xlabel("Within-compound Spearman ρ")
    style_axis(ax)
    fig.subplots_adjust(left=0.35, right=0.97, top=0.96, bottom=0.14)
    save(fig, 3)


def fig4():
    d = table(4, "Source_S4_Additive_Architecture.csv")
    fig, ax = plt.subplots(figsize=(6.35, 5.7))
    terms = ["Positive-weight variance", "2 × signed covariance", "Negative-weight variance"]
    ax.set_yticks([2, 1, 0], terms)
    ax.set_ylim(-0.65, 2.65)
    ax.set_xlim(-16, 11)
    ax.axvline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    for y, term in zip([2, 1, 0], ["Positive-weight variance", "Twice signed covariance", "Negative-weight variance"]):
        vals = d[d.term == term]
        for _, r in vals.iterrows():
            org = r.system == "Organoid"
            point_interval(ax, r.scaled_by_total_variance, r.ci_lo, r.ci_hi,
                           y + (0.085 if org else -0.085), BLUE if org else VIOLET,
                           marker="o" if org else "D", ms=5.0)
    ax.set_xlabel("Variance term / total-score variance")
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color=BLUE, lw=0, markersize=5.3, label="Organoid dataset"),
        Line2D([0], [0], marker="D", color=VIOLET, lw=0, markersize=5.0, label="Cardiomyocyte dataset"),
    ], loc="upper right", frameon=True, framealpha=1, facecolor=WHITE,
       edgecolor=PALE, fontsize=7.3, borderpad=0.45, handletextpad=0.5)
    style_axis(ax)
    fig.subplots_adjust(left=0.35, right=0.98, top=0.91, bottom=0.15)
    save(fig, 4)


def fig5():
    axes_data = table(5, "Source_S5a_Context_Axis.csv")
    culture = table(5, "Source_S5b_Culture_Contrasts.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.45, 7.75),
                               gridspec_kw={"height_ratios": [1.08, 0.9], "hspace": 0.68})
    section(a, "a", "Organoid vs cardiomyocyte axis associations")
    ax_rows = axes_data.iloc[::-1].reset_index(drop=True)
    label_wrap = {
        "Sarcomere contractile": "Sarcomere\ncontractile",
        "Oxidative mitochondrial": "Oxidative\nmitochondrial",
        "Calcium excitation coupling": "Calcium excitation\ncoupling",
        "Non-cardiomyocyte marker context": "Non-cardiomyocyte\nmarker context",
    }
    display_labels = [label_wrap.get(label, label) for label in ax_rows.axis_label]
    a.set_yticks(np.arange(len(ax_rows)), display_labels)
    a.set_ylim(-0.6, 4.6)
    a.set_xlim(-0.55, 0.95)
    a.axvline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    for y, r in enumerate(ax_rows.itertuples()):
        point_interval(a, r.cardiomyocyte_minus_organoid, r.ci_lo, r.ci_hi,
                       y, GRAY, marker="s", ms=4.9)
    a.set_xlabel("Δ Spearman ρ (cardiomyocyte − organoid)")
    style_axis(a)

    section(b, "b", "3D minus 2D culture contrasts")
    groups = ["Sarcomere\ncontractile", "Oxidative\nmitochondrial"]
    b.set_yticks([1, 0], groups)
    b.set_ylim(-0.6, 1.6)
    b.set_xlim(-0.12, 0.015)
    b.axvline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    for y, axis_name in [(1, "Sarcomere contractile"), (0, "Oxidative mitochondrial")]:
        rows = culture[culture.axis_label == axis_name]
        for _, r in rows.iterrows():
            wt = r.contrast.startswith("WT")
            point_interval(b, r.estimate, r.ci_lo, r.ci_hi,
                           y + (0.09 if wt else -0.09), BLUE if wt else VIOLET,
                           marker="o" if wt else "D", ms=5.0)
    b.set_xlabel("Δ Spearman ρ (3D − 2D)")
    b.legend(handles=[
        Line2D([0], [0], marker="o", color=BLUE, lw=0, markersize=5.2, label="Wild-type EHT"),
        Line2D([0], [0], marker="D", color=VIOLET, lw=0, markersize=4.9, label="Dystrophin-deficient EHT"),
    ], loc="lower right", frameon=True, framealpha=1, facecolor=WHITE,
       edgecolor=PALE, fontsize=7.2, borderpad=0.4, handletextpad=0.45)
    style_axis(b)
    fig.subplots_adjust(left=0.34, right=0.98, top=0.94, bottom=0.09)
    save(fig, 5)

def fig6():
    genes = table(6, "Source_S6a_Adult_LV_Gene_Effects.csv")
    concord = table(6, "Source_S6b_Sign_Concordance.csv")
    etiology = table(6, "Source_S6c_Etiology_Estimates.csv")
    fig = plt.figure(figsize=(8.2, 5.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 0.95],
                          height_ratios=[0.9, 1.1], wspace=0.65, hspace=0.78)
    a = fig.add_subplot(gs[:, 0])
    b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 1])

    section(a, "a", "Gene-level association")
    a.scatter(genes.frozen_mol_score_weight, genes.HF_minus_NF_mean_log2_FPKM,
              s=4, c=BLUE, alpha=0.28, linewidths=0)
    a.set_xlabel("Frozen projection weight (a.u.)")
    a.set_ylabel("Adult-LV HF expression effect")
    a.text(0.035, 0.965,
           "Spearman ρ = 0.2385\n(95% CI 0.1508–0.2674)\nn = 12,929 genes",
           transform=a.transAxes, ha="left", va="top", fontsize=8.1,
           color=NAVY, linespacing=1.3)
    style_axis(a)

    section(b, "b", "Sign concordance")
    names = ["Sign", "Absolute-weighted sign"]
    b.set_yticks([1, 0], names)
    b.set_ylim(-0.45, 1.55)
    b.set_xlim(0.45, 0.70)
    for y, r in zip([1, 0], concord.itertuples()):
        color = BLUE if y == 1 else VIOLET
        point_interval(b, r.estimate, r.ci_lo, r.ci_hi, y, color, marker="o" if y == 1 else "D", ms=5.0)
        b.text(r.estimate, y + 0.16, f"{r.estimate*100:.1f}%", color=color,
               va="bottom", ha="center", fontsize=8, fontweight="bold", clip_on=False)
    b.set_xlabel("Concordance")
    b.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x:.0%}"))
    style_axis(b)

    section(c, "c", "Etiology-specific associations")
    order = ["Ischemic cardiomyopathy", "Dilated cardiomyopathy", "ICM minus DCM"]
    labels = ["Ischemic CM", "Dilated CM", "ICM − DCM"]
    c.set_yticks([2, 1, 0], labels)
    c.set_ylim(-0.6, 2.6)
    c.set_xlim(-0.1, 0.36)
    c.axvline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
    for y, key in zip([2, 1, 0], order):
        r = etiology[etiology.contrast == key].iloc[0]
        clr = BLUE if y == 2 else (VIOLET if y == 1 else GRAY)
        point_interval(c, r.estimate, r.ci_lo, r.ci_hi, y, clr,
                       marker="o" if y == 2 else ("D" if y == 1 else "s"), ms=4.9)
    c.set_xlabel("Spearman ρ")
    style_axis(c)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.92, bottom=0.12)
    save(fig, 6)


def fig7():
    d = table(7, "Source_S7_Expression_Technical.csv")
    technical = ["Observed cardiomyocyte distance", "Count-depth matched simulation",
                 "Detection matched simulation"]
    high_detection = ["High-detection organoid reference", "High-detection cardiomyocyte"]
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.35, 5.85),
                               gridspec_kw={"height_ratios": [1.05, 0.78], "hspace": 0.72})
    section(a, "a", "Technical matching sensitivity")
    a.set_yticks([2, 1, 0], ["Observed cardiomyocyte", "Count-depth matched simulation",
                            "Detection matched simulation"])
    a.set_ylim(-0.55, 2.55)
    a.set_xlim(0, 0.72)
    for i, label in enumerate(technical):
        r = d[d.comparison == label].iloc[0]
        y = 2 - i
        observed = i == 0
        point_interval(a, r.distance, r.lo, r.hi, y,
                       BLUE if observed else GRAY,
                       marker="D" if observed else "o", ms=5.1 if observed else 4.6)
        a.text(r.distance + 0.014, y, f"{r.distance:.4f}", ha="left", va="center",
               color=BLUE if observed else NAVY, fontsize=7.3, clip_on=False)
    a.set_xlabel("Fixed expression-space distance")
    style_axis(a)

    section(b, "b", "High-detection subset")
    b.set_yticks([1, 0], ["Organoid reference", "Cardiomyocyte"])
    b.set_ylim(-0.5, 1.5)
    b.set_xlim(0, 0.72)
    for i, label in enumerate(high_detection):
        r = d[d.comparison == label].iloc[0]
        y = 1 - i
        point_interval(b, r.distance, r.lo, r.hi, y, GRAY, marker="o", ms=4.8)
        b.text(r.distance + 0.014, y, f"{r.distance:.4f}", ha="left", va="center",
               color=NAVY, fontsize=7.3, clip_on=False)
    b.set_xlabel("Fixed expression-space distance")
    style_axis(b)
    fig.subplots_adjust(left=0.43, right=0.94, top=0.93, bottom=0.09)
    save(fig, 7)


def fig8():
    d = table(8, "Source_S8_Temporal_Associations.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.35, 7.3),
                               gridspec_kw={"height_ratios": [1, 1], "hspace": 0.70})
    endpoints = ["4 h Force", "24 h Force", "24 h − 4 h Force"]
    x = np.arange(3)
    for ax, panel, prefix in [(a, "a", "Pooled"), (b, "b", "Pan")]:
        section(ax, panel, f"{prefix} component associations")
        ax.set_xticks(x, ["4 h", "24 h", "Change\n(24 h − 4 h)"])
        ax.set_xlim(-0.35, 2.35)
        ax.set_ylim(-0.9, 0.68)
        ax.axhline(0, color=GRAY, lw=0.85, ls=(0, (3, 3)), zorder=0)
        for comp, color, marker, offset, legend_name in [
            ("concordant", BLUE, "o", -0.065, "HF-concordant"),
            ("discordant", VIOLET, "D", 0.065, "HF-discordant"),
        ]:
            rs = []
            for endpoint in endpoints:
                hit = d[(d.endpoint == endpoint) &
                        (d.component.str.lower().str.startswith(prefix.lower())) &
                        (d.component.str.lower().str.endswith(comp))]
                rs.append(hit.iloc[0])
            xx = x + offset
            yy = np.array([r.rho for r in rs], dtype=float)
            ax.plot(xx, yy, color=color, lw=1.15, zorder=2)
            for xi, r in zip(xx, rs):
                ax.errorbar(xi, r.rho, yerr=[[r.rho - r.ci_lo], [r.ci_hi - r.rho]],
                            fmt=marker, color=color, ecolor=color, markersize=4.9,
                            markerfacecolor=color, markeredgecolor=color,
                            markeredgewidth=0.7, elinewidth=0.9, capsize=2.0,
                            capthick=0.9, zorder=3)
        ax.set_ylabel("Component association (Spearman ρ)")
        style_axis(ax)
    b.legend(handles=[
        Line2D([0], [0], marker="o", color=BLUE, lw=1.0, markersize=5.0, label="HF-concordant"),
        Line2D([0], [0], marker="D", color=VIOLET, lw=1.0, markersize=4.8, label="HF-discordant"),
    ], loc="lower right", frameon=True, framealpha=1, facecolor=WHITE,
       edgecolor=PALE, fontsize=7.3, borderpad=0.4, handletextpad=0.5)
    fig.subplots_adjust(left=0.19, right=0.98, top=0.94, bottom=0.09)
    save(fig, 8)


if __name__ == "__main__":
    builders = {1: fig1, 2: fig2, 3: fig3, 4: fig4,
                5: fig5, 6: fig6, 7: fig7, 8: fig8}
    if BUNDLE_MODE:
        match = re.fullmatch(r"Figure S([1-8])", SCRIPT_DIR.name)
        number = int(sys.argv[1]) if len(sys.argv) > 1 else int(match.group(1))
        builders[number]()
        print(f"Wrote Figure S{number} PNG and SVG to {OUT}")
    else:
        for make in builders.values():
            make()
        print(f"Wrote eight figures to {OUT}")
