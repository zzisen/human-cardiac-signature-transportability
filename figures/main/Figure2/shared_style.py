"""Shared reference-matched styling and source paths for Figure 2 panels."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from matplotlib.text import Text


HERE = Path(__file__).resolve().parent
FIGURE_DIR = HERE.parent
DATA_DIR = FIGURE_DIR / "data"
if not DATA_DIR.exists():
    # The source bundle includes the accepted CSV snapshots needed by all panels.
    DATA_DIR = HERE / "data"
OUTPUT_DIR = HERE / "Output"

# Colors sampled/visually matched to Figure_Tem/Figure2.png.
PALETTE = {
    "coral": "#F07F70",
    "coral_dark": "#D84D2F",
    "blue": "#14558F",
    "blue_bright": "#2971CA",
    "blue_band": "#A9C7E4",
    "blue_light": "#9BC7F2",
    "teal": "#2C8781",
    "purple": "#7041AE",
    "grey": "#D0D2D4",
    "mid_grey": "#929DA6",
    "dark_grey": "#20262D",
    "axis": "#242A30",
    "neutral": "#8F9BA5",
    "white": "#FFFFFF",
}

FONT = "Arial"
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": [FONT, "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.edgecolor": PALETTE["axis"],
    "axes.labelcolor": PALETTE["dark_grey"],
    "axes.linewidth": 0.8,
    "xtick.color": PALETTE["dark_grey"],
    "ytick.color": PALETTE["dark_grey"],
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "text.color": PALETTE["dark_grey"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "figure.facecolor": PALETTE["white"],
    "savefig.facecolor": PALETTE["white"],
})


def style_axes(ax, *, grid=False):
    """Use the clean axes and tick treatment from the supplied reference."""
    ax.set_facecolor(PALETTE["white"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(PALETTE["axis"])
    ax.spines["bottom"].set_color(PALETTE["axis"])
    ax.tick_params(length=3.2, width=0.8, direction="out")
    if grid:
        ax.grid(True, color="#E6E9ED", linewidth=0.55, zorder=0)
        ax.set_axisbelow(True)
    else:
        ax.grid(False)


def scale_panel_fonts(ax, factor=2.0):
    """Scale all text artists in a panel, including ticks, annotations and legend."""
    for text in ax.findobj(match=Text):
        text.set_fontsize(text.get_fontsize() * factor)


def _tight_width_inches(fig, pad_inches):
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer())
    return bbox.width + 2 * pad_inches


def save_panel(fig, stem, *, target_png_width_px=None):
    """Write both a high-resolution raster and an editable vector source."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    png_path = OUTPUT_DIR / f"{stem}.png"
    dpi = 300
    pad_inches = 0.08
    if target_png_width_px is None:
        fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=pad_inches,
                    facecolor=PALETTE["white"], metadata={"Software": "Matplotlib"})
    else:
        # Tight cropping makes output width differ slightly from figsize. Adjust
        # the canvas until the actual PNG width matches the requested pixel width.
        for _ in range(12):
            fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=pad_inches,
                        facecolor=PALETTE["white"], metadata={"Software": "Matplotlib"})
            with Image.open(png_path) as rendered:
                actual_width = rendered.width
            if actual_width == target_png_width_px:
                break

            width, height = fig.get_size_inches()
            base_extent = _tight_width_inches(fig, pad_inches)
            probe = 0.25
            fig.set_size_inches(width + probe, height, forward=False)
            probe_extent = _tight_width_inches(fig, pad_inches)
            slope = (probe_extent - base_extent) / probe
            fig.set_size_inches(width, height, forward=False)
            if slope <= 0:
                raise RuntimeError("Could not determine tight-cropped panel width scaling.")
            correction = (target_png_width_px - actual_width) / dpi / slope
            fig.set_size_inches(width + correction, height, forward=False)
        else:
            raise RuntimeError(
                f"Could not match requested PNG width {target_png_width_px}px for {stem}."
            )

    fig.savefig(OUTPUT_DIR / f"{stem}.svg", format="svg", bbox_inches="tight",
                pad_inches=pad_inches, facecolor=PALETTE["white"], metadata={"Date": None})
    plt.close(fig)
