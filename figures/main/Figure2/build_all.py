"""Generate all four reference-styled Figure 2 panels from accepted data."""

from panel_b import build_panel_b
from panel_c import build_panel_c
from panel_d import build_panel_d
from panel_a import build_panel_a
from PIL import Image
from shared_style import OUTPUT_DIR


if __name__ == "__main__":
    build_panel_b()
    build_panel_c()
    build_panel_d()
    lower_widths = []
    for panel in "bcd":
        with Image.open(OUTPUT_DIR / f"Figure2_panel_{panel}.png") as rendered:
            lower_widths.append(rendered.width)
    target_width = sum(lower_widths)
    build_panel_a(target_width_px=target_width)
    print(
        "Wrote four Figure 2 real-data panels (PNG and SVG); "
        f"Panel A width matches the sum of Panels B–D ({target_width}px)."
    )
