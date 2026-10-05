# Final Figure 5 panels a–c source data

The CSV files required to reproduce the final panels are included in `data/`.
They were copied from the accepted Figure 5 source data and retain their
original values and partition order.

- **Panel a — Cross-half partition distributions:**
  `unadjusted_cross_half_partitions.csv`,
  `class_adjusted_cross_half_partitions.csv`,
  `plate_adjusted_cross_half_partitions.csv`,
  `burden_adjusted_cross_half_partitions.csv`, and
  `joint_burden_class_cross_half_partitions.csv`.
- **Panel b — Joint component associations:** `joint_component_effects.csv`.
- **Panel c — Temporal Force preferences:** `temporal_preferences.csv`.

The self-contained plotting source is `plot_figure5_final_abc.py`. It reads
from this directory's `data/` folder and writes the PNG and SVG panels beside
the script.
