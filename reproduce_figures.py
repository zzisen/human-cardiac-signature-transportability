from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
JOBS = [
    ("figures/main/Figure2", "build_all.py"),
    ("figures/main/Figure3", "make_figure3.py"),
    ("figures/main/Figure4", "build_figure4_panels.py"),
    ("figures/main/Figure5", "plot_figure5_final_abc.py"),
    *[(f"figures/supplement/Figure S{i}", "build_figure.py") for i in range(1, 9)],
]

for directory, script in JOBS:
    working_directory = ROOT / directory
    print(f"Building {directory}", flush=True)
    subprocess.run([sys.executable, script], cwd=working_directory, check=True)

print("All scripted figure exports were regenerated.")
