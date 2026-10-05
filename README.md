# Adult heart failure molecular anchoring and functional mapping diverge across cardiac systems

Reproducibility release v1.0.0 for the study by Zisen Zhou, Auckland Bioengineering Institute, University of Auckland.

## Contents

- `figures/main/`: Main Figures 1–5, plotting scripts, derived source tables, final exports, and the editable main-figure deck.
- `figures/supplement/`: Supplementary Figures 1–8, each with PNG, SVG, source CSV tables, and a build script.
- `analysis/validate_figure_data.py`: recomputes headline numerical checks from the released source tables.
- `DATA_SOURCES.md` and `THIRD_PARTY_NOTICES.md`: public accessions, source attribution, and file-specific license boundaries.

## Reproduce

Use Python 3.12 with the packages in `environment.yml`. From the repository root, run `conda env create -f environment.yml`, then `conda activate cardiac-signature-transportability`, then `python reproduce_figures.py`. To check headline statistics and cross-figure consistency, run `python analysis/validate_figure_data.py`.

The scripts rebuild Main Figures 2–5 and Supplementary Figures 1–8 from bundled figure-level source tables. Main Figure 1 is an explanatory diagram; its raster export and editable PowerPoint source are included. The CSVs are author-derived figure source data; raw expression and assay datasets are not redistributed. Their GEO accessions are listed in `DATA_SOURCES.md`. The upstream raw-data processing pipeline is not included.

The MIT license covers original code only. It does not license non-code figures or source tables. The two gene-level tables containing ReHeaT reference weights have a separate CC BY-NC 4.0 notice in `THIRD_PARTY_NOTICES.md`.

Repository: https://github.com/zzisen/human-cardiac-signature-transportability
Zenodo archive (v1.0.0): https://doi.org/10.5281/zenodo.23156479
