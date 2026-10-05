"""Recompute headline numerical checks from the released figure-source tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]


def close(actual: float, expected: float, *, atol: float = 5e-5) -> None:
    if not np.isclose(actual, expected, atol=atol, rtol=0):
        raise AssertionError(f"{actual:.10g} != {expected:.10g} (atol={atol})")


organoids = pd.read_csv(ROOT / "figures/main/Figure2/data/organoid_ligand_points.csv")
rho_organoid = float(spearmanr(organoids["mean_frozen_score"], organoids["median_force_24h"]).statistic)
assert len(organoids) == 82
assert int(organoids["n_paired_organoids"].sum()) == 382
close(rho_organoid, -0.7050532754323527)

pharm = pd.read_csv(ROOT / "figures/main/Figure2/data/pharmacology_primary_summary.csv")
primary = pharm.loc[pharm["analysis"].eq("Primary frozen score")].iloc[0]
assert int(primary["compounds"]) == 42 and int(primary["dose_groups"]) == 310
close(float(primary["rho"]), 0.035812)
close(float(primary["ci_lo"]), -0.100774, atol=1e-6)
close(float(primary["ci_hi"]), 0.176197, atol=1e-6)

genes = pd.read_csv(ROOT / "figures/main/Figure3/source_data/adult_lv_gene_effects.csv")
rho_adult_lv = float(spearmanr(genes["frozen_projection_weight"], genes["adult_hf_expression_effect"]).statistic)
eligible = genes["nonzero_sign_eligible"].astype(bool)
concordant = genes["sign_concordant"].astype(bool)
sign_concordance = float(concordant[eligible].mean())
weights = genes.loc[eligible, "frozen_projection_weight"].abs().to_numpy(float)
weighted_concordance = float(np.average(concordant[eligible].to_numpy(float), weights=weights))
assert len(genes) == 12_929
close(rho_adult_lv, 0.23846561112775183)
close(sign_concordance, 0.590920)
close(weighted_concordance, 0.621138)

hierarchy_main = pd.read_csv(ROOT / "figures/main/Figure2/data/hierarchy_associations.csv")
hierarchy_s1 = pd.read_csv(ROOT / "figures/supplement/Figure S1/Source_S1a_Hierarchy.csv")
columns = ["system", "hierarchy", "rho", "ci_lo", "ci_hi", "n_rows", "n_clusters"]
left = hierarchy_main[columns].sort_values(["system", "hierarchy"]).reset_index(drop=True)
right = hierarchy_s1[columns].copy()
right["system"] = right["system"].map({
    "B5": "Cardiopedia organoids",
    "B6": "hiPSC-CM pharmacology",
})
right = right[columns].sort_values(["system", "hierarchy"]).reset_index(drop=True)
pd.testing.assert_frame_equal(left, right, check_exact=True)

checks = {
    "organoid_ligands": int(len(organoids)),
    "paired_organoids": int(organoids["n_paired_organoids"].sum()),
    "organoid_spearman_rho": rho_organoid,
    "pharmacology_compounds": int(primary["compounds"]),
    "pharmacology_dose_groups": int(primary["dose_groups"]),
    "pharmacology_primary_rho": float(primary["rho"]),
    "adult_lv_shared_genes": int(len(genes)),
    "adult_lv_spearman_rho": rho_adult_lv,
    "sign_concordance": sign_concordance,
    "absolute_weighted_sign_concordance": weighted_concordance,
    "main_figure_2_hierarchy_equals_supplementary_figure_1": True,
    "status": "PASS",
}
print(json.dumps(checks, indent=2))
