"""Re-estimate external effects after droplet/doublet QC and common-gene harmonisation."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def independent_inference(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    effect = float(a.mean() - b.mean())
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se = math.sqrt(va / len(a) + vb / len(b))
    df = (va / len(a) + vb / len(b)) ** 2 / (
        (va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1)
    )
    critical = stats.t.ppf(0.975, df)
    pooled = np.concatenate([a, b])
    null = []
    for indices in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(indices)] = True
        null.append(pooled[mask].mean() - pooled[~mask].mean())
    null = np.asarray(null)
    return {
        "mean_difference": effect,
        "ci95_low": float(effect - critical * se),
        "ci95_high": float(effect + critical * se),
        "exact_two_sided_p": float(np.mean(np.abs(null) >= abs(effect) - 1e-12)),
        "exact_directional_p": float(np.mean(null >= effect - 1e-12)),
    }


def gse181316_qc_sensitivity() -> None:
    scores = pd.read_csv(TABLES / "external_scrna_GSE181316_frozen_cell_scores.csv.gz")
    qc = pd.read_csv(TABLES / "external_scrna_GSE181316_locked_fibroblast_qc.csv.gz")
    merge_keys = ["dataset_id", "sample_id", "donor_id", "condition", "cell_id"]
    data = scores.merge(qc, on=merge_keys, how="inner", validate="one_to_one")
    if len(data) != len(scores):
        raise RuntimeError("QC table does not cover every locked GSE181316 fibroblast")
    data["relaxed_marker_purity"] = (
        data.fibroblast_marker_detected_n.ge(3)
        & data.immune_marker_detected_n.le(1)
        & data.epithelial_marker_detected_n.le(1)
        & data.endothelial_marker_detected_n.le(1)
        & data.neural_marker_detected_n.le(1)
    )
    subsets = {
        "locked_all": np.ones(len(data), dtype=bool),
        "emptydrops_fdr01": data.emptydrops_called_fdr01.to_numpy(bool),
        "scDblFinder_singlet": data.scDblFinder_class.eq("singlet").to_numpy(),
        "emptydrops_and_singlet": (
            data.emptydrops_called_fdr01 & data.scDblFinder_class.eq("singlet")
        ).to_numpy(),
        "emptydrops_singlet_relaxed_marker_purity": (
            data.emptydrops_called_fdr01
            & data.scDblFinder_class.eq("singlet")
            & data.relaxed_marker_purity
        ).to_numpy(),
    }
    effect_rows, donor_rows = [], []
    for name, keep in subsets.items():
        subset = data.loc[keep].copy()
        donor = subset.groupby(["donor_id", "condition"], as_index=False).agg(
            fssi=("fssi_frozen_weighted", "mean"),
            retained_cells=("cell_id", "size"),
        )
        donor["qc_subset"] = name
        donor_rows.append(donor)
        a = donor.loc[donor.condition.eq("keloid"), "fssi"].to_numpy(float)
        b = donor.loc[donor.condition.eq("normal_scar"), "fssi"].to_numpy(float)
        inference = independent_inference(a, b)
        effect_rows.append({
            "dataset_id": "GSE181316", "contrast": "keloid - normal scar",
            "qc_subset": name, "retained_locked_fibroblasts": len(subset),
            "retained_fraction": len(subset) / len(data), "n_keloid_donors": len(a),
            "n_normal_scar_donors": len(b), **inference,
        })
    pd.concat(donor_rows, ignore_index=True).to_csv(
        TABLES / "external_scrna_GSE181316_qc_filtered_donor_scores.csv", index=False)
    pd.DataFrame(effect_rows).to_csv(
        TABLES / "external_scrna_GSE181316_qc_filtered_effects.csv", index=False)

    marker_columns = [
        "fibroblast_marker_detected_n", "immune_marker_detected_n",
        "epithelial_marker_detected_n", "endothelial_marker_detected_n",
        "neural_marker_detected_n",
    ]
    marker = data.groupby("condition")[marker_columns].agg(["median", "mean"])
    marker.columns = ["__".join(column) for column in marker.columns]
    marker = marker.reset_index()
    purity = data.groupby("condition", as_index=False).agg(
        cells=("cell_id", "size"),
        emptydrops_support_fraction=("emptydrops_called_fdr01", "mean"),
        singlet_fraction=("scDblFinder_class", lambda x: float(np.mean(x == "singlet"))),
        relaxed_marker_purity_fraction=("relaxed_marker_purity", "mean"),
    )
    purity.merge(marker, on="condition").to_csv(
        TABLES / "external_scrna_GSE181316_marker_purity_summary.csv", index=False)


def gse156326_common_gene_sensitivity() -> None:
    cells = pd.read_csv(TABLES / "external_scrna_GSE156326_frozen_cell_scores.csv.gz")
    model = pd.read_csv(TABLES / "fssi_frozen_model.csv").set_index("gene")
    z_columns = [column for column in cells if column.startswith("frozen_z__")]
    common_columns = [column for column in z_columns if cells[column].notna().all()]
    common_genes = [column.removeprefix("frozen_z__") for column in common_columns]
    weights = model.loc[common_genes, "weight"].to_numpy(float).copy()
    weights /= np.abs(weights).sum()
    cells["fssi_common_gene"] = cells[common_columns].to_numpy(float) @ weights
    donor = cells.groupby(["donor_id", "condition"], as_index=False).agg(
        fssi_common_gene=("fssi_common_gene", "mean"),
        fssi_original=("fssi_frozen_weighted", "mean"),
        fibroblast_cells=("cell_id", "size"),
    )
    donor.to_csv(TABLES / "external_scrna_GSE156326_common_gene_donor_scores.csv", index=False)
    rows = []
    for score in ["fssi_original", "fssi_common_gene"]:
        a = donor.loc[donor.condition.eq("hypertrophic_scar"), score].to_numpy(float)
        b = donor.loc[donor.condition.eq("normal_skin"), score].to_numpy(float)
        rows.append({
            "dataset_id": "GSE156326", "contrast": "hypertrophic scar - normal skin",
            "score": score, "genes_used": 18 if score == "fssi_original" else len(common_genes),
            "common_genes": ";".join(common_genes) if score == "fssi_common_gene" else "mixed 17/18 by sample",
            **independent_inference(a, b),
        })
    pd.DataFrame(rows).to_csv(
        TABLES / "external_scrna_GSE156326_common_gene_sensitivity.csv", index=False)
    coverage = cells.groupby(["sample_id", "condition"], as_index=False).agg(
        fibroblast_cells=("cell_id", "size"), genes_present=("genes_present", "first")
    )
    coverage.to_csv(TABLES / "external_scrna_GSE156326_gene_coverage_audit.csv", index=False)


def main() -> None:
    gse181316_qc_sensitivity()
    gse156326_common_gene_sensitivity()
    print(pd.read_csv(TABLES / "external_scrna_GSE181316_qc_filtered_effects.csv").to_string(index=False))
    print(pd.read_csv(TABLES / "external_scrna_GSE156326_common_gene_sensitivity.csv").to_string(index=False))


if __name__ == "__main__":
    main()
