"""Derive and project the externally defined Greene et al. susceptibility anchor.

The source-study thresholds and tissues were locked before this script was run.
The resulting genes do not alter the FSSI model or its candidate-level tests.
"""

from __future__ import annotations

import importlib.util
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(os.environ.get("SCARVCELL_ROOT", Path(__file__).resolve().parents[1]))
TABLES = ROOT / "results" / "tables"
SOURCE = (ROOT / ".stage" / "genetics" / "Greene2025_supplement" /
          "41467_2025_62945_MOESM10_ESM.xlsx")
SHEET = "SD8a. TWAS coloc Multi"
RELEVANT_TISSUES = {
    "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg",
    "Cells_Cultured_fibroblasts",
}
GPGE_THRESHOLD = 1.5e-6
COLOC_THRESHOLD = 0.90


def load_gse_module():
    path = ROOT / "scripts" / "167_analyse_gse158395.py"
    spec = importlib.util.spec_from_file_location("gse158395", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def derive_anchor() -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_excel(SOURCE, sheet_name=SHEET, skiprows=1)
    source["p-value"] = pd.to_numeric(source["p-value"], errors="coerce")
    source["PP.H4.abf"] = pd.to_numeric(source["PP.H4.abf"], errors="coerce")
    source["Z-score"] = pd.to_numeric(source["Z-score"], errors="coerce")
    relevant = source[
        source.Tissue.isin(RELEVANT_TISSUES)
        & source["p-value"].lt(GPGE_THRESHOLD)
        & source["PP.H4.abf"].gt(COLOC_THRESHOLD)
    ].copy()
    relevant["source_study_threshold"] = "GPGE P<1.5e-6 and PP.H4>0.90"
    relevant.to_csv(TABLES / "Greene2025_relevant_tissue_colocalised_rows.csv", index=False)

    rows = []
    for gene, data in relevant.groupby("GeneName", sort=True):
        directions = set(np.sign(data["Z-score"].dropna()).astype(int)) - {0}
        selected = data.loc[data["Z-score"].abs().idxmax()]
        rows.append({
            "gene": gene,
            "selected_tissue": selected.Tissue,
            "gpge_z": selected["Z-score"],
            "gpge_p": selected["p-value"],
            "pp_h4": selected["PP.H4.abf"],
            "qualifying_tissues": data.Tissue.nunique(),
            "direction_conflict": len(directions) > 1,
        })
    anchor = pd.DataFrame(rows)
    anchor["included_in_signed_score"] = ~anchor.direction_conflict
    denominator = anchor.loc[anchor.included_in_signed_score, "gpge_z"].abs().sum()
    anchor["signed_weight"] = np.where(
        anchor.included_in_signed_score, anchor.gpge_z / denominator, np.nan)
    fssi = set(pd.read_csv(TABLES / "fssi_frozen_model.csv").gene)
    extended = {"RUNX2", "POSTN", "TGFB1", "CSF1", "CXCL12", "MIF", "APOD", "PI16", "PTGDS", "PDGFB"}
    anchor["in_fssi"] = anchor.gene.isin(fssi)
    anchor["in_extended_candidate_panel"] = anchor.gene.isin(extended)
    anchor.to_csv(TABLES / "Greene2025_locked_susceptibility_anchor.csv", index=False)
    return relevant, anchor


def paired_interval(values: np.ndarray) -> tuple[float, float, float, float]:
    values = np.asarray(values, float)
    effect = float(values.mean())
    se = values.std(ddof=1) / math.sqrt(len(values))
    q = stats.t.ppf(0.975, len(values) - 1)
    p_t = 2 * stats.t.sf(abs(effect / se), len(values) - 1) if se else 1.0
    observed = abs(effect)
    null = [abs(np.mean(values * np.asarray(signs)))
            for signs in __import__("itertools").product([-1, 1], repeat=len(values))]
    return effect, effect - q * se, effect + q * se, float(np.mean(np.asarray(null) >= observed - 1e-12))


def project_gse158395(anchor: pd.DataFrame) -> None:
    module = load_gse_module()
    expression = module.build_expression()
    meta = module.metadata()
    selected = anchor[anchor.included_in_signed_score & anchor.gene.isin(expression.index)].copy()
    genes = selected.gene.tolist()
    x = expression.loc[genes].T
    z = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1).replace(0, 1).fillna(1)
    weight = selected.set_index("gene").loc[genes, "signed_weight"]
    weight = weight / weight.abs().sum()
    scores = pd.DataFrame({
        "sample_id": z.index,
        "genetic_susceptibility_expression_score": z.to_numpy() @ weight.to_numpy(),
        "anchor_genes_present": len(genes),
        "anchor_genes_total": int(anchor.included_in_signed_score.sum()),
    }).merge(meta, on="sample_id", validate="one_to_one")
    fssi = pd.read_csv(TABLES / "GSE158395_locked_sample_scores.csv")
    scores = scores.merge(fssi[["sample_id", "fixed_weight_response"]], on="sample_id", validate="one_to_one")
    scores.to_csv(TABLES / "GSE158395_genetic_anchor_scores.csv", index=False)

    paired = scores[scores.condition.isin(["lesion", "nonlesion"])]
    wide = paired.pivot(index="participant_id", columns="condition",
                        values="genetic_susceptibility_expression_score").dropna()
    delta = wide.lesion - wide.nonlesion
    effect, low, high, exact_p = paired_interval(delta.to_numpy())
    fssi_wide = paired.pivot(index="participant_id", columns="condition",
                             values="fixed_weight_response").dropna()
    fssi_delta = fssi_wide.lesion - fssi_wide.nonlesion
    rho_sample, p_sample = stats.spearmanr(
        scores.genetic_susceptibility_expression_score, scores.fixed_weight_response)
    rho_delta, p_delta = stats.spearmanr(delta, fssi_delta.loc[delta.index])
    pd.DataFrame([{
        "dataset_id": "GSE158395",
        "contrast": "lesion_minus_matched_nonlesion",
        "participant_pairs": len(delta),
        "anchor_genes_present": len(genes),
        "mean_difference": effect,
        "ci95_low": low,
        "ci95_high": high,
        "exact_signflip_p_two_sided": exact_p,
        "positive_pairs": int((delta > 0).sum()),
        "sample_level_spearman_with_fssi": rho_sample,
        "sample_level_spearman_p": p_sample,
        "paired_delta_spearman_with_fssi_delta": rho_delta,
        "paired_delta_spearman_p": p_delta,
    }]).to_csv(TABLES / "GSE158395_genetic_anchor_effect.csv", index=False)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    _, anchor = derive_anchor()
    project_gse158395(anchor)
    print(anchor.to_string(index=False))
    print(pd.read_csv(TABLES / "GSE158395_genetic_anchor_effect.csv").to_string(index=False))


if __name__ == "__main__":
    main()
