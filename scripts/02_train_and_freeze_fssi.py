"""Train and freeze FSSI on development datasets at the biological-unit level."""

from __future__ import annotations

import json
import math
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
H5AD = ROOT / "data" / "processed" / "virtual_cell_inputs" / "fibroblast_reference_counts.h5ad"
SPEC = ROOT / "config" / "fssi_frozen_spec.yaml"
OUT = ROOT / "results" / "tables"
SD_FLOOR = 0.25


def role(disease: str) -> str:
    if disease in {"keloid", "hypertrophic_scar"}:
        return "pathological"
    if disease in {"normal_scar", "normal_or_adjacent"}:
        return "reference"
    return "excluded"


def normalized_expression(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    symbols = adata.var["gene_symbol"].astype(str).str.upper().to_numpy()
    lookup = {g: i for i, g in enumerate(symbols)}
    present = [g for g in genes if g in lookup]
    idx = [lookup[g] for g in present]
    totals = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    totals[totals <= 0] = 1.0
    selected = adata.X[:, idx]
    if sparse.issparse(selected):
        selected = selected.toarray()
    return np.log1p(np.asarray(selected, float) / totals[:, None] * 10_000), present


def sample_means(x: np.ndarray, obs: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for keys, idx in obs.groupby(["dataset_id", "sample_id", "disease_group"], sort=False).groups.items():
        loc = obs.index.get_indexer(idx)
        row = dict(zip(genes, np.mean(x[loc], axis=0)))
        row.update(dataset_id=keys[0], sample_id=keys[1], disease_group=keys[2], group=role(keys[2]), cells=len(loc))
        rows.append(row)
    return pd.DataFrame(rows)


def collapse_units(samples: pd.DataFrame, nf_pair: tuple[str, str] | None) -> pd.DataFrame:
    x = samples.copy()
    x["donor_id"] = x["dataset_id"] + "_" + x["sample_id"]
    if nf_pair:
        x.loc[x["sample_id"].isin(nf_pair), "donor_id"] = "GSE163973_NF_COLLAPSED_" + "_".join(nf_pair)
    x.loc[x["dataset_id"].eq("GSE243716"), "donor_id"] = "GSE243716_P1"
    numeric = [c for c in x.columns if c not in {"dataset_id", "sample_id", "disease_group", "group", "donor_id"}]
    units = x.groupby(["dataset_id", "donor_id", "group"], as_index=False)[numeric].mean()
    return units


def dataset_effects(units: pd.DataFrame, genes: list[str], scenario: str) -> pd.DataFrame:
    rows = []
    for dataset, sub in units.groupby("dataset_id"):
        a, b = sub[sub.group.eq("pathological")], sub[sub.group.eq("reference")]
        if a.empty or b.empty:
            continue
        for gene in genes:
            av, bv = a[gene].to_numpy(float), b[gene].to_numpy(float)
            pooled = np.concatenate([av, bv])
            scale = max(float(np.std(pooled, ddof=1)) if len(pooled) > 1 else 0.0, SD_FLOOR)
            rows.append({"scenario": scenario, "dataset_id": dataset, "gene": gene,
                         "n_pathological": len(av), "n_reference": len(bv),
                         "standardized_difference": float(np.clip((av.mean() - bv.mean()) / scale, -5, 5))})
    return pd.DataFrame(rows)


def main() -> None:
    spec = yaml.safe_load(SPEC.read_text(encoding="utf-8"))
    path_genes = [g.upper() for g in spec["pathological_candidate_universe"]]
    repair_genes = [g.upper() for g in spec["repair_candidate_universe"]]
    genes = path_genes + repair_genes
    adata = ad.read_h5ad(H5AD)
    obs = adata.obs.reset_index(drop=True).copy()
    x, present = normalized_expression(adata, genes)
    samples = sample_means(x, obs, present)

    scenarios = {
        "all_matrices_distinct": None,
        "collapse_NF1_NF2": ("NF1_matrix", "NF2_matrix"),
        "collapse_NF1_NF3": ("NF1_matrix", "NF3_matrix"),
        "collapse_NF2_NF3": ("NF2_matrix", "NF3_matrix"),
    }
    effect_tables = [dataset_effects(collapse_units(samples, pair), present, name) for name, pair in scenarios.items()]
    effects = pd.concat(effect_tables, ignore_index=True)
    expected = {g: 1.0 for g in path_genes} | {g: -1.0 for g in repair_genes}
    effects["expected_direction"] = effects.gene.map(expected)
    effects["directional_strength"] = effects.standardized_difference * effects.expected_direction

    stability = effects.groupby("gene", as_index=False).agg(
        min_directional_strength=("directional_strength", "min"),
        mean_standardized_difference=("standardized_difference", "mean"),
        informative_dataset_scenarios=("dataset_id", "size"),
    )
    stability["program"] = np.where(stability.gene.isin(path_genes), "pathological", "repair")
    eligible = stability[stability.min_directional_strength.gt(0)].copy()
    selected = pd.concat([
        eligible[eligible.program.eq("pathological")].nlargest(spec["training_rule"]["pathological_cap"], "min_directional_strength"),
        eligible[eligible.program.eq("repair")].nlargest(spec["training_rule"]["repair_cap"], "min_directional_strength"),
    ], ignore_index=True)
    if selected.program.nunique() != 2:
        raise RuntimeError("Strict selection did not retain both state programs.")

    selected["raw_weight"] = selected.mean_standardized_difference
    selected["weight"] = selected.raw_weight / selected.raw_weight.abs().sum()
    base_units = collapse_units(samples, None)
    train_mean = base_units[selected.gene].mean()
    train_sd = base_units[selected.gene].std(ddof=1).replace(0, 1.0).fillna(1.0)
    selected["training_mean"] = selected.gene.map(train_mean)
    selected["training_sd"] = selected.gene.map(train_sd)
    selected["equal_weight"] = selected.program.map({"pathological": 1.0, "repair": -1.0})
    selected["equal_weight"] /= selected.equal_weight.abs().sum()
    selected.insert(0, "model_version", spec["model_name"])
    selected.insert(1, "lock_date", str(spec["lock_date"]))

    arr = (base_units[selected.gene].to_numpy(float) - selected.training_mean.to_numpy()[None, :]) / selected.training_sd.to_numpy()[None, :]
    base_units["fssi_weighted"] = arr @ selected.weight.to_numpy(float)
    base_units["fssi_equal_weight"] = arr @ selected.equal_weight.to_numpy(float)
    base_units["inference_status"] = np.where(
        base_units.dataset_id.eq("GSE163973") & base_units.group.eq("reference"),
        "normal_donor_mapping_uncertain", "development_descriptive")

    OUT.mkdir(parents=True, exist_ok=True)
    selected.to_csv(OUT / "fssi_frozen_model.csv", index=False)
    effects.to_csv(OUT / "fssi_training_scenario_gene_effects.csv", index=False)
    stability.to_csv(OUT / "fssi_feature_stability.csv", index=False)
    base_units.to_csv(OUT / "fssi_development_donor_scores.csv", index=False)
    summary = {
        "model_version": spec["model_name"], "lock_date": str(spec["lock_date"]),
        "selected_genes": selected.gene.tolist(), "pathological_n": int((selected.program == "pathological").sum()),
        "repair_n": int((selected.program == "repair").sum()),
        "primary_external_validation_datasets": spec["locked_external_validation_datasets"],
        "validation_endpoint": "within-dataset donor-level effect; no correlation to the legacy SRI",
    }
    (OUT / "fssi_frozen_model_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
