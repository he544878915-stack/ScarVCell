"""Condition-blind matched-gene nulls and unified external score benchmarks."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
SEED = 20260711
PERMUTATIONS = 5000

CONTRASTS = {
    "GSE156326": ("hypertrophic_scar", "normal_skin", "hypertrophic scar - normal skin"),
    "GSE181316": ("keloid", "normal_scar", "keloid - normal scar"),
}


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    n1, n0 = len(a), len(b)
    pooled = np.sqrt(((n1 - 1) * np.var(a, ddof=1) + (n0 - 1) * np.var(b, ddof=1)) / (n1 + n0 - 2))
    if pooled == 0:
        return np.nan
    d = (a.mean() - b.mean()) / pooled
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(correction * d)


def load_dataset(dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = TABLES / f"external_scrna_{dataset}_locked_fibroblast_donor_gene_expression.csv.gz"
    long = pd.read_csv(path)
    metadata = long[["donor_id", "condition"]].drop_duplicates().set_index("donor_id")
    matrix = long.pivot(index="donor_id", columns="gene", values="mean_log1p_cp10k")
    matrix = matrix.reindex(metadata.index)
    return matrix, metadata


def donor_z(matrix: pd.DataFrame) -> pd.DataFrame:
    sd = matrix.std(axis=0, ddof=1)
    keep = sd[(sd > 0) & np.isfinite(sd)].index
    matrix = matrix[keep].fillna(0)
    return (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=1)


def signed_score(z: pd.DataFrame, genes: list[str], weights: np.ndarray) -> np.ndarray:
    present = [gene for gene in genes if gene in z.columns]
    position = {gene: i for i, gene in enumerate(genes)}
    selected_weights = np.asarray([weights[position[gene]] for gene in present], float)
    selected_weights /= np.abs(selected_weights).sum()
    return z[present].to_numpy(float) @ selected_weights


def effect(score: np.ndarray, metadata: pd.DataFrame, positive: str, reference: str) -> tuple[float, float]:
    labels = metadata.condition.to_numpy()
    a, b = score[labels == positive], score[labels == reference]
    return float(a.mean() - b.mean()), hedges_g(a, b)


def matched_null(dataset: str, matrix: pd.DataFrame, metadata: pd.DataFrame,
                 model: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    positive, reference, contrast = CONTRASTS[dataset]
    z = donor_z(matrix)
    model = model[model.gene.isin(z.columns)].copy()
    genes = model.gene.tolist()
    weights = model.weight.to_numpy(float)
    observed_score = signed_score(z, genes, weights)
    observed_effect, observed_g = effect(observed_score, metadata, positive, reference)

    feature = pd.DataFrame({
        "mean": matrix.mean(axis=0),
        "sd": matrix.std(axis=0, ddof=1),
    }).replace([np.inf, -np.inf], np.nan).dropna()
    candidates = feature[(feature["mean"] > 0.05) & (feature["sd"] > 0)]
    candidates = candidates.loc[~candidates.index.isin(model.gene)].copy()
    scale = candidates.std(axis=0, ddof=1).replace(0, 1)
    rng = np.random.default_rng(SEED + (1 if dataset == "GSE181316" else 0))
    nearest: dict[str, np.ndarray] = {}
    for gene in genes:
        target = feature.loc[gene]
        distance = (((candidates - target) / scale) ** 2).sum(axis=1)
        nearest[gene] = distance.nsmallest(min(200, len(distance))).index.to_numpy()

    null_effect = np.empty(PERMUTATIONS)
    null_g = np.empty(PERMUTATIONS)
    rows = []
    for iteration in range(PERMUTATIONS):
        sampled: list[str] = []
        used: set[str] = set()
        for gene in genes:
            pool = np.asarray([item for item in nearest[gene] if item not in used])
            if not len(pool):
                pool = candidates.index.difference(list(used)).to_numpy()
            chosen = str(rng.choice(pool))
            sampled.append(chosen)
            used.add(chosen)
        score = signed_score(z, sampled, weights)
        null_effect[iteration], null_g[iteration] = effect(score, metadata, positive, reference)
        rows.append({"dataset_id": dataset, "signature_id": iteration + 1,
                     "matched_genes_in_frozen_weight_order": ";".join(sampled),
                     "mean_difference": null_effect[iteration], "hedges_g": null_g[iteration]})
    summary = {
        "dataset_id": dataset,
        "contrast": contrast,
        "null_type": "condition_blind_expression_variability_matched_random_gene_sets",
        "matching_pool": "200 nearest genes by pooled donor mean and SD; disease labels unused",
        "genes_common": len(genes),
        "permutations": PERMUTATIONS,
        "observed_benchmark_scale_effect": observed_effect,
        "observed_hedges_g": observed_g,
        "effect_percentile": float(np.mean(null_effect <= observed_effect)),
        "hedges_g_percentile": float(np.mean(null_g <= observed_g)),
        "empirical_one_sided_p_effect": float((1 + np.sum(null_effect >= observed_effect)) / (PERMUTATIONS + 1)),
        "empirical_one_sided_p_hedges_g": float((1 + np.sum(null_g >= observed_g)) / (PERMUTATIONS + 1)),
        "seed": int(rng.bit_generator._seed_seq.entropy),
        "scope": "post_lock_condition_blind_sensitivity_benchmark",
    }
    return summary, pd.DataFrame(rows)


def unified_baselines(dataset: str, matrix: pd.DataFrame, metadata: pd.DataFrame,
                      model: pd.DataFrame) -> list[dict]:
    positive, reference, contrast = CONTRASTS[dataset]
    z = donor_z(matrix)
    definitions = {
        "FSSI": (model.gene.tolist(), model.weight.to_numpy(float)),
        "Pathological component": (model.loc[model.program.eq("pathological"), "gene"].tolist(), np.ones(10)),
        "Repair reversed": (model.loc[model.program.eq("repair"), "gene"].tolist(), -np.ones(8)),
        "Fib_5 markers": (["ADAM12", "COMP", "POSTN"], np.ones(3)),
        "Collagen/ECM": (["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL12A1", "POSTN", "FN1"], np.ones(7)),
    }
    rows = []
    for score_name, (genes, weights) in definitions.items():
        present = [gene for gene in genes if gene in z.columns]
        values = signed_score(z, genes, weights)
        mean_difference, g = effect(values, metadata, positive, reference)
        rows.append({
            "dataset_id": dataset, "contrast": contrast, "score": score_name,
            "genes_present": len(present), "mean_difference_benchmark_scale": mean_difference,
            "hedges_g": g, "direction_supported": bool(mean_difference > 0),
        })
    return rows


def main() -> None:
    model = pd.read_csv(TABLES / "fssi_frozen_model.csv")
    summaries, nulls, baselines = [], [], []
    for dataset in CONTRASTS:
        matrix, metadata = load_dataset(dataset)
        summary, null = matched_null(dataset, matrix, metadata, model)
        summaries.append(summary)
        nulls.append(null)
        baselines.extend(unified_baselines(dataset, matrix, metadata, model))
    pd.DataFrame(summaries).to_csv(TABLES / "external_scrna_matched_random_gene_set_summary.csv", index=False)
    pd.concat(nulls, ignore_index=True).to_csv(
        TABLES / "external_scrna_matched_random_gene_set_null.csv.gz", index=False, compression="gzip")
    baseline = pd.DataFrame(baselines)
    baseline.to_csv(TABLES / "external_scrna_unified_score_benchmark.csv", index=False)
    stability = baseline.groupby("score", as_index=False).agg(
        datasets_supported=("direction_supported", "sum"),
        minimum_hedges_g=("hedges_g", "min"),
        mean_hedges_g=("hedges_g", "mean"),
        hedges_g_range=("hedges_g", lambda x: float(x.max() - x.min())),
    )
    stability["datasets_total"] = len(CONTRASTS)
    stability.to_csv(TABLES / "external_scrna_unified_score_stability.csv", index=False)
    print(pd.DataFrame(summaries).to_string(index=False))
    print(baseline.to_string(index=False))
    print(stability.to_string(index=False))


if __name__ == "__main__":
    main()
