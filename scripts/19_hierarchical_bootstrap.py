"""Hierarchical development bootstrap and descriptive PIEZO2 context analysis."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "tables"
SPEC = ROOT / "config" / "fssi_frozen_spec.yaml"
PIEZO = ROOT / ".stage" / "clinical" / "GSE266338" / "GSE266338_230313_FPKM.csv.gz"
SD_FLOOR = 0.25
N_BOOT = 2000
SEED = 20260712


def effect(sub: pd.DataFrame, genes: list[str]) -> np.ndarray:
    a = sub[sub.group.eq("pathological")][genes].to_numpy(float)
    b = sub[sub.group.eq("reference")][genes].to_numpy(float)
    scale = np.maximum(np.std(np.vstack([a, b]), axis=0, ddof=1), SD_FLOOR)
    return np.clip((a.mean(axis=0) - b.mean(axis=0)) / scale, -5, 5)


def bootstrap_mode(units: pd.DataFrame, genes: list[str], program: dict[str, str],
                   mode: str, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    informative = [d for d, x in units.groupby("dataset_id") if set(x.group) == {"pathological", "reference"}]
    records, iteration_rows = [], []
    gene_index = {g: i for i, g in enumerate(genes)}
    expected = np.array([1 if program[g] == "pathological" else -1 for g in genes])
    original = set(pd.read_csv(T / "fssi_frozen_model.csv").gene)
    for iteration in range(N_BOOT):
        selected_datasets = informative if mode == "donor" else list(rng.choice(informative, len(informative), replace=True))
        strengths, raw = [], []
        for dataset in selected_datasets:
            source = units[units.dataset_id.eq(dataset)]
            sampled = []
            for group, block in source.groupby("group"):
                idx = rng.choice(block.index.to_numpy(), len(block), replace=True)
                draw = block.loc[idx].copy(); draw["group"] = group; sampled.append(draw)
            values = effect(pd.concat(sampled), genes)
            raw.append(values); strengths.append(values * expected)
        strengths = np.vstack(strengths); raw = np.vstack(raw)
        minimum = strengths.min(axis=0); mean_raw = raw.mean(axis=0)
        eligible = minimum > 0
        selected = []
        for label, cap in [("pathological", 10), ("repair", 8)]:
            candidates = [g for g in genes if program[g] == label and eligible[gene_index[g]]]
            candidates.sort(key=lambda g: minimum[gene_index[g]], reverse=True)
            selected.extend(candidates[:cap])
        if selected:
            denom = sum(abs(mean_raw[gene_index[g]]) for g in selected)
            weights = {g: mean_raw[gene_index[g]] / denom for g in selected} if denom else {}
        else:
            weights = {}
        for gene in genes:
            records.append({"mode": mode, "iteration": iteration + 1, "gene": gene,
                            "selected": gene in selected, "weight": weights.get(gene, 0.0),
                            "weight_expected_sign": bool(weights.get(gene, 0.0) * expected[gene_index[gene]] > 0) if gene in selected else np.nan})
        iteration_rows.append({"mode": mode, "iteration": iteration + 1, "selected_gene_n": len(selected),
                               "original_gene_jaccard": len(set(selected) & original) / len(set(selected) | original) if selected else 0,
                               "all_selected_weights_expected_sign": all(weights[g] * expected[gene_index[g]] > 0 for g in selected) if selected else False,
                               "sampled_dataset_pattern": "+".join(selected_datasets)})
    return pd.DataFrame(records), pd.DataFrame(iteration_rows)


def piezo_context(model: pd.DataFrame) -> None:
    expression = pd.read_csv(PIEZO, compression="gzip").set_index("gene_id")
    expression.index = expression.index.astype(str).str.upper()
    expression = expression.groupby(level=0).median()
    genes = [g for g in model.gene if g in expression.index]
    logx = np.log2(expression.loc[genes] + 0.1)
    delta = logx["PIEZO2_Posi"] - logx["PIEZO2_Nega"]
    weights = model.set_index("gene").loc[genes, "weight"]
    contribution = delta * weights
    rows = pd.DataFrame({"gene": genes, "program": model.set_index("gene").loc[genes, "program"],
                         "log2_fpkm_pos_minus_neg": delta, "frozen_weight": weights,
                         "weighted_directional_contribution": contribution})
    rows["expected_direction_support"] = np.where(rows.program.eq("pathological"), rows.log2_fpkm_pos_minus_neg > 0,
                                                   rows.log2_fpkm_pos_minus_neg < 0)
    rows.to_csv(T / "GSE266338_PIEZO2_fibroblast_gene_context.csv", index=False)
    pd.DataFrame([{"dataset_id": "GSE266338", "profiles_per_group": 1, "genes_present": len(genes),
                   "weighted_log2_difference": float(contribution.sum()),
                   "expected_direction_gene_fraction": float(rows.expected_direction_support.mean()),
                   "interpretation": "descriptive_single_profile_per_group_not_validation"}]).to_csv(
                       T / "GSE266338_PIEZO2_fibroblast_summary.csv", index=False)


def main() -> None:
    spec = yaml.safe_load(SPEC.read_text(encoding="utf-8"))
    path_genes = [g.upper() for g in spec["pathological_candidate_universe"]]
    repair_genes = [g.upper() for g in spec["repair_candidate_universe"]]
    genes = path_genes + repair_genes
    program = {g: "pathological" for g in path_genes} | {g: "repair" for g in repair_genes}
    units = pd.read_csv(T / "fssi_development_donor_scores.csv")
    rng = np.random.default_rng(SEED)
    all_records, all_iterations = [], []
    for mode in ["donor", "hierarchical_dataset_and_donor"]:
        records, iterations = bootstrap_mode(units, genes, program, mode, rng)
        all_records.append(records); all_iterations.append(iterations)
    records = pd.concat(all_records, ignore_index=True)
    iterations = pd.concat(all_iterations, ignore_index=True)
    summary = records.groupby(["mode", "gene"], as_index=False).agg(
        selection_frequency=("selected", "mean"), median_weight=("weight", "median"),
        weight_q025=("weight", lambda x: x.quantile(.025)), weight_q975=("weight", lambda x: x.quantile(.975)))
    summary["original_frozen_gene"] = summary.gene.isin(set(pd.read_csv(T / "fssi_frozen_model.csv").gene))
    summary["program"] = summary.gene.map(program)
    summary.to_csv(T / "hierarchical_bootstrap_gene_stability.csv", index=False)
    iterations.to_csv(T / "hierarchical_bootstrap_iteration_summary.csv", index=False)
    records.to_csv(T / "hierarchical_bootstrap_gene_iterations.csv.gz", index=False, compression="gzip")
    piezo_context(pd.read_csv(T / "fssi_frozen_model.csv"))
    print(summary[summary.original_frozen_gene].sort_values(["mode", "selection_frequency"]).to_string(index=False))


if __name__ == "__main__":
    main()
