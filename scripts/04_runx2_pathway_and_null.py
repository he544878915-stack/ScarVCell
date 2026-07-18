"""Extend measured siRUNX2 validation with gene, pathway and matched-null analyses.

The experiment is small (three paired cultures), so pathway statistics are
reported as sensitivity evidence and are never used to refit the frozen FSSI.
"""

from __future__ import annotations

import itertools
import os
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
INPUT = Path(os.environ.get("SCARVCELL_STAGE", str(ROOT / ".stage"))) / "GSE293677" / "GSE293677_TPM.tab.gz"
TABLES = ROOT / "results" / "tables"
GMT_DIR = ROOT / "data" / "reference" / "msigdb"
SEED = 20260711


def paired_gene_statistics(expression: pd.DataFrame) -> pd.DataFrame:
    log_expr = np.log2(expression + 0.1)
    ctrl = log_expr[[f"CTRL{i}.TPM" for i in (1, 2, 3)]].to_numpy()
    kd = log_expr[[f"KD{i}.TPM" for i in (1, 2, 3)]].to_numpy()
    change = kd - ctrl
    mean = change.mean(axis=1)
    sd = change.std(axis=1, ddof=1)
    t_stat = np.divide(mean, sd / np.sqrt(3), out=np.zeros_like(mean), where=sd > 0)
    p = 2 * stats.t.sf(np.abs(t_stat), df=2)
    return pd.DataFrame({
        "gene": expression.index,
        "mean_log2_change_siRUNX2_minus_control": mean,
        "paired_t_statistic": t_stat,
        "paired_t_p": p,
        "mean_control_tpm": expression[[f"CTRL{i}.TPM" for i in (1, 2, 3)]].mean(axis=1),
        "mean_siRUNX2_tpm": expression[[f"KD{i}.TPM" for i in (1, 2, 3)]].mean(axis=1),
    }).reset_index(drop=True).sort_values("paired_t_statistic")


def score_reduction(expression: pd.DataFrame, genes: list[str], weights: np.ndarray) -> tuple[float, np.ndarray]:
    genes = [g for g in genes if g in expression.index]
    weights = np.asarray(weights, float)[: len(genes)]
    weights = weights / np.abs(weights).sum()
    x = np.log1p(expression.loc[genes].T)
    z = (x - x.mean()) / x.std(ddof=1).replace(0, 1).fillna(1)
    score = z.to_numpy() @ weights
    paired_reduction = score[:3] - score[3:]
    return float(paired_reduction.mean()), paired_reduction


def component_and_baseline_scores(expression: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    definitions = {
        "frozen_FSSI_weighted": (model.gene.tolist(), model.weight.to_numpy(float)),
        "frozen_FSSI_equal_weight": (model.gene.tolist(), model.equal_weight.to_numpy(float)),
        "pathological_component": (model.loc[model.program.eq("pathological"), "gene"].tolist(), np.ones(10)),
        "repair_component_reversed": (model.loc[model.program.eq("repair"), "gene"].tolist(), -np.ones(8)),
        "Fib_5_three_marker_baseline": (["ADAM12", "COMP", "POSTN"], np.ones(3)),
        "collagen_ECM_baseline": (["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL12A1", "POSTN", "FN1"], np.ones(7)),
    }
    rows = []
    for name, (genes, weights) in definitions.items():
        present = [g for g in genes if g in expression.index]
        reduction, pair_values = score_reduction(expression, present, weights)
        rows.append({
            "score": name,
            "genes_present": len(present),
            "mean_reduction_control_minus_siRUNX2": reduction,
            "pair_1_reduction": pair_values[0],
            "pair_2_reduction": pair_values[1],
            "pair_3_reduction": pair_values[2],
            "direction_support_fraction": float(np.mean(pair_values > 0)),
        })
    return pd.DataFrame(rows)


def matched_random_null(expression: pd.DataFrame, model: pd.DataFrame, permutations: int = 5000) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    log_expr = np.log1p(expression)
    feature = pd.DataFrame({"mean": log_expr.mean(axis=1), "sd": log_expr.std(axis=1, ddof=1)})
    feature = feature[(feature["mean"] > 0) & (feature["sd"] > 0)].copy()
    feature["mean_bin"] = pd.qcut(feature["mean"], 10, labels=False, duplicates="drop")
    feature["sd_bin"] = pd.qcut(feature["sd"], 5, labels=False, duplicates="drop")
    excluded = set(model.gene)
    candidates = feature.loc[~feature.index.isin(excluded)]
    target_bins = feature.reindex(model.gene)[["mean_bin", "sd_bin"]]
    weights = model.weight.to_numpy(float)
    observed, _ = score_reduction(expression, model.gene.tolist(), weights)
    values = np.empty(permutations)
    for i in range(permutations):
        sampled = []
        used: set[str] = set()
        for gene, row in target_bins.iterrows():
            pool = candidates[(candidates.mean_bin.eq(row.mean_bin)) & (candidates.sd_bin.eq(row.sd_bin))]
            pool = pool.loc[~pool.index.isin(used)]
            if pool.empty:
                pool = candidates.loc[~candidates.index.isin(used)]
            chosen = str(rng.choice(pool.index.to_numpy()))
            sampled.append(chosen)
            used.add(chosen)
        values[i], _ = score_reduction(expression, sampled, weights)
    summary = pd.DataFrame([{
        "null_type": "expression_and_variability_matched_18_gene_signed_weight_null",
        "permutations": permutations,
        "observed_frozen_FSSI_reduction": observed,
        "null_mean": float(values.mean()),
        "null_sd": float(values.std(ddof=1)),
        "empirical_one_sided_p": float((1 + np.sum(values >= observed)) / (permutations + 1)),
        "observed_percentile": float(np.mean(values <= observed)),
        "seed": SEED,
        "scope": "within_experiment_sensitivity_not_external_validation",
    }])
    pd.DataFrame({"random_signature_id": np.arange(1, permutations + 1),
                  "mean_reduction_control_minus_siRUNX2": values}).to_csv(
        TABLES / "runx2_matched_random_signature_null.csv", index=False)
    return summary


def run_preranked_gsea(gene_stats: pd.DataFrame) -> pd.DataFrame:
    # Positive ranks mean higher after siRUNX2; negative NES therefore indicates
    # suppression after knockdown. Duplicate symbols have already been collapsed.
    ranking = gene_stats[["gene", "paired_t_statistic"]].dropna().copy()
    # Deterministic, negligible tie breaking prevents arbitrary ordering of the
    # many genes with exactly zero paired variance and zero mean change.
    ranking = ranking.sort_values(["paired_t_statistic", "gene"], ascending=[False, True]).reset_index(drop=True)
    ranking["paired_t_statistic"] += np.linspace(1e-10, 0, len(ranking), endpoint=False)
    outputs = []
    for gmt in [
        GMT_DIR / "h.all.v2025.1.Hs.symbols.gmt",
        GMT_DIR / "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
        GMT_DIR / "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    ]:
        if not gmt.exists() or gmt.stat().st_size == 0:
            continue
        result = gp.prerank(
            rnk=ranking,
            gene_sets=str(gmt),
            permutation_num=1000,
            min_size=10,
            max_size=500,
            seed=SEED,
            threads=4,
            outdir=None,
            verbose=False,
        ).res2d.copy()
        result["collection_file"] = gmt.name
        outputs.append(result)
    all_results = pd.concat(outputs, ignore_index=True)
    all_results.columns = [str(c).lower().replace(" ", "_").replace("%", "pct") for c in all_results.columns]
    all_results.to_csv(TABLES / "runx2_preranked_pathway_all.csv", index=False)
    term = all_results["term"].astype(str).str.upper()
    focus = all_results[term.str.contains("COLLAGEN|EXTRACELLULAR_MATRIX|ECM_|FOCAL_ADHESION|TGF_BETA|EPITHELIAL_MESENCHYMAL|INTEGRIN", regex=True)].copy()
    focus.to_csv(TABLES / "runx2_preranked_pathway_focused.csv", index=False)
    return focus


def main() -> None:
    raw = pd.read_csv(INPUT, sep="\t", compression="gzip")
    raw["gene"] = raw["FEATURE_NAME"].astype(str).str.split("|").str[-1].str.upper()
    sample_cols = [f"CTRL{i}.TPM" for i in (1, 2, 3)] + [f"KD{i}.TPM" for i in (1, 2, 3)]
    expression = raw.groupby("gene")[sample_cols].median()
    model = pd.read_csv(TABLES / "fssi_frozen_model.csv")

    gene_stats = paired_gene_statistics(expression)
    gene_stats["frozen_FSSI_gene"] = gene_stats.gene.isin(model.gene)
    gene_stats.to_csv(TABLES / "runx2_paired_gene_response_all.csv", index=False)
    gene_stats[gene_stats.frozen_FSSI_gene].merge(
        model[["gene", "program", "weight"]], on="gene", how="left"
    ).to_csv(TABLES / "runx2_frozen_gene_response.csv", index=False)

    components = component_and_baseline_scores(expression, model)
    components.to_csv(TABLES / "runx2_component_and_simple_baseline_scores.csv", index=False)
    null_summary = matched_random_null(expression, model)
    null_summary.to_csv(TABLES / "runx2_matched_random_signature_summary.csv", index=False)
    focus = run_preranked_gsea(gene_stats)
    print(components.to_string(index=False))
    print(null_summary.to_string(index=False))
    print(focus.sort_values("fdr_q-val").head(20).to_string(index=False))


if __name__ == "__main__":
    main()
