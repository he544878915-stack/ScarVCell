"""Validate the frozen programme in the metadata-locked siPOSTN experiment."""

from __future__ import annotations

import itertools
import math
import os
from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
INPUT = Path(os.environ.get("SCARVCELL_STAGE", str(ROOT / ".stage"))) / "postlock_locked" / "GSE307210" / "GSE307210_gene_fpkm.txt.gz"
TABLES = ROOT / "results" / "tables"
GMT_DIR = ROOT / "data" / "reference" / "msigdb"
SEED = 20260712


def exact_unpaired(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    observed = float(a.mean() - b.mean())
    pooled = np.concatenate([a, b])
    null = []
    for idx in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        null.append(pooled[mask].mean() - pooled[~mask].mean())
    null = np.asarray(null)
    return float(np.mean(null >= observed - 1e-12)), float(np.mean(np.abs(null) >= abs(observed) - 1e-12))


def welch_interval(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    effect = float(a.mean() - b.mean())
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se2 = va / len(a) + vb / len(b)
    df = se2**2 / ((va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
    q = stats.t.ppf(0.975, df)
    return effect, float(effect - q * math.sqrt(se2)), float(effect + q * math.sqrt(se2))


def score_matrix(expression: pd.DataFrame, genes: list[str], weights: np.ndarray) -> pd.Series:
    present = [gene for gene in genes if gene in expression.index]
    selected = expression.loc[present]
    x = np.log1p(selected.T)
    z = (x - x.mean()) / x.std(ddof=1).replace(0, 1).fillna(1)
    weights = np.asarray(weights, dtype=float)[: len(present)].copy()
    weights /= np.abs(weights).sum()
    return pd.Series(z.to_numpy() @ weights, index=x.index)


def matched_null(expression: pd.DataFrame, model: pd.DataFrame, permutations: int = 5000) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    log_expression = np.log1p(expression)
    feature = pd.DataFrame({"mean": log_expression.mean(axis=1), "sd": log_expression.std(axis=1, ddof=1)})
    feature = feature[(feature["mean"] > 0) & (feature["sd"] > 0)].copy()
    feature["mean_bin"] = pd.qcut(feature["mean"], 10, labels=False, duplicates="drop")
    feature["sd_bin"] = pd.qcut(feature["sd"], 5, labels=False, duplicates="drop")
    candidates = feature.loc[~feature.index.isin(model.gene)]
    target_bins = feature.reindex(model.gene)[["mean_bin", "sd_bin"]]
    observed_scores = score_matrix(expression, model.gene.tolist(), model.weight.to_numpy(float))
    observed = float(observed_scores.iloc[:3].mean() - observed_scores.iloc[3:].mean())
    values = np.empty(permutations)
    for i in range(permutations):
        selected: list[str] = []
        used: set[str] = set()
        for _, row in target_bins.iterrows():
            pool = candidates[candidates.mean_bin.eq(row.mean_bin) & candidates.sd_bin.eq(row.sd_bin)]
            pool = pool.loc[~pool.index.isin(used)]
            if pool.empty:
                pool = candidates.loc[~candidates.index.isin(used)]
            gene = str(rng.choice(pool.index.to_numpy()))
            selected.append(gene)
            used.add(gene)
        random_scores = score_matrix(expression, selected, model.weight.to_numpy(float))
        values[i] = random_scores.iloc[:3].mean() - random_scores.iloc[3:].mean()
    summary = pd.DataFrame([{
        "null_type": "expression_and_variability_matched_18_gene_signed_weight_null",
        "permutations": permutations,
        "observed_control_minus_siPOSTN": observed,
        "null_mean": float(values.mean()),
        "null_sd": float(values.std(ddof=1)),
        "empirical_one_sided_p": float((1 + np.sum(values >= observed)) / (permutations + 1)),
        "observed_percentile": float(np.mean(values <= observed)),
        "seed": SEED,
    }])
    values_df = pd.DataFrame({"random_signature_id": np.arange(1, permutations + 1), "control_minus_siPOSTN": values})
    return summary, values_df


def pathway_response(expression: pd.DataFrame) -> pd.DataFrame:
    log_expression = np.log2(expression + 0.1)
    ctrl = log_expression.iloc[:, :3]
    kd = log_expression.iloc[:, 3:]
    statistic = stats.ttest_ind(ctrl.to_numpy(), kd.to_numpy(), axis=1, equal_var=False).statistic
    ranking = pd.DataFrame({"gene": expression.index, "welch_t_control_minus_siPOSTN": np.nan_to_num(statistic)})
    ranking = ranking.sort_values(["welch_t_control_minus_siPOSTN", "gene"], ascending=[False, True]).reset_index(drop=True)
    ranking["welch_t_control_minus_siPOSTN"] += np.linspace(1e-10, 0, len(ranking), endpoint=False)
    outputs = []
    for gmt in [
        GMT_DIR / "h.all.v2025.1.Hs.symbols.gmt",
        GMT_DIR / "c2.cp.reactome.v2025.1.Hs.symbols.gmt",
        GMT_DIR / "c2.cp.kegg_legacy.v2025.1.Hs.symbols.gmt",
    ]:
        if not gmt.exists() or not gmt.stat().st_size:
            continue
        result = gp.prerank(rnk=ranking, gene_sets=str(gmt), permutation_num=1000,
                            min_size=10, max_size=500, seed=SEED, threads=4,
                            outdir=None, verbose=False).res2d.copy()
        result["collection_file"] = gmt.name
        outputs.append(result)
    all_results = pd.concat(outputs, ignore_index=True)
    all_results.columns = [str(c).lower().replace(" ", "_").replace("%", "pct") for c in all_results.columns]
    all_results.to_csv(TABLES / "postn_preranked_pathway_all.csv", index=False)
    term = all_results.term.astype(str).str.upper()
    focus = all_results[term.str.contains("COLLAGEN|EXTRACELLULAR_MATRIX|ECM_|FOCAL_ADHESION|TGF_BETA|INTEGRIN", regex=True)].copy()
    focus.to_csv(TABLES / "postn_preranked_pathway_focused.csv", index=False)
    return focus


def main() -> None:
    raw = pd.read_csv(INPUT, sep="\t", compression="gzip")
    sample_cols = ["A1", "A2", "A3", "B1", "B2", "B3"]
    raw["gene"] = raw.gene_name.astype(str).str.upper()
    expression = raw.groupby("gene")[sample_cols].median()
    model = pd.read_csv(TABLES / "fssi_frozen_model.csv")
    present = [gene for gene in model.gene if gene in expression.index]

    weighted = score_matrix(expression, present, model.set_index("gene").loc[present, "weight"].to_numpy(float))
    equal = score_matrix(expression, present, model.set_index("gene").loc[present, "equal_weight"].to_numpy(float))
    scores = pd.DataFrame({
        "dataset_id": "GSE307210",
        "sample_id": sample_cols,
        "treatment": ["siControl"] * 3 + ["siPOSTN"] * 3,
        "fixed_weight_response_score": weighted.to_numpy(),
        "equal_weight_response_score": equal.to_numpy(),
        "genes_present": len(present),
        "inference_unit": "unpaired_profile",
    })

    summaries = []
    for column in ["fixed_weight_response_score", "equal_weight_response_score"]:
        a = scores.loc[scores.treatment.eq("siControl"), column].to_numpy(float)
        b = scores.loc[scores.treatment.eq("siPOSTN"), column].to_numpy(float)
        effect, lo, hi = welch_interval(a, b)
        p1, p2 = exact_unpaired(a, b)
        summaries.append({
            "dataset_id": "GSE307210", "target": "POSTN", "perturbation": "siPOSTN",
            "score_variant": column, "n_control": len(a), "n_siPOSTN": len(b),
            "mean_control_minus_siPOSTN": effect, "ci95_low": lo, "ci95_high": hi,
            "exact_one_sided_p": p1, "exact_two_sided_p": p2,
            "inference": "measured_target_specific_unpaired_expression_response_small_n",
        })

    log_expression = np.log2(expression + 0.1)
    gene_rows = []
    for gene in expression.index:
        a = log_expression.loc[gene, sample_cols[:3]].to_numpy(float)
        b = log_expression.loc[gene, sample_cols[3:]].to_numpy(float)
        stat = stats.ttest_ind(a, b, equal_var=False)
        gene_rows.append({
            "gene": gene, "mean_log2_change_siPOSTN_minus_control": float(b.mean() - a.mean()),
            "welch_t_control_minus_siPOSTN": float(stat.statistic) if np.isfinite(stat.statistic) else 0.0,
            "welch_p": float(stat.pvalue) if np.isfinite(stat.pvalue) else 1.0,
            "mean_control_fpkm": float(expression.loc[gene, sample_cols[:3]].mean()),
            "mean_siPOSTN_fpkm": float(expression.loc[gene, sample_cols[3:]].mean()),
            "frozen_FSSI_gene": gene in set(model.gene),
        })
    genes = pd.DataFrame(gene_rows)

    definitions = {
        "frozen_FSSI_weighted": (model.gene.tolist(), model.weight.to_numpy(float)),
        "frozen_FSSI_equal_weight": (model.gene.tolist(), model.equal_weight.to_numpy(float)),
        "pathological_component": (model.loc[model.program.eq("pathological"), "gene"].tolist(), np.ones(10)),
        "repair_component_reversed": (model.loc[model.program.eq("repair"), "gene"].tolist(), -np.ones(8)),
        "Fib_5_three_marker_baseline": (["ADAM12", "COMP", "POSTN"], np.ones(3)),
        "collagen_ECM_baseline": (["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL12A1", "POSTN", "FN1"], np.ones(7)),
    }
    component_rows = []
    for name, (geneset, weights) in definitions.items():
        present_set = [gene for gene in geneset if gene in expression.index]
        component_score = score_matrix(expression, present_set, np.asarray(weights)[: len(present_set)])
        a, b = component_score.iloc[:3].to_numpy(), component_score.iloc[3:].to_numpy()
        component_rows.append({"score": name, "genes_present": len(present_set),
                               "mean_control_minus_siPOSTN": float(a.mean() - b.mean()),
                               "control_values": ";".join(f"{x:.6g}" for x in a),
                               "siPOSTN_values": ";".join(f"{x:.6g}" for x in b)})

    null_summary, null_values = matched_null(expression, model)
    focus = pathway_response(expression)
    scores.to_csv(TABLES / "postn_measured_perturbation_scores.csv", index=False)
    pd.DataFrame(summaries).to_csv(TABLES / "postn_measured_perturbation_summary.csv", index=False)
    genes.to_csv(TABLES / "postn_gene_response_all.csv", index=False)
    genes[genes.frozen_FSSI_gene].merge(model[["gene", "program", "weight"]], on="gene", how="left").to_csv(
        TABLES / "postn_frozen_gene_response.csv", index=False)
    pd.DataFrame(component_rows).to_csv(TABLES / "postn_component_and_simple_baseline_scores.csv", index=False)
    null_summary.to_csv(TABLES / "postn_matched_random_signature_summary.csv", index=False)
    null_values.to_csv(TABLES / "postn_matched_random_signature_null.csv.gz", index=False, compression="gzip")
    print(pd.DataFrame(summaries).to_string(index=False))
    print(genes[genes.gene.eq("POSTN")].to_string(index=False))
    print(pd.DataFrame(component_rows).to_string(index=False))
    print(null_summary.to_string(index=False))
    print(focus.sort_values("fdr_q-val").head(12).to_string(index=False))


if __name__ == "__main__":
    main()
