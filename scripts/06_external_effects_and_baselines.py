"""Summarise frozen external scRNA projection and benchmark simple/null scores."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "tables"
SEED = 20260711


def ci_and_exact(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float]:
    effect = float(a.mean() - b.mean())
    if len(a) < 2 or len(b) < 2:
        return effect, np.nan, np.nan, np.nan
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se = math.sqrt(va / len(a) + vb / len(b))
    df = (va / len(a) + vb / len(b)) ** 2 / (
        (va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
    q = stats.t.ppf(.975, df)
    pooled = np.concatenate([a, b])
    null = []
    for idx in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), bool); mask[list(idx)] = True
        null.append(pooled[mask].mean() - pooled[~mask].mean())
    p = float(np.mean(np.asarray(null) >= effect - 1e-12))
    return effect, float(effect - q * se), float(effect + q * se), p


def donor_from_cells(dataset: str) -> pd.DataFrame:
    cells = pd.read_csv(T / f"external_scrna_{dataset}_frozen_cell_scores.csv.gz")
    numeric = [c for c in cells.columns if c.startswith("frozen_z__") or c in [
        "fssi_frozen_weighted", "fssi_frozen_equal_weight", "pathological_component",
        "repair_component_reversed", "Fib5_three_marker_baseline"]]
    donor = cells.groupby(["dataset_id", "donor_id", "condition"], as_index=False)[numeric].mean()
    counts = cells.groupby(["dataset_id", "donor_id", "condition"], as_index=False).size().rename(columns={"size": "fibroblast_cell_n"})
    return donor.merge(counts, on=["dataset_id", "donor_id", "condition"])


def main() -> None:
    g156 = donor_from_cells("GSE156326")
    g181 = donor_from_cells("GSE181316")
    g130 = pd.read_csv(T / "external_scrna_GSE130973_frozen_donor_scores.csv")
    donor = pd.concat([g130, g156, g181], ignore_index=True, sort=False)
    donor.to_csv(T / "external_scrna_frozen_donor_scores_all.csv", index=False)

    contrasts = [
        ("GSE156326", "hypertrophic_scar", "normal_skin", "hypertrophic scar - normal skin"),
        ("GSE181316", "keloid", "normal_scar", "keloid - normal scar"),
        ("GSE181316", "keloid", "healthy_skin", "keloid - healthy skin (descriptive one-donor reference)"),
    ]
    scores = ["fssi_frozen_weighted", "fssi_frozen_equal_weight", "pathological_component",
              "repair_component_reversed", "Fib5_three_marker_baseline"]
    rows = []
    for dataset, positive, reference, label in contrasts:
        d = donor[donor.dataset_id.eq(dataset)]
        for score in [s for s in scores if s in d.columns]:
            a = d.loc[d.condition.eq(positive), score].dropna().to_numpy(float)
            b = d.loc[d.condition.eq(reference), score].dropna().to_numpy(float)
            effect, lo, hi, p = ci_and_exact(a, b)
            rows.append({"dataset_id": dataset, "contrast": label, "score": score,
                         "n_positive": len(a), "n_reference": len(b), "mean_difference": effect,
                         "ci95_low": lo, "ci95_high": hi, "exact_one_sided_p": p,
                         "inference": "descriptive" if len(b) < 2 else "donor_level_within_dataset"})
    effects = pd.DataFrame(rows)
    effects.to_csv(T / "external_scrna_frozen_effects_and_simple_baselines.csv", index=False)

    model = pd.read_csv(T / "fssi_frozen_model.csv")
    genes = model.gene.tolist()
    zcols = [f"frozen_z__{g}" for g in genes]
    weights = model.weight.to_numpy(float).copy(); weights /= np.abs(weights).sum()
    rng = np.random.default_rng(SEED)
    null_rows = []
    null_values = []
    for dataset, positive, reference, label in contrasts[:2]:
        d = donor[donor.dataset_id.eq(dataset)].copy()
        common_zcols = [c for c in zcols if c in d.columns and d[c].notna().all()]
        common_genes = [c.removeprefix("frozen_z__") for c in common_zcols]
        common_weights = model.set_index("gene").loc[common_genes, "weight"].to_numpy(float).copy()
        common_weights /= np.abs(common_weights).sum()
        matrix = d[common_zcols].to_numpy(float)
        group = d.condition.to_numpy()
        observed_score = matrix @ common_weights
        observed = observed_score[group == positive].mean() - observed_score[group == reference].mean()
        null = np.empty(5000)
        for i in range(len(null)):
            random_weights = rng.permutation(common_weights)
            score = matrix @ random_weights
            null[i] = score[group == positive].mean() - score[group == reference].mean()
        null_rows.append({"dataset_id": dataset, "contrast": label,
                          "null_type": "same_18_genes_random_weight_to_gene_assignment",
                          "permutations": len(null), "observed_effect": observed,
                          "genes_common": len(common_genes),
                          "null_mean": null.mean(), "null_sd": null.std(ddof=1),
                          "observed_percentile": np.mean(null <= observed),
                          "empirical_one_sided_p": (1 + np.sum(null >= observed)) / (len(null) + 1),
                          "seed": SEED})
        null_values.append(pd.DataFrame({"dataset_id": dataset, "contrast": label,
                                         "permutation": np.arange(1, len(null) + 1),
                                         "null_effect": null}))
    pd.DataFrame(null_rows).to_csv(T / "external_scrna_random_signature_benchmark.csv", index=False)
    pd.concat(null_values, ignore_index=True).to_csv(T / "external_scrna_random_signature_null_values.csv.gz",
                                                     index=False, compression="gzip")

    audits = []
    for dataset in ["GSE130973", "GSE151177", "GSE156326", "GSE181316"]:
        a = pd.read_csv(T / f"external_scrna_{dataset}_projection_audit.csv")
        audits.append(a)
    pd.concat(audits, ignore_index=True).to_csv(T / "external_scrna_projection_audit_all.csv", index=False)
    print(effects.to_string(index=False))
    print(pd.DataFrame(null_rows).to_string(index=False))


if __name__ == "__main__":
    main()
