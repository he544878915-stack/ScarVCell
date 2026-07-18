"""Summarise the locked post-lock scRNA projections without cross-dataset pooling."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
SEED = 20260712
SCORES = [
    "fssi_frozen_weighted",
    "fssi_frozen_equal_weight",
    "pathological_component",
    "repair_component_reversed",
    "Fib5_three_marker_baseline",
]


def unpaired_summary(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float, float]:
    effect = float(a.mean() - b.mean())
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se2 = va / len(a) + vb / len(b)
    df = se2**2 / ((va / len(a)) ** 2 / (len(a) - 1) + (vb / len(b)) ** 2 / (len(b) - 1))
    q = stats.t.ppf(0.975, df)
    pooled = np.concatenate([a, b])
    null = []
    for idx in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx)] = True
        null.append(pooled[mask].mean() - pooled[~mask].mean())
    null = np.asarray(null)
    p1 = float(np.mean(null >= effect - 1e-12))
    p2 = float(np.mean(np.abs(null) >= abs(effect) - 1e-12))
    return effect, float(effect - q * math.sqrt(se2)), float(effect + q * math.sqrt(se2)), p1, p2


def paired_summary(difference: np.ndarray) -> tuple[float, float, float, float, float]:
    effect = float(difference.mean())
    se = stats.sem(difference)
    q = stats.t.ppf(0.975, len(difference) - 1)
    null = np.asarray([
        np.mean(difference * np.asarray(signs))
        for signs in itertools.product([-1, 1], repeat=len(difference))
    ])
    p1 = float(np.mean(null >= effect - 1e-12))
    p2 = float(np.mean(np.abs(null) >= abs(effect) - 1e-12))
    return effect, float(effect - q * se), float(effect + q * se), p1, p2


def random_weight_null(data: pd.DataFrame, contrast: str, permutations: int = 5000) -> tuple[dict[str, object], pd.DataFrame]:
    model = pd.read_csv(TABLES / "fssi_frozen_model.csv")
    zcols = [
        f"frozen_z__{gene}" for gene in model.gene
        if f"frozen_z__{gene}" in data.columns and data[f"frozen_z__{gene}"].notna().all()
    ]
    weights = model.set_index("gene").loc[[column.removeprefix("frozen_z__") for column in zcols], "weight"].to_numpy(float)
    weights = weights / np.abs(weights).sum()
    matrix = data[zcols].to_numpy(float)
    rng = np.random.default_rng(SEED)

    if contrast == "GSE282885":
        positive = data.condition.eq("keloid").to_numpy()
        effect = lambda score: float(score[positive].mean() - score[~positive].mean())
    else:
        hyper = data[data.condition.eq("hypercellular_zone")].set_index("participant_id")
        infiltrating = data[data.condition.eq("infiltrating_zone")].set_index("participant_id")
        order = sorted(set(hyper.index) & set(infiltrating.index))
        hyper_idx = [data.index[data.sample_id.eq(hyper.loc[p, "sample_id"])][0] for p in order]
        infiltrating_idx = [data.index[data.sample_id.eq(infiltrating.loc[p, "sample_id"])][0] for p in order]
        positions = {idx: pos for pos, idx in enumerate(data.index)}
        hpos = np.asarray([positions[idx] for idx in hyper_idx])
        ipos = np.asarray([positions[idx] for idx in infiltrating_idx])
        effect = lambda score: float(np.mean(score[hpos] - score[ipos]))

    observed = effect(matrix @ weights)
    values = np.empty(permutations)
    for i in range(permutations):
        values[i] = effect(matrix @ rng.permutation(weights))
    summary = {
        "dataset_id": contrast,
        "null_type": "same_18_genes_random_weight_to_gene_assignment",
        "permutations": permutations,
        "genes_common": len(zcols),
        "observed_effect": observed,
        "null_mean": float(values.mean()),
        "null_sd": float(values.std(ddof=1)),
        "empirical_one_sided_p": float((1 + np.sum(values >= observed)) / (permutations + 1)),
        "observed_percentile": float(np.mean(values <= observed)),
        "seed": SEED,
    }
    return summary, pd.DataFrame({"dataset_id": contrast, "permutation": np.arange(1, permutations + 1), "null_effect": values})


def main() -> None:
    g282 = pd.read_csv(TABLES / "postlock_GSE282885_frozen_unit_scores.csv")
    g335 = pd.read_csv(TABLES / "postlock_GSE335482_frozen_unit_scores.csv")
    rows: list[dict[str, object]] = []

    for score in SCORES:
        a = g282.loc[g282.condition.eq("keloid"), score].dropna().to_numpy(float)
        b = g282.loc[g282.condition.eq("normal_skin"), score].dropna().to_numpy(float)
        effect, lo, hi, p1, p2 = unpaired_summary(a, b)
        rows.append({
            "dataset_id": "GSE282885", "contrast": "keloid - normal skin",
            "score": score, "n_positive": len(a), "n_reference": len(b),
            "mean_difference": effect, "ci95_low": lo, "ci95_high": hi,
            "exact_one_sided_p": p1, "exact_two_sided_p": p2,
            "direction_support_fraction": np.nan,
            "inference": "library_level_within_dataset_donor_independence_unverified",
        })

    for score in SCORES:
        wide = g335.pivot(index="participant_id", columns="condition", values=score)
        difference = (wide.hypercellular_zone - wide.infiltrating_zone).dropna().to_numpy(float)
        effect, lo, hi, p1, p2 = paired_summary(difference)
        rows.append({
            "dataset_id": "GSE335482", "contrast": "hypercellular - infiltrating zone",
            "score": score, "n_positive": len(difference), "n_reference": len(difference),
            "mean_difference": effect, "ci95_low": lo, "ci95_high": hi,
            "exact_one_sided_p": p1, "exact_two_sided_p": p2,
            "direction_support_fraction": float(np.mean(difference > 0)),
            "inference": "paired_participant_level_within_dataset_gradient_context",
        })

    null_summaries, null_values = [], []
    for dataset, data in [("GSE282885", g282), ("GSE335482", g335)]:
        summary, values = random_weight_null(data.reset_index(drop=True), dataset)
        null_summaries.append(summary)
        null_values.append(values)

    decisions = pd.DataFrame([
        {
            "dataset_id": "GSE282885", "locked_role": "independent disease-control scRNA projection",
            "locked_expectation": "keloid FSSI > normal-skin FSSI",
            "observed_result": "directionally concordant but near-null among five analysable libraries; one normal library yielded no fibroblasts",
            "decision": "retain as a transparent transfer-boundary analysis; do not count as confirmatory donor-level validation",
        },
        {
            "dataset_id": "GSE335482", "locked_role": "paired within-keloid zonal gradient context",
            "locked_expectation": "hypercellular-zone FSSI > infiltrating-zone FSSI",
            "observed_result": "three of four participants concordant; one participant reversed",
            "decision": "retain as paired mechanistic context with small-n uncertainty",
        },
        {
            "dataset_id": "GSE307210", "locked_role": "measured POSTN-knockdown response",
            "locked_expectation": "lower fixed-weight response score and pathological component after siPOSTN",
            "observed_result": "POSTN and Fib_5 module decreased, but the full fixed-weight score and pathological component did not",
            "decision": "retain as target-specific falsification of global state reversal; classify as local-module response",
        },
    ])

    pd.DataFrame(rows).to_csv(TABLES / "postlock_locked_scrna_effects.csv", index=False)
    pd.DataFrame(null_summaries).to_csv(TABLES / "postlock_locked_scrna_random_weight_summary.csv", index=False)
    pd.concat(null_values, ignore_index=True).to_csv(TABLES / "postlock_locked_scrna_random_weight_null.csv.gz", index=False, compression="gzip")
    decisions.to_csv(TABLES / "postlock_locked_validation_decisions.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(pd.DataFrame(null_summaries).to_string(index=False))
    print(decisions.to_string(index=False))


if __name__ == "__main__":
    main()
