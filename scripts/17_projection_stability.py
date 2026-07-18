"""Summarise final release sample effects, fixed-cell downsampling and model stability."""

from __future__ import annotations

import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "tables"
SEED = 20260712
ITERATIONS = 1000


def exact_unpaired(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    pooled = np.concatenate([a, b])
    observed = float(a.mean() - b.mean())
    null = []
    for idx in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), bool); mask[list(idx)] = True
        null.append(pooled[mask].mean() - pooled[~mask].mean())
    null = np.asarray(null)
    return float(np.mean(null >= observed - 1e-12)), float(np.mean(np.abs(null) >= abs(observed) - 1e-12))


def welch(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float, float]:
    effect = float(a.mean() - b.mean())
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se2 = va / len(a) + vb / len(b)
    df = se2**2 / ((va / len(a))**2 / (len(a) - 1) + (vb / len(b))**2 / (len(b) - 1))
    q = stats.t.ppf(.975, df)
    p1, p2 = exact_unpaired(a, b)
    return effect, effect - q * math.sqrt(se2), effect + q * math.sqrt(se2), p1, p2


def paired(d: np.ndarray) -> tuple[float, float, float, float, float]:
    effect = float(d.mean()); se = stats.sem(d); q = stats.t.ppf(.975, len(d) - 1)
    null = np.asarray([np.mean(d * np.asarray(s)) for s in itertools.product([-1, 1], repeat=len(d))])
    return effect, effect - q * se, effect + q * se, float(np.mean(null >= effect - 1e-12)), float(np.mean(np.abs(null) >= abs(effect) - 1e-12))


def summarise_g191() -> pd.DataFrame:
    units = pd.read_csv(T / "GSE191067_frozen_sample_scores_sensitivity.csv")
    contrasts = [
        ("keloid", "normal_scar", "keloid - normal scar", True),
        ("keloid", "normal_skin", "keloid - normal skin", False),
        ("keloid", "perilesional_skin", "keloid - perilesional skin", False),
        ("normal_scar", "normal_skin", "normal scar - normal skin", False),
    ]
    rows = []
    for (threshold, annotation), frame in units.groupby(["threshold", "annotation_sensitivity"]):
        for score in ["fssi_frozen_weighted", "fssi_frozen_equal_weight", "pathological_component",
                      "repair_component_reversed", "Fib5_three_marker_baseline", "collagen_ECM_raw"]:
            for positive, reference, label, primary in contrasts:
                a = frame.loc[frame.condition.eq(positive), score].dropna().to_numpy(float)
                b = frame.loc[frame.condition.eq(reference), score].dropna().to_numpy(float)
                if min(len(a), len(b)) < 2:
                    continue
                effect, lo, hi, p1, p2 = welch(a, b)
                rows.append({"dataset_id": "GSE191067", "contrast": label, "score": score,
                             "threshold": threshold, "annotation_sensitivity": annotation,
                             "n_positive": len(a), "n_reference": len(b), "mean_difference": effect,
                             "ci95_low": lo, "ci95_high": hi, "exact_one_sided_p": p1,
                             "exact_two_sided_p": p2, "locked_primary": primary and threshold == .5 and annotation == "reference_transfer"})
    return pd.DataFrame(rows)


def downsample_dataset(path: Path, dataset: str, contrasts: list[tuple[str, str, str]], paired_design: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    use = ["sample_id", "condition", "fssi_frozen_weighted"]
    data = pd.read_csv(path, usecols=lambda c: c in use + ["donor_id", "participant_id", "prediction_score"])
    if dataset == "GSE191067":
        data = data[data.prediction_score >= .5].copy()
        data["unit_id"] = data.sample_id
    elif "participant_id" in data:
        data["unit_id"] = data.participant_id
    elif "donor_id" in data:
        data["unit_id"] = data.donor_id
    else:
        data["unit_id"] = data.sample_id
    per_sample_n = data.groupby("sample_id").size()
    fixed_n = int(per_sample_n.min())
    rng = np.random.default_rng(SEED)
    iteration_rows = []
    for iteration in range(1, ITERATIONS + 1):
        sampled = data.groupby("sample_id", group_keys=False).sample(n=fixed_n, replace=False, random_state=int(rng.integers(1, 2**31 - 1)))
        sample_means = sampled.groupby(["sample_id", "unit_id", "condition"], as_index=False).fssi_frozen_weighted.mean()
        unit_means = sample_means.groupby(["unit_id", "condition"], as_index=False).fssi_frozen_weighted.mean()
        for positive, reference, label in contrasts:
            if paired_design:
                wide = unit_means.pivot(index="unit_id", columns="condition", values="fssi_frozen_weighted")
                effect = float((wide[positive] - wide[reference]).dropna().mean())
            else:
                a = unit_means.loc[unit_means.condition.eq(positive), "fssi_frozen_weighted"]
                b = unit_means.loc[unit_means.condition.eq(reference), "fssi_frozen_weighted"]
                effect = float(a.mean() - b.mean())
            iteration_rows.append({"dataset_id": dataset, "contrast": label, "iteration": iteration,
                                   "fixed_cells_per_library": fixed_n, "mean_difference": effect})
    values = pd.DataFrame(iteration_rows)
    summary = values.groupby(["dataset_id", "contrast", "fixed_cells_per_library"], as_index=False).agg(
        downsample_median_effect=("mean_difference", "median"),
        monte_carlo_95_low=("mean_difference", lambda x: np.quantile(x, .025)),
        monte_carlo_95_high=("mean_difference", lambda x: np.quantile(x, .975)),
        prespecified_direction_fraction=("mean_difference", lambda x: np.mean(x > 0)),
    )
    summary["iterations"] = ITERATIONS
    return summary, values


def model_stability() -> tuple[pd.DataFrame, pd.DataFrame]:
    effects = pd.read_csv(T / "fssi_training_scenario_gene_effects.csv")
    original = pd.read_csv(T / "fssi_frozen_model.csv")
    original_genes = set(original.gene)
    rows, selections = [], []
    datasets = sorted(effects.dataset_id.unique())
    for heldout in datasets:
        sub = effects[effects.dataset_id.ne(heldout)].copy()
        stability = sub.groupby("gene", as_index=False).agg(
            min_directional_strength=("directional_strength", "min"),
            mean_standardized_difference=("standardized_difference", "mean"),
            scenarios=("dataset_id", "size"))
        program = original.set_index("gene").program.to_dict()
        # Candidate universes are encoded by expected direction in the locked training table.
        expected = effects.drop_duplicates("gene").set_index("gene").expected_direction
        stability["program"] = stability.gene.map(lambda g: "pathological" if expected[g] > 0 else "repair")
        eligible = stability[stability.min_directional_strength > 0]
        chosen = pd.concat([eligible[eligible.program.eq("pathological")].nlargest(10, "min_directional_strength"),
                            eligible[eligible.program.eq("repair")].nlargest(8, "min_directional_strength")])
        selected = set(chosen.gene)
        rows.append({"heldout_dataset": heldout, "training_datasets_remaining": ";".join(sorted(set(sub.dataset_id))),
                     "selected_genes_n": len(selected), "original_genes_retained": len(selected & original_genes),
                     "jaccard_with_frozen_model": len(selected & original_genes) / len(selected | original_genes),
                     "both_programmes_retained": chosen.program.nunique() == 2})
        for _, row in stability.iterrows():
            selections.append({"heldout_dataset": heldout, "gene": row.gene, "program": row.program,
                               "selected": row.gene in selected, "mean_standardized_difference": row.mean_standardized_difference,
                               "weight_sign_expected": (row.mean_standardized_difference > 0) == (row.program == "pathological")})
    detail = pd.DataFrame(selections)
    gene = detail.groupby(["gene", "program"], as_index=False).agg(
        lodo_selection_fraction=("selected", "mean"), weight_sign_stability=("weight_sign_expected", "mean"))
    return pd.DataFrame(rows), gene


def communication_summary() -> pd.DataFrame:
    data = pd.read_csv(T / "GSE191067_target_blind_communication_by_sample.csv")
    rows = []
    for axis, sub in data.groupby("axis_id"):
        support = int(sub.expected_direction_support.sum()); n = len(sub)
        rows.append({"dataset_id": "GSE191067", "axis_id": axis, "samples_evaluable": n,
                     "support_n": support, "support_fraction": support / n,
                     "exact_one_sided_sign_p": stats.binomtest(support, n, .5, alternative="greater").pvalue,
                     "median_oriented_delta": float(np.median(np.where(sub.expected_delta.eq("positive"), sub.delta_high_minus_low, -sub.delta_high_minus_low)))})
    out = pd.DataFrame(rows)
    out["fdr_bh"] = stats.false_discovery_control(out.exact_one_sided_sign_p)
    return out


def main() -> None:
    effects = summarise_g191()
    effects.to_csv(T / "GSE191067_sample_level_effects.csv", index=False)
    configs = [
        (T / "external_scrna_GSE156326_frozen_cell_scores.csv.gz", "GSE156326", [("hypertrophic_scar", "normal_skin", "hypertrophic scar - normal skin")], False),
        (T / "external_scrna_GSE181316_frozen_cell_scores.csv.gz", "GSE181316", [("keloid", "normal_scar", "keloid - normal scar")], False),
        (T / "postlock_GSE282885_frozen_cell_scores.csv.gz", "GSE282885", [("keloid", "normal_skin", "keloid - normal skin")], False),
        (T / "postlock_GSE335482_frozen_cell_scores.csv.gz", "GSE335482", [("hypercellular_zone", "infiltrating_zone", "hypercellular - infiltrating zone")], True),
        (T / "GSE191067_frozen_fibroblast_cell_scores.csv.gz", "GSE191067", [("keloid", "normal_scar", "keloid - normal scar"), ("keloid", "normal_skin", "keloid - normal skin")], False),
    ]
    summaries, iterations = [], []
    for path, dataset, contrasts, is_paired in configs:
        summary, values = downsample_dataset(path, dataset, contrasts, is_paired)
        summaries.append(summary); iterations.append(values)
    pd.concat(summaries, ignore_index=True).to_csv(T / "fixed_cell_downsampling_summary.csv", index=False)
    pd.concat(iterations, ignore_index=True).to_csv(T / "fixed_cell_downsampling_iterations.csv.gz", index=False, compression="gzip")
    folds, genes = model_stability()
    folds.to_csv(T / "fssi_leave_one_development_dataset_out_summary.csv", index=False)
    genes.to_csv(T / "fssi_lodo_gene_stability.csv", index=False)
    comm = communication_summary()
    comm.to_csv(T / "GSE191067_target_blind_communication_summary.csv", index=False)
    print(effects.query("locked_primary").to_string(index=False))
    print(pd.concat(summaries, ignore_index=True).to_string(index=False))
    print(folds.to_string(index=False)); print(comm.to_string(index=False))


if __name__ == "__main__":
    main()
