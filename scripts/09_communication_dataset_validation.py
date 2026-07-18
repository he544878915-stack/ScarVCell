"""Dataset-stratified and leave-one-dataset-out communication validation."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
SOURCE = TABLES / "communication_target_blind_sample_axis_scores.csv"


def donor_ids(frame: pd.DataFrame, scenario: str, collapsed_pair: tuple[str, str] | None) -> pd.Series:
    ids = frame["dataset_id"].astype(str) + ":" + frame["sample_id"].astype(str)
    ids.loc[frame.dataset_id.eq("GSE243716")] = "GSE243716:P1"
    if collapsed_pair is not None:
        mask = frame.dataset_id.eq("GSE163973") & frame.sample_id.isin(collapsed_pair)
        ids.loc[mask] = "GSE163973:" + scenario
    return ids


def summarize(values: pd.Series) -> dict[str, float | int]:
    values = values.dropna().astype(float)
    support = int((values > 0).sum())
    n = int(len(values))
    return {
        "donor_n": n,
        "support_n": support,
        "support_fraction": support / n if n else np.nan,
        "exact_one_sided_sign_p": binomtest(support, n, .5, alternative="greater").pvalue if n else np.nan,
        "oriented_mean": values.mean() if n else np.nan,
        "oriented_median": values.median() if n else np.nan,
    }


def main() -> None:
    sample = pd.read_csv(SOURCE)
    sample["oriented_delta"] = sample.delta_high_minus_low * np.where(
        sample.expected_delta.eq("positive"), 1.0, -1.0
    )
    scenarios: dict[str, tuple[str, str] | None] = {"all_matrices_distinct": None}
    for pair in combinations(["NF1_matrix", "NF2_matrix", "NF3_matrix"], 2):
        scenarios["collapse_" + "_".join(item.replace("_matrix", "") for item in pair)] = pair

    dataset_rows: list[dict] = []
    donor_rows: list[pd.DataFrame] = []
    lodo_rows: list[dict] = []
    concordance_rows: list[dict] = []
    for scenario, collapsed_pair in scenarios.items():
        current = sample.copy()
        current["scenario"] = scenario
        current["donor_id"] = donor_ids(current, scenario, collapsed_pair)
        donor = current.groupby(
            ["scenario", "axis_id", "ligand", "receptor", "expected_delta", "dataset_id", "donor_id"],
            as_index=False,
        ).oriented_delta.mean()
        donor_rows.append(donor)

        for keys, subset in donor.groupby(["axis_id", "ligand", "receptor", "expected_delta", "dataset_id"]):
            axis, ligand, receptor, expected, dataset = keys
            dataset_rows.append({"scenario": scenario, "axis_id": axis, "ligand": ligand,
                                 "receptor": receptor, "expected_delta": expected, "dataset_id": dataset,
                                 **summarize(subset.oriented_delta)})

        dataset_means = donor.groupby(["axis_id", "dataset_id"], as_index=False).oriented_delta.mean()
        for axis, subset in dataset_means.groupby("axis_id"):
            values = subset.oriented_delta
            concordance_rows.append({"scenario": scenario, "axis_id": axis,
                                     "dataset_n": len(values), "dataset_support_n": int((values > 0).sum()),
                                     "dataset_support_fraction": float((values > 0).mean()),
                                     "exact_one_sided_study_sign_p": binomtest(int((values > 0).sum()), len(values), .5,
                                                                               alternative="greater").pvalue,
                                     "minimum_dataset_oriented_mean": values.min(),
                                     "maximum_dataset_oriented_mean": values.max()})
            for omitted in sorted(subset.dataset_id.unique()):
                kept_datasets = set(subset.dataset_id) - {omitted}
                kept = donor[(donor.axis_id.eq(axis)) & donor.dataset_id.isin(kept_datasets)]
                lodo_rows.append({"scenario": scenario, "axis_id": axis, "omitted_dataset": omitted,
                                  "remaining_dataset_n": len(kept_datasets), **summarize(kept.oriented_delta)})

    donor_all = pd.concat(donor_rows, ignore_index=True)
    dataset_all = pd.DataFrame(dataset_rows)
    lodo_all = pd.DataFrame(lodo_rows)
    concordance_all = pd.DataFrame(concordance_rows)

    dataset_robust = dataset_all.groupby(["axis_id", "ligand", "receptor", "expected_delta", "dataset_id"], as_index=False).agg(
        donor_n_min=("donor_n", "min"), donor_n_max=("donor_n", "max"),
        support_fraction_min=("support_fraction", "min"), support_fraction_max=("support_fraction", "max"),
        worst_exact_one_sided_sign_p=("exact_one_sided_sign_p", "max"),
        oriented_mean_min=("oriented_mean", "min"), oriented_mean_max=("oriented_mean", "max"),
    )
    concordance_robust = concordance_all.groupby("axis_id", as_index=False).agg(
        dataset_n_min=("dataset_n", "min"), dataset_support_fraction_min=("dataset_support_fraction", "min"),
        worst_exact_one_sided_study_sign_p=("exact_one_sided_study_sign_p", "max"),
        minimum_dataset_oriented_mean=("minimum_dataset_oriented_mean", "min"),
        maximum_dataset_oriented_mean=("maximum_dataset_oriented_mean", "max"),
    )
    concordance_robust["study_level_fdr"] = concordance_robust.worst_exact_one_sided_study_sign_p.rank(method="max")
    order = concordance_robust.worst_exact_one_sided_study_sign_p.sort_values().index
    p = concordance_robust.loc[order, "worst_exact_one_sided_study_sign_p"].to_numpy()
    adjusted = np.minimum.accumulate((p * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    concordance_robust.loc[order, "study_level_fdr"] = np.minimum(adjusted, 1)

    donor_all.to_csv(TABLES / "communication_axis_dataset_stratified_donor_values.csv", index=False)
    dataset_all.to_csv(TABLES / "communication_axis_dataset_stratified_scenarios.csv", index=False)
    dataset_robust.to_csv(TABLES / "communication_axis_dataset_stratified_robust.csv", index=False)
    lodo_all.to_csv(TABLES / "communication_axis_leave_one_dataset_out.csv", index=False)
    concordance_all.to_csv(TABLES / "communication_axis_dataset_concordance_scenarios.csv", index=False)
    concordance_robust.to_csv(TABLES / "communication_axis_dataset_concordance_robust.csv", index=False)
    print("Wrote dataset-stratified communication validation tables")


if __name__ == "__main__":
    main()
