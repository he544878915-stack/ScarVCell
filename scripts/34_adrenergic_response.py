from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

WORKBOOK = STAGE / "41467_2026_72823_MOESM5_ESM.xlsx"
MODEL = pd.read_csv(TABLES / "fssi_frozen_model.csv")
MODEL["gene"] = MODEL.gene.str.upper()
GENES = MODEL.gene.tolist()
WEIGHTS = MODEL.set_index("gene").weight
MECHANISTIC_GENES = ["RUNX2", "IBSP", "ADRB1", "TSN", "CREB1"]

GROUPS = {
    "control": ["Ctr1", "Ctr2", "Ctr3"],
    "epinephrine": ["AE1", "AE2", "AE3"],
    "epinephrine_metoprolol": ["Meto1", "Meto2", "Meto3"],
}


def fixed_score(expression: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    available = [gene for gene in GENES if gene in expression.index]
    matrix = expression.loc[available].astype(float)
    gene_sd = matrix.std(axis=1, ddof=1).replace(0, np.nan)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(gene_sd, axis=0).fillna(0)
    weights = WEIGHTS.loc[available]
    score = z.mul(weights, axis=0).sum(axis=0) / weights.abs().sum()
    return score, z


def welch_difference(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    estimate = float(np.mean(left) - np.mean(right))
    left_var = float(np.var(left, ddof=1))
    right_var = float(np.var(right, ddof=1))
    se2 = left_var / len(left) + right_var / len(right)
    se = math.sqrt(se2)
    numerator = se2**2
    denominator = (left_var / len(left)) ** 2 / (len(left) - 1) + (right_var / len(right)) ** 2 / (len(right) - 1)
    df = numerator / denominator if denominator > 0 else math.nan
    margin = float(stats.t.ppf(0.975, df) * se) if np.isfinite(df) else math.nan
    test = stats.ttest_ind(left, right, equal_var=False)
    return {
        "estimate": estimate,
        "ci95_low": estimate - margin,
        "ci95_high": estimate + margin,
        "welch_df": df,
        "welch_p_two_sided": float(test.pvalue),
        "left_mean": float(np.mean(left)),
        "right_mean": float(np.mean(right)),
        "left_sd": float(np.std(left, ddof=1)),
        "right_sd": float(np.std(right, ddof=1)),
    }


def indexed_difference(left: np.ndarray, right: np.ndarray) -> dict[str, float | int]:
    """Descriptive sensitivity only; replicate suffixes are not asserted to be paired."""
    values = np.asarray(left, dtype=float) - np.asarray(right, dtype=float)
    return {
        "indexed_mean_difference": float(np.mean(values)),
        "indexed_positive_differences": int(np.sum(values > 0)),
        "indexed_total_differences": int(len(values)),
    }


def read_timepoint(sheet: str) -> pd.DataFrame:
    frame = pd.read_excel(WORKBOOK, sheet_name=sheet)
    frame["Symbol"] = frame.Symbol.astype(str).str.upper().str.strip()
    sample_columns = sum(GROUPS.values(), [])
    expression = frame.groupby("Symbol")[sample_columns].mean()
    return np.log2(expression.clip(lower=0) + 1)


def analyse_timepoint(sheet: str, timepoint: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = read_timepoint(sheet)
    score, z = fixed_score(expression)
    sample_rows: list[dict[str, object]] = []
    for group, samples in GROUPS.items():
        for sample in samples:
            row: dict[str, object] = {
                "source": "Lou et al., Nature Communications, 2026",
                "timepoint": timepoint,
                "sample_id": sample,
                "group": group,
                "frozen_score": float(score[sample]),
                "genes_detected": int(sum(gene in expression.index for gene in GENES)),
                "genes_total": len(GENES),
            }
            for gene in MECHANISTIC_GENES:
                row[f"{gene}_log2_abundance"] = float(expression.loc[gene, sample]) if gene in expression.index else math.nan
            sample_rows.append(row)

    contrast_rows: list[dict[str, object]] = []
    contrasts = [
        ("epinephrine", "control", "epinephrine_minus_control", "state induction"),
        (
            "epinephrine",
            "epinephrine_metoprolol",
            "epinephrine_minus_epinephrine_metoprolol",
            "metoprolol-associated reversal",
        ),
    ]
    for left_group, right_group, contrast, interpretation in contrasts:
        left = score.loc[GROUPS[left_group]].to_numpy(float)
        right = score.loc[GROUPS[right_group]].to_numpy(float)
        result: dict[str, object] = {
            "source": "Lou et al., Nature Communications, 2026",
            "timepoint": timepoint,
            "contrast": contrast,
            "interpretation": interpretation,
            "n_left": len(left),
            "n_right": len(right),
            "independent_unit": "cultured-fibroblast replicate",
            "pairing_assumed_for_primary_inference": False,
        }
        result.update(welch_difference(left, right))
        result.update(indexed_difference(left, right))
        contrast_rows.append(result)

    contribution_rows: list[dict[str, object]] = []
    for gene in [gene for gene in GENES if gene in z.index]:
        weight = float(WEIGHTS.loc[gene])
        for left_group, right_group, contrast, _ in contrasts:
            left_mean = float(z.loc[gene, GROUPS[left_group]].mean())
            right_mean = float(z.loc[gene, GROUPS[right_group]].mean())
            contribution_rows.append(
                {
                    "timepoint": timepoint,
                    "contrast": contrast,
                    "gene": gene,
                    "program": MODEL.set_index("gene").loc[gene, "program"],
                    "frozen_weight": weight,
                    "mean_z_difference": left_mean - right_mean,
                    "weighted_contribution": (left_mean - right_mean) * weight / WEIGHTS.loc[z.index].abs().sum(),
                }
            )
    return pd.DataFrame(sample_rows), pd.DataFrame(contrast_rows), pd.DataFrame(contribution_rows)


def main() -> None:
    sample_frames = []
    contrast_frames = []
    contribution_frames = []
    for sheet, timepoint in [("6 h_RNA-seq", "6 h"), ("24 h_RNA-seq", "24 h")]:
        samples, contrasts, contributions = analyse_timepoint(sheet, timepoint)
        sample_frames.append(samples)
        contrast_frames.append(contrasts)
        contribution_frames.append(contributions)

    samples = pd.concat(sample_frames, ignore_index=True)
    contrasts = pd.concat(contrast_frames, ignore_index=True)
    contributions = pd.concat(contribution_frames, ignore_index=True)
    samples.to_csv(TABLES / "Lou2026_adrenergic_sample_scores.csv", index=False)
    contrasts.to_csv(TABLES / "Lou2026_adrenergic_contrasts.csv", index=False)
    contributions.to_csv(TABLES / "Lou2026_adrenergic_gene_contributions.csv", index=False)

    summary = {
        "source": "Lou et al., Nature Communications, 2026; DOI: 10.1038/s41467-026-72823-9",
        "source_file": WORKBOOK.name,
        "source_group_mapping": {
            "Ctr": "control",
            "AE": "epinephrine",
            "Meto": "epinephrine and metoprolol",
        },
        "primary_inference": "Welch comparison of three cultured-fibroblast replicates per group",
        "indexed_difference": "descriptive sensitivity only; replicate suffixes were not treated as biological pairing",
        "score_definition": "frozen 18-gene weights applied after log2(x+1) and within-timepoint gene standardisation",
        "contrasts": contrasts.to_dict(orient="records"),
    }
    (TABLES / "Lou2026_adrenergic_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    print(contrasts.to_string(index=False))


if __name__ == "__main__":
    main()
