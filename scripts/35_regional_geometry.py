from __future__ import annotations

import gzip
import itertools
import json
import math
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion"
TABLES = ROOT / "results" / "tables"
MODEL = pd.read_csv(TABLES / "fssi_frozen_model.csv")
MODEL["gene"] = MODEL.gene.str.upper()
GENES = MODEL.gene.tolist()
WEIGHTS = MODEL.set_index("gene").weight


def read_series_matrix(path: Path) -> pd.DataFrame:
    lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                lines.append(line)
    frame = pd.read_csv(StringIO("".join(lines)), sep="\t", quotechar='"')
    frame.columns = [str(column).strip('"') for column in frame.columns]
    frame.iloc[:, 0] = frame.iloc[:, 0].astype(str).str.strip('"')
    return frame


def fixed_score(expression: pd.DataFrame) -> pd.Series:
    available = [gene for gene in GENES if gene in expression.index]
    matrix = expression.loc[available].astype(float)
    gene_sd = matrix.std(axis=1, ddof=1).replace(0, np.nan)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(gene_sd, axis=0).fillna(0)
    weights = WEIGHTS.loc[available]
    return z.mul(weights, axis=0).sum(axis=0) / weights.abs().sum()


def t_interval(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return math.nan, math.nan
    margin = float(stats.t.ppf(0.975, len(values) - 1) * stats.sem(values))
    return float(np.mean(values) - margin), float(np.mean(values) + margin)


def exact_sign_flip(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    observed = abs(float(np.mean(values)))
    null = [abs(float(np.mean(values * signs))) for signs in itertools.product([-1, 1], repeat=len(values))]
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def paired_summary(contrast: str, values: pd.Series) -> dict[str, object]:
    vector = values.dropna().to_numpy(float)
    low, high = t_interval(vector)
    return {
        "dataset": "GSE218007",
        "contrast": contrast,
        "n_patients": len(vector),
        "mean_within_patient_difference": float(np.mean(vector)),
        "ci95_low": low,
        "ci95_high": high,
        "positive_patients": int(np.sum(vector > 0)),
        "exact_signflip_p_two_sided": exact_sign_flip(vector),
    }


def welch_summary(left: pd.Series, right: pd.Series) -> dict[str, float | int]:
    left_values = left.dropna().to_numpy(float)
    right_values = right.dropna().to_numpy(float)
    estimate = float(np.mean(left_values) - np.mean(right_values))
    left_var = float(np.var(left_values, ddof=1))
    right_var = float(np.var(right_values, ddof=1))
    se2 = left_var / len(left_values) + right_var / len(right_values)
    df = se2**2 / (
        (left_var / len(left_values)) ** 2 / (len(left_values) - 1)
        + (right_var / len(right_values)) ** 2 / (len(right_values) - 1)
    )
    margin = float(stats.t.ppf(0.975, df) * math.sqrt(se2))
    test = stats.ttest_ind(left_values, right_values, equal_var=False)
    return {
        "n_keloid_patients": len(left_values),
        "n_normal_donors": len(right_values),
        "mean_keloid_minus_normal": estimate,
        "ci95_low": estimate - margin,
        "ci95_high": estimate + margin,
        "welch_df": df,
        "welch_p_two_sided": float(test.pvalue),
    }


def sample_metadata() -> pd.DataFrame:
    samples = [f"GSM{value}" for value in range(6732515, 6732544)]
    keloid = [
        ("SOU22", "nodular", "CP"), ("SOU22", "nodular", "CR"), ("SOU22", "nodular", "PP"), ("SOU22", "nodular", "PR"),
        ("CAM31", "nodular", "CP"), ("CAM31", "nodular", "CR"), ("CAM31", "nodular", "PP"), ("CAM31", "nodular", "PR"),
        ("KAR25", "nodular", "CP"), ("KAR25", "nodular", "CR"), ("KAR25", "nodular", "PP"), ("KAR25", "nodular", "PR"),
        ("DUA21", "extensive", "CP"), ("DUA21", "extensive", "CR"), ("DUA21", "extensive", "PP"), ("DUA21", "extensive", "PR"),
        ("OSS69", "extensive", "CP"), ("OSS69", "extensive", "CR"), ("OSS69", "extensive", "PP"),
        ("NKO27", "extensive", "CP"), ("NKO27", "extensive", "CR"), ("NKO27", "extensive", "PP"), ("NKO27", "extensive", "PR"),
    ]
    normal = [
        ("MAJ23", "normal", "PAP"), ("MAJ23", "normal", "RET"),
        ("MOU19", "normal", "PAP"), ("MOU19", "normal", "RET"),
        ("OUL25", "normal", "PAP"), ("OUL25", "normal", "RET"),
    ]
    rows = []
    for sample, (donor, morphology, region) in zip(samples, keloid + normal):
        rows.append({
            "sample_id": sample,
            "donor_id": donor,
            "disease": "keloid" if morphology != "normal" else "normal",
            "morphology": morphology,
            "region": region,
        })
    return pd.DataFrame(rows)


def main() -> None:
    matrix = read_series_matrix(STAGE / "GSE218007_series_matrix.txt.gz")
    mapping = pd.read_csv(TABLES / "GPL23126_frozen_probe_symbol_map.csv")
    mapping["SYMBOL"] = mapping.SYMBOL.str.upper()
    observed = matrix.merge(mapping[["PROBEID", "SYMBOL"]], left_on="ID_REF", right_on="PROBEID", how="inner")
    sample_columns = matrix.columns[1:].tolist()
    expression = observed.groupby("SYMBOL")[sample_columns].median()
    score = fixed_score(expression)

    metadata = sample_metadata()
    metadata["frozen_score"] = metadata.sample_id.map(score)
    metadata["genes_detected"] = len(expression.index.intersection(GENES))
    metadata["genes_total"] = len(GENES)
    metadata.to_csv(TABLES / "GSE218007_regional_sample_scores.csv", index=False)

    donor = metadata.groupby(["donor_id", "disease", "morphology"], as_index=False).frozen_score.mean()
    donor.to_csv(TABLES / "GSE218007_donor_scores.csv", index=False)
    disease_summary = {
        "dataset": "GSE218007",
        "contrast": "keloid_patient_mean_minus_normal_donor_mean",
        **welch_summary(
            donor.loc[donor.disease.eq("keloid"), "frozen_score"],
            donor.loc[donor.disease.eq("normal"), "frozen_score"],
        ),
    }

    keloid = metadata.loc[metadata.disease.eq("keloid")].copy()
    keloid["axis_centre_periphery"] = keloid.region.str[0].map({"C": "centre", "P": "periphery"})
    keloid["axis_depth"] = keloid.region.str[1].map({"P": "superficial", "R": "deep"})
    centre = keloid.groupby(["donor_id", "axis_centre_periphery"]).frozen_score.mean().unstack()
    depth = keloid.groupby(["donor_id", "axis_depth"]).frozen_score.mean().unstack()
    centre_difference = centre.centre - centre.periphery
    depth_difference = depth.deep - depth.superficial
    regional = pd.concat(
        [
            centre_difference.rename("centre_minus_periphery"),
            depth_difference.rename("deep_minus_superficial"),
        ],
        axis=1,
    ).reset_index()
    regional.to_csv(TABLES / "GSE218007_within_patient_regional_responses.csv", index=False)

    summaries = [
        disease_summary,
        paired_summary("centre_minus_periphery", centre_difference),
        paired_summary("deep_minus_superficial", depth_difference),
    ]
    summary = {
        "source": "GSE218007",
        "biological_unit": "six keloid patients and three normal-skin donors",
        "regions": "CP/CR/PP/PR; OSS69 has no released PR sample and no value was imputed",
        "genes_detected": len(expression.index.intersection(GENES)),
        "genes_total": len(GENES),
        "summaries": summaries,
    }
    (TABLES / "GSE218007_regional_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
