from __future__ import annotations

import gzip
import itertools
import json
import math
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)
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


def read_platform_annotation(path: Path) -> pd.DataFrame:
    lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                lines.append(line)
    return pd.read_csv(StringIO("".join(lines)), sep="\t", dtype=str, low_memory=False)


def t_interval(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return math.nan, math.nan
    mean = float(np.mean(values))
    sem = float(stats.sem(values))
    margin = float(stats.t.ppf(0.975, len(values) - 1) * sem)
    return mean - margin, mean + margin


def exact_sign_flip(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    observed = abs(float(np.mean(values)))
    permutations = np.array(
        [abs(float(np.mean(values * signs))) for signs in itertools.product([-1, 1], repeat=len(values))]
    )
    return float(np.mean(permutations >= observed - 1e-12))


def summarise_response(dataset: str, contrast: str, values: pd.Series, reference: str) -> dict[str, object]:
    vector = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if len(vector) == 0:
        return {
            "dataset": dataset,
            "contrast": contrast,
            "reference_tissue": reference,
            "n_units": 0,
            "mean_response": math.nan,
            "ci95_low": math.nan,
            "ci95_high": math.nan,
            "positive_units": 0,
            "exact_signflip_p_two_sided": math.nan,
            "response_sd": math.nan,
        }
    low, high = t_interval(vector)
    return {
        "dataset": dataset,
        "contrast": contrast,
        "reference_tissue": reference,
        "n_units": len(vector),
        "mean_response": float(np.mean(vector)),
        "ci95_low": low,
        "ci95_high": high,
        "positive_units": int(np.sum(vector > 0)),
        "exact_signflip_p_two_sided": exact_sign_flip(vector),
        "response_sd": float(np.std(vector, ddof=1)) if len(vector) > 1 else math.nan,
    }


def fixed_score(expression: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame]:
    available = [gene for gene in GENES if gene in expression.index]
    matrix = expression.loc[available].astype(float)
    gene_sd = matrix.std(axis=1, ddof=1).replace(0, np.nan)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(gene_sd, axis=0).fillna(0)
    weights = WEIGHTS.loc[available]
    score = z.mul(weights, axis=0).sum(axis=0) / weights.abs().sum()
    return score, z


def analyse_gse90051() -> tuple[pd.DataFrame, dict[str, object]]:
    matrix = read_series_matrix(STAGE / "GSE90051_series_matrix.txt.gz")
    annotation = read_platform_annotation(STAGE / "GPL6480.annot.gz")
    annotation = annotation[["ID", "Gene symbol"]].dropna()
    annotation["Gene symbol"] = annotation["Gene symbol"].str.upper().str.split("///").str[0].str.strip()
    merged = matrix.merge(annotation, left_on=matrix.columns[0], right_on="ID", how="left")
    observed = merged[merged["Gene symbol"].isin(GENES)].copy()
    sample_columns = [column for column in matrix.columns[1:]]
    observed[[matrix.columns[0], "Gene symbol"] + sample_columns].to_csv(
        TABLES / "GSE90051_frozen_probe_values.csv", index=False
    )
    gene_ratio = observed.groupby("Gene symbol")[sample_columns].median()
    rank_matrix = matrix.copy()
    for sample in sample_columns:
        values = pd.to_numeric(rank_matrix[sample], errors="coerce")
        valid = values.notna()
        transformed = pd.Series(np.nan, index=values.index, dtype=float)
        ranks = values.loc[valid].rank(method="average").to_numpy(float)
        transformed.loc[valid] = stats.norm.ppf((ranks - 0.375) / (valid.sum() + 0.25))
        rank_matrix[sample] = transformed
    rank_merged = rank_matrix.merge(annotation, left_on=rank_matrix.columns[0], right_on="ID", how="left")
    rank_gene = rank_merged.loc[rank_merged["Gene symbol"].isin(GENES)].groupby("Gene symbol")[sample_columns].median()
    available = [gene for gene in GENES if gene in gene_ratio.index]
    weights = WEIGHTS.loc[available]
    responses = gene_ratio.loc[available].mul(weights, axis=0).sum(axis=0) / weights.abs().sum()
    rank_responses = rank_gene.loc[available].mul(weights, axis=0).sum(axis=0) / weights.abs().sum()
    pathological = MODEL.loc[MODEL.program.eq("pathological"), "gene"].tolist()
    repair = MODEL.loc[MODEL.program.eq("repair"), "gene"].tolist()
    rows = []
    for index, sample in enumerate(sample_columns, 1):
        rows.append(
            {
                "dataset": "GSE90051",
                "patient_id": f"Patient_{index}",
                "gsm": sample,
                "fixed_weight_keloid_minus_adjacent_response": responses[sample],
                "rank_normalised_fixed_response": rank_responses[sample],
                "pathological_component": gene_ratio.loc[[g for g in pathological if g in available], sample].mul(
                    WEIGHTS.loc[[g for g in pathological if g in available]]
                ).sum() / WEIGHTS.loc[[g for g in pathological if g in available]].abs().sum(),
                "repair_component_reversed": gene_ratio.loc[[g for g in repair if g in available], sample].mul(
                    WEIGHTS.loc[[g for g in repair if g in available]]
                ).sum() / WEIGHTS.loc[[g for g in repair if g in available]].abs().sum(),
                "genes_detected": len(available),
                "genes_total": len(GENES),
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(TABLES / "GSE90051_patient_responses.csv", index=False)
    gene_ratio.loc[available].reset_index().to_csv(TABLES / "GSE90051_frozen_gene_ratios.csv", index=False)
    summary = summarise_response(
        "GSE90051",
        "active_keloid_minus_adjacent_normal_skin",
        output.fixed_weight_keloid_minus_adjacent_response,
        "same-patient adjacent normal skin",
    )
    summary["genes_detected"] = len(available)
    rank_summary = summarise_response(
        "GSE90051",
        "active_keloid_minus_adjacent_normal_skin_rank_sensitivity",
        output.rank_normalised_fixed_response,
        "same-patient adjacent normal skin",
    )
    summary.update({f"rank_sensitivity_{key}": value for key, value in rank_summary.items() if key not in {"dataset", "reference_tissue"}})
    summary["source_scale_deviation"] = "GEO labels VALUE as normalised log10 ratio, but released values extend far beyond a plausible log10 range"
    return output, summary


def gse83286_feature_map(path: Path) -> pd.DataFrame:
    header: list[str] | None = None
    rows: list[tuple[str, str]] = []
    in_features = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("FEATURES\t"):
                header = line.rstrip("\n").split("\t")[1:]
                in_features = True
                continue
            if in_features and line.startswith("DATA\t") and header is not None:
                fields = line.rstrip("\n").split("\t")[1:]
                if len(fields) >= len(header):
                    record = dict(zip(header, fields))
                    rows.append((record.get("ProbeName", ""), record.get("GeneName", "")))
            elif in_features and line.startswith("*"):
                break
    mapping = pd.DataFrame(rows, columns=["ID_REF", "gene"])
    mapping["gene"] = mapping.gene.astype(str).str.upper().str.split("///").str[0].str.strip()
    return mapping.drop_duplicates()


def analyse_gse83286() -> tuple[pd.DataFrame, dict[str, object]]:
    matrix = read_series_matrix(STAGE / "GSE83286_series_matrix.txt.gz")
    mapping = gse83286_feature_map(STAGE / "GSM2198196_N1.txt.gz")
    mapping.to_csv(TABLES / "GSE83286_probe_gene_map.csv.gz", index=False, compression="gzip")
    merged = matrix.merge(mapping, on="ID_REF", how="left")
    sample_columns = matrix.columns[1:].tolist()
    expression = merged.loc[merged.gene.isin(GENES)].groupby("gene")[sample_columns].median()
    if expression.empty:
        output = pd.DataFrame(
            columns=["dataset", "patient_id", "normal_gsm", "keloid_gsm", "normal_score", "keloid_score", "keloid_minus_normal_response", "genes_detected", "genes_total"]
        )
        output.to_csv(TABLES / "GSE83286_paired_responses.csv", index=False)
        summary = summarise_response(
            "GSE83286", "earlobe_keloid_minus_normal_skin", pd.Series(dtype=float), "reported paired normal specimen"
        )
        summary["genes_detected"] = 0
        summary["status"] = "not_evaluable_without_public_probe_to_symbol_annotation"
        return output, summary
    score, _ = fixed_score(expression)
    rows = []
    for index in range(3):
        normal = sample_columns[index]
        keloid = sample_columns[index + 3]
        rows.append(
            {
                "dataset": "GSE83286",
                "patient_id": f"Pair_{index + 1}",
                "normal_gsm": normal,
                "keloid_gsm": keloid,
                "normal_score": score[normal],
                "keloid_score": score[keloid],
                "keloid_minus_normal_response": score[keloid] - score[normal],
                "genes_detected": len(expression.index.intersection(GENES)),
                "genes_total": len(GENES),
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(TABLES / "GSE83286_paired_responses.csv", index=False)
    summary = summarise_response(
        "GSE83286",
        "earlobe_keloid_minus_normal_skin",
        output.keloid_minus_normal_response,
        "reported paired normal specimen",
    )
    summary["genes_detected"] = int(output.genes_detected.iloc[0])
    return output, summary


def analyse_gse212954() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frame = pd.read_excel(STAGE / "GSE212954_Expression_Gene.xlsx", sheet_name="gene_exp", header=8)
    frame["Gene_Name"] = frame.Gene_Name.astype(str).str.upper()
    sample_columns = [column for column in frame.columns if str(column).startswith(("C-", "M-", "N-"))]
    expression = frame.groupby("Gene_Name")[sample_columns].median()
    expression = np.log2(expression.clip(lower=0) + 1)
    score, _ = fixed_score(expression)
    records = []
    for sample in sample_columns:
        zone, patient = sample.split("-", 1)
        records.append(
            {
                "dataset": "GSE212954",
                "sample_id": sample,
                "patient_id": patient,
                "zone": {"C": "centre", "M": "margin", "N": "normal"}[zone],
                "fixed_weight_response": score[sample],
                "genes_detected": sum(gene in expression.index for gene in GENES),
                "genes_total": len(GENES),
            }
        )
    samples = pd.DataFrame(records)
    samples.to_csv(TABLES / "GSE212954_sample_scores.csv", index=False)
    wide = samples.pivot(index="patient_id", columns="zone", values="fixed_weight_response")
    contrasts: list[dict[str, object]] = []
    paired_rows: list[dict[str, object]] = []
    for left, right, label in [
        ("centre", "normal", "centre_minus_normal"),
        ("margin", "normal", "margin_minus_normal"),
        ("centre", "margin", "centre_minus_margin"),
    ]:
        paired = wide[[left, right]].dropna()
        response = paired[left] - paired[right]
        for patient, value in response.items():
            paired_rows.append({"dataset": "GSE212954", "patient_id": patient, "contrast": label, "response": value})
        contrasts.append(
            summarise_response(
                "GSE212954",
                label,
                response,
                "same-patient normal skin" if right == "normal" else "same-patient keloid margin",
            )
        )
    pd.DataFrame(paired_rows).to_csv(TABLES / "GSE212954_paired_responses.csv", index=False)
    return samples, contrasts


def analyse_gse282479() -> tuple[pd.DataFrame, list[dict[str, object]]]:
    counts = pd.read_csv(STAGE / "GSE282479_VitD_counts.csv.gz", index_col=0)
    counts.index = counts.index.astype(str).str.split(".").str[0]
    mapping = pd.read_csv(STAGE.parent / "external" / "GSE191067" / "all_ensembl_symbol_map.csv")
    mapping.ENSEMBL = mapping.ENSEMBL.astype(str).str.split(".").str[0]
    merged = counts.reset_index(names="ENSEMBL").merge(mapping, on="ENSEMBL", how="left")
    merged["SYMBOL"] = merged.SYMBOL.astype(str).str.upper()
    sample_columns = counts.columns.tolist()
    expression = merged.groupby("SYMBOL")[sample_columns].sum()
    library_sizes = expression.sum(axis=0)
    log_cpm = np.log2(expression.div(library_sizes, axis=1) * 1_000_000 + 1)
    score, _ = fixed_score(log_cpm)
    sample_map = []
    for prefix, context, donors in [("K", "keloid", 3), ("ST", "normal", 4)]:
        for donor in range(1, donors + 1):
            untreated = f"{prefix}{2 * donor - 1}_sorted.bam"
            treated = f"{prefix}{2 * donor}_sorted.bam"
            sample_map.extend(
                [
                    (untreated, context, str(donor), "untreated"),
                    (treated, context, str(donor), "paricalcitol"),
                ]
            )
    records = []
    for sample, context, donor, treatment in sample_map:
        records.append(
            {
                "dataset": "GSE282479",
                "sample_id": sample,
                "context": context,
                "donor_id": f"{context}_{donor}",
                "treatment": treatment,
                "fixed_weight_response": score[sample],
                "RXRA_log2CPM": float(log_cpm.loc["RXRA", sample]),
                "VDR_log2CPM": float(log_cpm.loc["VDR", sample]),
                "genes_detected": sum(gene in expression.index for gene in GENES),
                "genes_total": len(GENES),
            }
        )
    samples = pd.DataFrame(records)
    samples.to_csv(TABLES / "GSE282479_paricalcitol_sample_scores.csv", index=False)
    responses = []
    summaries = []
    for context, subset in samples.groupby("context"):
        wide = subset.pivot(index="donor_id", columns="treatment", values="fixed_weight_response")
        reduction = wide.untreated - wide.paricalcitol
        for donor, value in reduction.items():
            responses.append({"dataset": "GSE282479", "context": context, "donor_id": donor, "untreated_minus_paricalcitol": value})
        summaries.append(
            summarise_response(
                "GSE282479",
                f"{context}_untreated_minus_paricalcitol",
                reduction,
                "same-donor untreated fibroblasts",
            )
        )
    pd.DataFrame(responses).to_csv(TABLES / "GSE282479_paricalcitol_paired_responses.csv", index=False)
    return samples, summaries


def hedges_one_sample(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    n = len(values)
    d = float(np.mean(values) / np.std(values, ddof=1))
    correction = float(special.gamma((n - 1) / 2) / (math.sqrt((n - 1) / 2) * special.gamma((n - 2) / 2))) if n > 2 else 1.0
    g = correction * d
    variance = correction**2 * (1 / n + d**2 / (2 * (n - 1)))
    return g, variance


def random_effects_dl(effects: pd.DataFrame) -> dict[str, float]:
    y = effects.hedges_g.to_numpy(float)
    variance = effects.variance.to_numpy(float)
    fixed_weight = 1 / variance
    fixed_mean = float(np.sum(fixed_weight * y) / np.sum(fixed_weight))
    q = float(np.sum(fixed_weight * (y - fixed_mean) ** 2))
    df = len(y) - 1
    c = float(np.sum(fixed_weight) - np.sum(fixed_weight**2) / np.sum(fixed_weight))
    tau2 = max(0.0, (q - df) / c)
    weight = 1 / (variance + tau2)
    estimate = float(np.sum(weight * y) / np.sum(weight))
    se = float(math.sqrt(1 / np.sum(weight)))
    return {
        "method": "DerSimonian-Laird random effects",
        "studies": len(y),
        "standardised_mean_response": estimate,
        "ci95_low": estimate - 1.96 * se,
        "ci95_high": estimate + 1.96 * se,
        "tau2": tau2,
        "Q": q,
        "Q_df": df,
        "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 else 0.0,
    }


def build_paired_tissue_summary(gse90051: pd.DataFrame, gse212954: pd.DataFrame, gse83286: pd.DataFrame) -> dict[str, object]:
    responses: list[tuple[str, str, np.ndarray]] = []
    responses.append(
        (
            "GSE90051",
            "active keloid vs adjacent normal skin",
            gse90051.fixed_weight_keloid_minus_adjacent_response.to_numpy(float),
        )
    )
    if not gse83286.empty:
        responses.append(
            (
                "GSE83286",
                "earlobe keloid vs reported paired normal skin",
                gse83286.keloid_minus_normal_response.to_numpy(float),
            )
        )
    paired_212 = pd.read_csv(TABLES / "GSE212954_paired_responses.csv")
    responses.append(
        (
            "GSE212954",
            "keloid centre vs same-patient normal skin",
            paired_212.loc[paired_212.contrast.eq("centre_minus_normal"), "response"].to_numpy(float),
        )
    )
    existing = pd.read_csv(TABLES / "GSE158395_locked_sample_scores.csv")
    lesion = existing.loc[existing.condition.eq("lesion")].set_index("participant_id").fixed_weight_response
    nonlesion = existing.loc[existing.condition.eq("nonlesion")].set_index("participant_id").fixed_weight_response
    responses.append(
        (
            "GSE158395",
            "keloid lesion vs same-patient distant non-lesional skin",
            (lesion - nonlesion).dropna().to_numpy(float),
        )
    )
    rows = []
    for dataset, contrast, vector in responses:
        g, variance = hedges_one_sample(vector)
        rows.append(
            {
                "dataset": dataset,
                "contrast": contrast,
                "n_patients": len(vector),
                "positive_patients": int(np.sum(vector > 0)),
                "mean_dataset_scale_response": float(np.mean(vector)),
                "hedges_g": g,
                "variance": variance,
                "se": math.sqrt(variance),
                "ci95_low": g - 1.96 * math.sqrt(variance),
                "ci95_high": g + 1.96 * math.sqrt(variance),
                "exact_signflip_p_two_sided": exact_sign_flip(vector),
            }
        )
    effects = pd.DataFrame(rows)
    effects.to_csv(TABLES / "patient_paired_tissue_standardised_effects.csv", index=False)
    meta = random_effects_dl(effects)
    pd.DataFrame([meta]).to_csv(TABLES / "patient_paired_tissue_random_effects.csv", index=False)
    return meta


def main() -> None:
    gse90051, summary_90051 = analyse_gse90051()
    gse212954, summaries_212954 = analyse_gse212954()
    gse83286, summary_83286 = analyse_gse83286()
    _, summaries_282479 = analyse_gse282479()
    meta = build_paired_tissue_summary(gse90051, gse212954, gse83286)
    all_summaries = pd.DataFrame([summary_90051, summary_83286] + summaries_212954 + summaries_282479)
    all_summaries.to_csv(TABLES / "new_transcriptomic_contrast_summary.csv", index=False)
    payload = {
        "GSE90051": summary_90051,
        "GSE212954": summaries_212954,
        "GSE83286": summary_83286,
        "GSE282479": summaries_282479,
        "paired_tissue_random_effects": meta,
    }
    (TABLES / "transcriptomic_expansion_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
