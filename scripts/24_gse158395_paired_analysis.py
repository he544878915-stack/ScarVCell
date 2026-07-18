"""Locked final release analysis of the African-American GSE158395 skin cohort.

The FSSI genes and weights are never refitted. Duplicate lesion libraries from
JKR782-003 are averaged before experiment-wide standardisation. The primary
contrast is the three participant-paired lesional versus non-lesional effects;
lesional versus healthy skin and composition adjustment are sensitivities.
"""

from __future__ import annotations

import gzip
import itertools
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls


ROOT = Path(os.environ.get("SCARVCELL_ROOT", Path(__file__).resolve().parents[1]))
INPUT = ROOT / ".stage" / "genetics" / "GSE158395" / "GSE158395_exp_fin_geo.csv.gz"
TABLES = ROOT / "results" / "tables"
MODEL_FILE = TABLES / "fssi_frozen_model.csv"
SOFT_FILE = ROOT / ".stage" / "genetics" / "GSE158395" / "GSE158395_family.soft.gz"


def score_matrix(expr: pd.DataFrame, model: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = [gene for gene in model.gene if gene in expr.index]
    x = expr.loc[genes].T
    z = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1).replace(0, 1).fillna(1)
    m = model.set_index("gene").loc[genes]
    weight = m.weight.to_numpy(float).copy()
    weight /= np.abs(weight).sum()
    equal = m.equal_weight.to_numpy(float).copy()
    equal /= np.abs(equal).sum()
    pathological = [gene for gene in genes if m.loc[gene, "program"] == "pathological"]
    repair = [gene for gene in genes if m.loc[gene, "program"] == "repair"]

    score = pd.DataFrame(index=x.index)
    score["fixed_weight_response"] = z.to_numpy() @ weight
    score["equal_weight_response"] = z.to_numpy() @ equal
    score["pathological_component"] = z[pathological].mean(axis=1)
    score["repair_component_reversed"] = -z[repair].mean(axis=1)
    score["genes_present"] = len(genes)
    score["genes_total"] = len(model)
    score.index.name = "sample_id"

    contributions = z.mul(pd.Series(weight, index=genes), axis=1)
    contributions.index.name = "sample_id"
    return score.reset_index(), contributions


def paired_ci(values: np.ndarray) -> tuple[float, float, float, float]:
    values = np.asarray(values, float)
    effect = float(values.mean())
    if len(values) < 2 or np.std(values, ddof=1) == 0:
        return effect, effect, effect, 1.0
    se = np.std(values, ddof=1) / math.sqrt(len(values))
    q = stats.t.ppf(0.975, len(values) - 1)
    p = 2 * stats.t.sf(abs(effect / se), len(values) - 1)
    return effect, effect - q * se, effect + q * se, float(p)


def exact_signflip(values: np.ndarray) -> float:
    values = np.asarray(values, float)
    observed = abs(values.mean())
    null = [abs(np.mean(values * np.asarray(signs)))
            for signs in itertools.product([-1, 1], repeat=len(values))]
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def welch_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    effect = float(a.mean() - b.mean())
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se2 = va / len(a) + vb / len(b)
    if se2 <= 0:
        return effect, effect, effect, 1.0
    df = se2**2 / ((va / len(a))**2 / (len(a) - 1) + (vb / len(b))**2 / (len(b) - 1))
    q = stats.t.ppf(0.975, df)
    p = 2 * stats.t.sf(abs(effect / math.sqrt(se2)), df)
    return effect, effect - q * math.sqrt(se2), effect + q * math.sqrt(se2), float(p)


def exact_independent_permutation(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    pooled = np.concatenate([a, b])
    observed = abs(a.mean() - b.mean())
    null = []
    for selected in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(selected)] = True
        null.append(abs(pooled[mask].mean() - pooled[~mask].mean()))
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def deconvolve(expr: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    donor = pd.read_csv(TABLES / "scrna_broad_celltype_donor_logcpm_reference.csv", index_col=0)
    meta = pd.read_csv(TABLES / "scrna_broad_celltype_donor_reference_metadata.csv")
    markers = pd.read_csv(TABLES / "scrna_broad_celltype_specificity_markers.csv")
    selected = markers[(markers.rank_within_type <= 50) & ~markers.gene.isin(set(model.gene))]
    genes = sorted(set(selected.gene) & set(donor.index) & set(expr.index))
    cell_types = sorted(meta.broad_cell_type.unique())
    signature = np.column_stack([
        donor.loc[genes, meta.loc[meta.broad_cell_type.eq(cell_type), "reference_unit"]]
        .mean(axis=1).to_numpy()
        for cell_type in cell_types
    ])
    signature = np.column_stack([
        stats.rankdata(signature[:, index]) / len(genes) for index in range(len(cell_types))
    ])
    rows = []
    for sample_id in expr.columns:
        vector = stats.rankdata(expr.loc[genes, sample_id].to_numpy()) / len(genes)
        coefficient, residual = nnls(signature, vector)
        coefficient = coefficient / coefficient.sum() if coefficient.sum() else coefficient
        rows.extend({
            "sample_id": sample_id,
            "broad_cell_type": cell_type,
            "proportion": value,
            "fit_rmse": residual / math.sqrt(len(genes)),
            "marker_genes": len(genes),
        } for cell_type, value in zip(cell_types, coefficient))
    return pd.DataFrame(rows)


def metadata() -> pd.DataFrame:
    rows = [
        ("JKR782-001_lesion", "JKR782-001", "lesion", 1, "GSM4798871"),
        ("JKR782-001_nonlesion", "JKR782-001", "nonlesion", 1, "GSM4798879"),
        ("JKR782-003_lesion", "JKR782-003", "lesion", 2, "GSM4798868;GSM4798869"),
        ("JKR782-003_nonlesion", "JKR782-003", "nonlesion", 1, "GSM4798873"),
        ("JKR782-004_lesion", "JKR782-004", "lesion", 1, "GSM4798878"),
        ("JKR782-004_nonlesion", "JKR782-004", "nonlesion", 1, "GSM4798870"),
        ("JKR491-163_healthy", "JKR491-163", "healthy", 1, "GSM4798880"),
        ("N1_healthy", "N1", "healthy", 1, "GSM4798877"),
        ("N2_healthy", "N2", "healthy", 1, "GSM4798872"),
        ("N3_healthy", "N3", "healthy", 1, "GSM4798875"),
        ("N19_healthy", "N19", "healthy", 1, "GSM4798876"),
        ("N20_healthy", "N20", "healthy", 1, "GSM4798874"),
    ]
    return pd.DataFrame(rows, columns=["sample_id", "participant_id", "condition",
                                       "source_profiles", "source_gsm"])


def geo_soft_metadata() -> pd.DataFrame:
    records, current = [], None
    with gzip.open(SOFT_FILE, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["sample_title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, content = value.split(": ", 1)
                    current[key] = content
        if current is not None:
            records.append(current)
    audit = pd.DataFrame(records)
    audit["paper_aligned_healthy_control"] = audit.gsm.isin(
        ["GSM4798872", "GSM4798874", "GSM4798875", "GSM4798876", "GSM4798877"])
    audit["additional_geo_normal_matrix"] = audit.gsm.eq("GSM4798880")
    audit["same_participant_multiple_lesional_profiles"] = (
        audit.ptid.eq("JKR782-003") & audit.lesion.eq("LS"))
    audit["emerging_lesion_profile_mapping"] = np.where(
        audit.same_participant_multiple_lesional_profiles, "not_resolved_between_K4_and_K5", "not_applicable")
    return audit


def build_expression() -> pd.DataFrame:
    raw = pd.read_csv(INPUT, compression="gzip", index_col=0)
    raw.index = raw.index.astype(str).str.upper()
    raw = raw.apply(pd.to_numeric, errors="coerce")
    raw = raw.groupby(level=0).median()
    columns = {
        "JKR782-001_lesion": ["JKR782-001_Keloid_LS_K2"],
        "JKR782-001_nonlesion": ["JKR782-001_Keloid_NL_K1"],
        "JKR782-003_lesion": ["JKR782-003_Keloid_LS_K5", "JKR782-003_Keloid_LS_K4"],
        "JKR782-003_nonlesion": ["JKR782-003_Keloid_NL_K3"],
        "JKR782-004_lesion": ["JKR782-004_Keloid_LS_K7"],
        "JKR782-004_nonlesion": ["JKR782-004_Keloid_NL_K6"],
        "JKR491-163_healthy": ["JKR491-163/BA_Normal_Normal_410"],
        "N1_healthy": ["N1_Normal_Normal_N1"],
        "N2_healthy": ["N2_Normal_Normal_N2"],
        "N3_healthy": ["N3_Normal_Normal_N3"],
        "N19_healthy": ["N19_Normal_Normal_N19"],
        "N20_healthy": ["N20_Normal_Normal_N20"],
    }
    missing = sorted({column for values in columns.values() for column in values} - set(raw.columns))
    if missing:
        raise ValueError(f"Missing locked source columns: {missing}")
    return pd.DataFrame({unit: raw[source].mean(axis=1) for unit, source in columns.items()})


def contrast_tables(scores: pd.DataFrame, endpoints: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    paired_rows, independent_rows = [], []
    paired = scores[scores.condition.isin(["lesion", "nonlesion"])]
    for endpoint in endpoints:
        wide = paired.pivot(index="participant_id", columns="condition", values=endpoint).dropna()
        differences = wide.lesion - wide.nonlesion
        effect, low, high, p_t = paired_ci(differences.to_numpy())
        paired_rows.append({
            "dataset_id": "GSE158395", "endpoint": endpoint,
            "contrast": "lesion_minus_matched_nonlesion", "analysis_role": "primary",
            "participant_pairs": len(differences), "mean_difference": effect,
            "ci95_low": low, "ci95_high": high, "paired_t_p": p_t,
            "exact_signflip_p_two_sided": exact_signflip(differences.to_numpy()),
            "positive_pairs": int((differences > 0).sum()),
        })
        lesion = scores.loc[scores.condition.eq("lesion"), endpoint].to_numpy(float)
        healthy = scores.loc[scores.condition.eq("healthy"), endpoint].to_numpy(float)
        effect, low, high, p_t = welch_ci(lesion, healthy)
        independent_rows.append({
            "dataset_id": "GSE158395", "endpoint": endpoint,
            "contrast": "lesion_minus_healthy", "analysis_role": "secondary_reference_sensitivity",
            "n_lesion": len(lesion), "n_healthy": len(healthy), "mean_difference": effect,
            "ci95_low": low, "ci95_high": high, "welch_p": p_t,
            "exact_permutation_p_two_sided": exact_independent_permutation(lesion, healthy),
        })
        publication_healthy = scores.loc[
            scores.condition.eq("healthy") & ~scores.participant_id.eq("JKR491-163"), endpoint
        ].to_numpy(float)
        effect, low, high, p_t = welch_ci(lesion, publication_healthy)
        independent_rows.append({
            "dataset_id": "GSE158395", "endpoint": endpoint,
            "contrast": "lesion_minus_publication_aligned_healthy",
            "analysis_role": "secondary_reference_metadata_sensitivity",
            "n_lesion": len(lesion), "n_healthy": len(publication_healthy),
            "mean_difference": effect, "ci95_low": low, "ci95_high": high,
            "welch_p": p_t,
            "exact_permutation_p_two_sided": exact_independent_permutation(lesion, publication_healthy),
        })
    return pd.DataFrame(paired_rows), pd.DataFrame(independent_rows)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    model = pd.read_csv(MODEL_FILE)
    geo_soft_metadata().to_csv(TABLES / "GSE158395_GEO_sample_metadata_audit.csv", index=False)
    expr = build_expression()
    meta = metadata()
    scores, contributions = score_matrix(expr, model)
    scores = scores.merge(meta, on="sample_id", validate="one_to_one")
    scores.to_csv(TABLES / "GSE158395_locked_sample_scores.csv", index=False)

    endpoints = ["fixed_weight_response", "equal_weight_response",
                 "pathological_component", "repair_component_reversed"]
    paired, independent = contrast_tables(scores, endpoints)
    paired.to_csv(TABLES / "GSE158395_paired_effects.csv", index=False)
    independent.to_csv(TABLES / "GSE158395_healthy_reference_effects.csv", index=False)

    contribution_rows = []
    contribution = contributions.reset_index().merge(meta, on="sample_id")
    for gene in model.gene:
        wide = contribution[contribution.condition.isin(["lesion", "nonlesion"])].pivot(
            index="participant_id", columns="condition", values=gene).dropna()
        differences = wide.lesion - wide.nonlesion
        contribution_rows.append({
            "gene": gene, "participant_pairs": len(differences),
            "mean_weighted_contribution_difference": differences.mean(),
            "positive_pairs": int((differences > 0).sum()),
            "exact_signflip_p_two_sided": exact_signflip(differences.to_numpy()),
        })
    pd.DataFrame(contribution_rows).to_csv(
        TABLES / "GSE158395_gene_contribution_differences.csv", index=False)

    proportions = deconvolve(expr, model)
    proportions.to_csv(TABLES / "GSE158395_rank_nnls_deconvolution.csv", index=False)
    fibroblast = proportions.loc[proportions.broad_cell_type.eq("fibroblast"),
                                 ["sample_id", "proportion"]].rename(
                                     columns={"proportion": "fibroblast_proportion"})
    adjusted = scores.merge(fibroblast, on="sample_id", validate="one_to_one")
    design = np.column_stack([np.ones(len(adjusted)), adjusted.fibroblast_proportion])
    for endpoint in endpoints:
        fitted = design @ np.linalg.lstsq(design, adjusted[endpoint], rcond=None)[0]
        adjusted[f"{endpoint}_composition_residual"] = adjusted[endpoint] - fitted
    adjusted.to_csv(TABLES / "GSE158395_composition_sensitivity_scores.csv", index=False)
    residual_endpoints = [f"{endpoint}_composition_residual" for endpoint in endpoints]
    residual_paired, residual_independent = contrast_tables(adjusted, residual_endpoints)
    residual_paired["analysis_role"] = "composition_sensitivity_primary_contrast"
    residual_independent["analysis_role"] = "composition_sensitivity_secondary_contrast"
    pd.concat([residual_paired, residual_independent], ignore_index=True, sort=False).to_csv(
        TABLES / "GSE158395_composition_sensitivity_effects.csv", index=False)

    audit = pd.DataFrame([{
        "dataset_id": "GSE158395", "public_profiles": 13,
        "analysis_units_after_duplicate_collapse": len(expr.columns),
        "keloid_participant_pairs": 3, "healthy_matrices_in_geo": 6,
        "healthy_controls_reported_in_source_paper": 5,
        "same_participant_lesional_profiles_averaged": 2,
        "additional_normal_matrix_not_counted_in_source_paper": "GSM4798880",
        "expression_genes_after_symbol_collapse": len(expr),
        "fssi_genes_present": int(scores.genes_present.min()),
        "fssi_genes_total": len(model),
        "primary_inference": "participant-paired lesion minus non-lesion",
        "secondary_inference": "lesion minus independent healthy skin",
        "composition_method": "Rank-NNLS; top-50 broad-cell markers excluding FSSI genes",
    }])
    audit.to_csv(TABLES / "GSE158395_analysis_audit.csv", index=False)
    print(paired.to_string(index=False))
    print(independent.to_string(index=False))


if __name__ == "__main__":
    main()
