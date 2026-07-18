from __future__ import annotations

import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "public_expansion"
OUT = ROOT / "results" / "tables"
MODEL_PATH = OUT / "fssi_frozen_model.csv"

# Reviewed human canonical accessions, frozen before inspecting either proteome.
UNIPROT = {
    "ADAM12": "O43184",
    "COMP": "P49747",
    "COL5A1": "P20908",
    "POSTN": "Q15063",
    "FAP": "Q12884",
    "ASPN": "Q9BXN1",
    "COL12A1": "Q99715",
    "COL16A1": "Q07092",
    "OGN": "P20774",
    "C1QTNF3": "Q9BXJ4",
    "PI16": "Q6UXB8",
    "DCN": "P07585",
    "APOD": "P05090",
    "ABCA8": "O94911",
    "C3": "P01024",
    "ABCA10": "Q8WWZ4",
    "GSTM5": "P46439",
    "NEGR1": "Q7Z3B1",
}

# Specimen classes transcribed from Supplementary Table 1 of the source article.
KELOID = {"SC167_4", "SC136_3", "SC141_7"}
NORMAL = {"NOR202_1", "NOR208_4", "NOR207_2", "NOR197_3", "NOR58_12", "NOR164_9", "NOR164_16", "NOR205_2"}
HYPERTROPHIC = {
    "SC211_9", "SC117_7", "SC116_7", "SC144_5", "SC168_4", "SC154_2", "SC158_3", "SC55_4",
    "SC70_6", "SC68_3", "SC215_4", "SC216_6", "SC171_2", "SC66_6", "SC203_3",
}


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=pvalues.index, dtype=float)
    valid = pvalues.dropna().sort_values()
    if valid.empty:
        return result
    n = len(valid)
    adjusted = np.minimum.accumulate((valid.to_numpy() * n / np.arange(1, n + 1))[::-1])[::-1]
    result.loc[valid.index] = np.minimum(adjusted, 1.0)
    return result


def welch_effect(case: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    case = np.asarray(case, dtype=float)
    reference = np.asarray(reference, dtype=float)
    case, reference = case[np.isfinite(case)], reference[np.isfinite(reference)]
    out = {"n_case": len(case), "n_reference": len(reference)}
    if len(case) < 2 or len(reference) < 2:
        return out | {"mean_case": np.nan, "mean_reference": np.nan, "difference": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan}
    difference = float(case.mean() - reference.mean())
    variance = case.var(ddof=1) / len(case) + reference.var(ddof=1) / len(reference)
    se = math.sqrt(variance)
    numerator = variance**2
    denominator = ((case.var(ddof=1) / len(case)) ** 2 / (len(case) - 1)) + ((reference.var(ddof=1) / len(reference)) ** 2 / (len(reference) - 1))
    df = numerator / denominator if denominator > 0 else np.inf
    critical = stats.t.ppf(0.975, df) if se > 0 else 0.0
    p_value = stats.ttest_ind(case, reference, equal_var=False).pvalue
    return out | {
        "mean_case": float(case.mean()), "mean_reference": float(reference.mean()), "difference": difference,
        "ci_low": difference - critical * se, "ci_high": difference + critical * se, "p_value": float(p_value),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sample_id(column: str) -> str:
    match = re.search(r"XB04182B1DA_(NOR\d+_\d+|SC\d+_\d+)_", column)
    if not match:
        raise ValueError(f"Cannot parse specimen ID from {column}")
    return match.group(1)


def analyse_pxd066833(model: pd.DataFrame) -> dict:
    path = RAW / "PXD066833" / "extracted" / "report.pg_matrix.tsv"
    matrix = pd.read_csv(path, sep="\t")
    quantity_columns = list(matrix.columns[5:])
    rename = {column: sample_id(column) for column in quantity_columns}
    observed = set(rename.values())
    expected = KELOID | NORMAL | HYPERTROPHIC
    if observed != expected:
        raise ValueError(f"PXD066833 specimen mismatch: missing={expected-observed}, extra={observed-expected}")

    manifest_rows = []
    for specimen in sorted(observed):
        group = "keloid" if specimen in KELOID else "hypertrophic_scar" if specimen in HYPERTROPHIC else "normal_skin"
        donor_guardrail = "possible same-donor repeat; collapsed in sensitivity" if specimen in {"NOR164_9", "NOR164_16"} else "donor identity not asserted"
        manifest_rows.append({"specimen_id": specimen, "group": group, "statistical_unit": "specimen", "donor_guardrail": donor_guardrail})
    pd.DataFrame(manifest_rows).to_csv(OUT / "PXD066833_sample_manifest.csv", index=False)

    selected_rows = []
    for record in model.itertuples(index=False):
        accession = UNIPROT[record.gene]
        hits = matrix[matrix["Protein.Ids"].fillna("").str.split(";").apply(lambda values: accession in values)].copy()
        if hits.empty:
            selected_rows.append({"gene": record.gene, "accession": accession, "selection": "not_detected"})
            continue
        exact = hits[hits["Protein.Group"].astype(str) == accession]
        candidate = exact if not exact.empty else hits
        medians = candidate[quantity_columns].apply(pd.to_numeric, errors="coerce").median(axis=1, skipna=True)
        row = candidate.loc[medians.idxmax()]
        values = pd.to_numeric(row[quantity_columns], errors="coerce")
        selected = {"gene": record.gene, "accession": accession, "selection": "canonical_exact" if not exact.empty else "canonical_containing_group", "protein_group": row["Protein.Group"]}
        selected.update({rename[column]: values[column] for column in quantity_columns})
        selected_rows.append(selected)

    selected = pd.DataFrame(selected_rows)
    selected.to_csv(OUT / "PXD066833_selected_protein_groups.csv", index=False)
    detected = selected[selected["selection"] != "not_detected"].copy()
    specimen_columns = sorted(observed)
    abundance = detected.set_index("gene")[specimen_columns].apply(pd.to_numeric, errors="coerce")
    log_abundance = np.log2(abundance.where(abundance > 0))

    effects = []
    comparisons = [("keloid_vs_normal", KELOID, NORMAL), ("hypertrophic_scar_vs_normal", HYPERTROPHIC, NORMAL)]
    for comparison, cases, refs in comparisons:
        for gene in model.gene:
            if gene not in log_abundance.index:
                effect = welch_effect(np.array([]), np.array([]))
                detection = 0
            else:
                effect = welch_effect(log_abundance.loc[gene, sorted(cases)].to_numpy(), log_abundance.loc[gene, sorted(refs)].to_numpy())
                detection = int(log_abundance.loc[gene].notna().sum())
            effects.append({"dataset": "PXD066833", "comparison": comparison, "gene": gene, "program": model.set_index("gene").loc[gene, "program"], "detected_specimens": detection, **effect})
    effects = pd.DataFrame(effects)
    effects["fdr_within_comparison"] = effects.groupby("comparison", group_keys=False)["p_value"].apply(bh_adjust)
    effects.to_csv(OUT / "PXD066833_fixed_protein_effects.csv", index=False)

    # A protein-domain projection is deliberately partial and renormalises only frozen weights with adequate coverage.
    coverage = log_abundance.notna().mean(axis=1)
    # Complete-case proteins keep every specimen in the projection and avoid
    # data-dependent missing-value imputation or specimen-specific weight sets.
    eligible = [gene for gene in model.gene if gene in coverage.index and coverage[gene] == 1.0]
    n_path = int(model.set_index("gene").loc[eligible, "program"].eq("pathological").sum()) if eligible else 0
    n_repair = int(model.set_index("gene").loc[eligible, "program"].eq("repair").sum()) if eligible else 0
    score_status = "computed" if n_path >= 4 and n_repair >= 3 else "not_evaluable"
    score_table = pd.DataFrame(manifest_rows).set_index("specimen_id")
    if score_status == "computed":
        values = log_abundance.loc[eligible].T
        standardised = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
        weights = model.set_index("gene").loc[eligible, "weight"]
        weights = weights / weights.abs().sum()
        score_table["partial_protein_projection"] = standardised.mul(weights, axis=1).sum(axis=1, min_count=len(eligible))
    else:
        score_table["partial_protein_projection"] = np.nan
    score_table["score_status"] = score_status
    score_table["eligible_fixed_proteins"] = len(eligible)
    score_table.reset_index().to_csv(OUT / "PXD066833_partial_protein_scores.csv", index=False)

    score_effects = []
    if score_status == "computed":
        for comparison, cases, refs in comparisons:
            score_effects.append({"dataset": "PXD066833", "comparison": comparison, **welch_effect(score_table.loc[sorted(cases), "partial_protein_projection"].to_numpy(), score_table.loc[sorted(refs), "partial_protein_projection"].to_numpy())})
        collapsed = score_table.copy()
        nor164 = collapsed.loc[["NOR164_9", "NOR164_16"], "partial_protein_projection"].mean()
        collapsed = collapsed.drop(index=["NOR164_9", "NOR164_16"])
        collapsed.loc["NOR164_collapsed", ["group", "partial_protein_projection"]] = ["normal_skin", nor164]
        ref_values = collapsed.loc[collapsed.group == "normal_skin", "partial_protein_projection"].to_numpy()
        for label, cases in [("keloid_vs_normal_NOR164_collapsed", KELOID), ("hypertrophic_scar_vs_normal_NOR164_collapsed", HYPERTROPHIC)]:
            score_effects.append({"dataset": "PXD066833", "comparison": label, **welch_effect(collapsed.loc[sorted(cases), "partial_protein_projection"].to_numpy(), ref_values)})
    pd.DataFrame(score_effects).to_csv(OUT / "PXD066833_partial_protein_score_effects.csv", index=False)
    return {"matrix_sha256": sha256(path), "detected_fixed_proteins": len(detected), "eligible_fixed_proteins": len(eligible), "eligible_pathological": n_path, "eligible_repair": n_repair, "partial_score_status": score_status}


def parse_mztab(path: Path, target_accessions: set[str]) -> tuple[dict[str, int], dict[str, int], int]:
    sequence_accessions: dict[str, set[str]] = defaultdict(set)
    accepted_records: list[tuple[str, str, str]] = []
    header = None
    total_spectra: set[str] = set()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("PSH\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if not line.startswith("PSM\t") or header is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                continue
            row = dict(zip(header, fields))
            if row.get("opt_global_pass_threshold") != "true" or row.get("opt_global_cv_MS:1002217_decoy_peptide") == "1":
                continue
            accession_field = row.get("accession", "")
            match = re.search(r"(?:sp|tr)\|([^|]+)\|", accession_field)
            accession = match.group(1).split("-")[0] if match else accession_field.split("-")[0]
            if accession not in target_accessions:
                continue
            sequence = row.get("sequence", "")
            spectrum = row.get("spectra_ref", "")
            sequence_accessions[sequence].add(accession)
            accepted_records.append((sequence, accession, spectrum))
            total_spectra.add(spectrum)
    unique_psm: dict[str, set[str]] = defaultdict(set)
    unique_peptides: dict[str, set[str]] = defaultdict(set)
    for sequence, accession, spectrum in accepted_records:
        if len(sequence_accessions[sequence]) == 1:
            unique_psm[accession].add(spectrum)
            unique_peptides[accession].add(sequence)
    return ({key: len(value) for key, value in unique_psm.items()}, {key: len(value) for key, value in unique_peptides.items()}, len(total_spectra))


def analyse_pxd039729(model: pd.DataFrame) -> dict:
    directory = RAW / "PXD039729"
    path = directory / "S0022202X24019729_mmc1.xlsx"
    if not path.exists():
        pd.DataFrame([{"dataset": "PXD039729", "status": "not_evaluable_missing_quantitative_supplement"}]).to_csv(OUT / "PXD039729_status.csv", index=False)
        return {"status": "not_evaluable_missing_quantitative_supplement"}
    source = pd.read_excel(path, sheet_name="Supplementary Table S2e", header=2)
    clinical = pd.read_excel(path, sheet_name="Supplementary Table S1", header=2)
    case_samples = ["KD24", "KD36", "KD37", "KD33", "KD28"]
    reference_samples = ["NS4", "NS6", "NS7", "NS9", "NS10"]
    expected = set(case_samples + reference_samples)
    available_participants = set(clinical["Sample keys"].dropna().astype(str))
    if not expected.issubset(available_participants):
        raise ValueError(f"PXD039729 participant mapping incomplete: {expected - available_participants}")

    selected_rows = []
    for record in model.itertuples(index=False):
        accession = UNIPROT[record.gene]
        hits = source[source["Accession number(s)"].fillna("").astype(str).str.contains(fr"(?:sp|tr)\|{re.escape(accession)}(?:[;-]|$)", regex=True)].copy()
        if hits.empty:
            selected_rows.append({"gene": record.gene, "accession": accession, "selection": "not_detected"})
            continue
        medians = hits[case_samples + reference_samples].apply(pd.to_numeric, errors="coerce").median(axis=1, skipna=True)
        row = hits.loc[medians.idxmax()]
        selected = {"gene": record.gene, "accession": accession, "selection": "canonical_containing_group", "reported_gene_symbol": row["Gene symbol"], "reported_accession_group": row["Accession number(s)"]}
        selected.update({sample: pd.to_numeric(row[sample], errors="coerce") for sample in case_samples + reference_samples})
        selected_rows.append(selected)
    selected = pd.DataFrame(selected_rows)
    selected.to_csv(OUT / "PXD039729_selected_WCL_proteins.csv", index=False)
    abundance = selected[selected.selection != "not_detected"].set_index("gene")[case_samples + reference_samples].apply(pd.to_numeric, errors="coerce")
    log_abundance = np.log2(abundance.fillna(0).clip(lower=0) + 1)
    long = log_abundance.reset_index().melt(id_vars="gene", var_name="sample", value_name="log2_abundance_plus1")
    long["group"] = np.where(long["sample"].str.startswith("KD"), "keloid_fibroblast", "normal_fibroblast")
    long.to_csv(OUT / "PXD039729_fixed_protein_WCL_abundance.csv", index=False)
    effects = []
    for gene in model.gene:
        if gene in log_abundance.index:
            effect = welch_effect(log_abundance.loc[gene, case_samples].to_numpy(), log_abundance.loc[gene, reference_samples].to_numpy())
            detected_case = int((abundance.loc[gene, case_samples] > 0).sum())
            detected_reference = int((abundance.loc[gene, reference_samples] > 0).sum())
        else:
            effect = welch_effect(np.array([]), np.array([]))
            detected_case = detected_reference = 0
        effects.append({"dataset": "PXD039729", "comparison": "keloid_vs_normal_fibroblast_WCL", "gene": gene, "program": model.set_index("gene").loc[gene, "program"], "detected_case": detected_case, "detected_reference": detected_reference, **effect})
    effects = pd.DataFrame(effects)
    effects["fdr_within_comparison"] = bh_adjust(effects.p_value)
    effects.to_csv(OUT / "PXD039729_fixed_protein_WCL_effects.csv", index=False)

    complete = [gene for gene in model.gene if gene in abundance.index and (abundance.loc[gene] > 0).all()]
    n_path = int(model.set_index("gene").loc[complete, "program"].eq("pathological").sum()) if complete else 0
    n_repair = int(model.set_index("gene").loc[complete, "program"].eq("repair").sum()) if complete else 0
    score_status = "computed" if n_path >= 4 and n_repair >= 3 else "not_evaluable"
    if score_status == "computed":
        values = log_abundance.loc[complete].T
        standardised = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
        weights = model.set_index("gene").loc[complete, "weight"]
        weights = weights / weights.abs().sum()
        scores = standardised.mul(weights, axis=1).sum(axis=1)
        score_table = pd.DataFrame({"sample": scores.index, "group": ["keloid_fibroblast" if sample.startswith("KD") else "normal_fibroblast" for sample in scores.index], "partial_protein_projection": scores.values})
        score_effect = welch_effect(score_table.loc[score_table.group == "keloid_fibroblast", "partial_protein_projection"].to_numpy(), score_table.loc[score_table.group == "normal_fibroblast", "partial_protein_projection"].to_numpy())
    else:
        score_table = pd.DataFrame(columns=["sample", "group", "partial_protein_projection"])
        score_effect = welch_effect(np.array([]), np.array([]))
    score_table["score_status"] = score_status
    score_table["eligible_fixed_proteins"] = len(complete)
    score_table.to_csv(OUT / "PXD039729_partial_protein_scores.csv", index=False)
    pd.DataFrame([{"dataset": "PXD039729", "comparison": "keloid_vs_normal_fibroblast_WCL", **score_effect}]).to_csv(OUT / "PXD039729_partial_protein_score_effect.csv", index=False)
    return {"status": "computed_sample_level_WCL_quantification", "supplement_sha256": sha256(path), "detected_fixed_proteins": len(abundance), "complete_case_fixed_proteins": len(complete), "eligible_pathological": n_path, "eligible_repair": n_repair, "partial_score_status": score_status}


def analyse_gse290035(model: pd.DataFrame) -> dict:
    path = RAW / "GSE290035" / "GSE290035_HDF_ZNF469_count.txt.gz"
    with gzip.open(path, "rt") as handle:
        data = pd.read_csv(handle, sep="\t")
    data = data.drop_duplicates("gene_name").set_index("gene_name")
    control = [f"HDF_shControl_count_{i}" for i in range(1, 4)]
    knockdown = [f"HDF_shZNF469_count_{i}" for i in range(1, 4)]
    counts = data[control + knockdown].apply(pd.to_numeric, errors="coerce")
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)
    rows = []
    for gene in [*model.gene, "ZNF469"]:
        if gene not in log_cpm.index:
            rows.append({"gene": gene, "status": "not_detected"})
            continue
        effect = welch_effect(log_cpm.loc[gene, knockdown].to_numpy(), log_cpm.loc[gene, control].to_numpy())
        rows.append({"gene": gene, "status": "detected", "program": model.set_index("gene").loc[gene, "program"] if gene in set(model.gene) else "perturbed_gene", **effect})
    gene_effects = pd.DataFrame(rows)
    gene_effects["fdr_within_fixed_panel_plus_znf469"] = bh_adjust(gene_effects.p_value)
    gene_effects.to_csv(OUT / "GSE290035_gene_response.csv", index=False)

    eligible = [gene for gene in model.gene if gene in log_cpm.index]
    values = log_cpm.loc[eligible].T
    standardised = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
    weights = model.set_index("gene").loc[eligible, "weight"]
    weights = weights / weights.abs().sum()
    scores = standardised.mul(weights, axis=1).sum(axis=1)
    sample_table = pd.DataFrame({"sample": scores.index, "condition": ["shControl" if "Control" in value else "shZNF469" for value in scores.index], "frozen_state_response_score": scores.values, "statistical_unit": "technical_RNA_library"})
    sample_table.to_csv(OUT / "GSE290035_frozen_response_scores.csv", index=False)
    score_effect = welch_effect(sample_table.loc[sample_table.condition == "shZNF469", "frozen_state_response_score"].to_numpy(), sample_table.loc[sample_table.condition == "shControl", "frozen_state_response_score"].to_numpy())
    pd.DataFrame([{"dataset": "GSE290035", "comparison": "shZNF469_vs_shControl", "biological_scope": "single CCD-1064Sk cell line; three technical RNA libraries per condition", **score_effect}]).to_csv(OUT / "GSE290035_frozen_response_effect.csv", index=False)
    return {"count_matrix_sha256": sha256(path), "fixed_genes_detected": len(eligible), "score_difference": score_effect["difference"], "score_p_value_technical": score_effect["p_value"]}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    model = pd.read_csv(MODEL_PATH)
    if list(model.gene) != list(UNIPROT):
        raise ValueError("Frozen model genes differ from the predeclared UniProt mapping")
    pd.DataFrame([{"gene": gene, "reviewed_human_canonical_accession": accession, "mapping_source": "UniProt reviewed human canonical entry; frozen 2026-07-14"} for gene, accession in UNIPROT.items()]).to_csv(ROOT / "config" / "fssi_uniprot_mapping.csv", index=False)
    summary = {
        "analysis_version": "final release",
        "preanalysis_lock": "docs/62_raw_proteomics_and_znf469_preanalysis_lock.md",
        "fixed_model_sha256": sha256(MODEL_PATH),
        "PXD066833": analyse_pxd066833(model),
        "PXD039729": analyse_pxd039729(model),
        "GSE290035": analyse_gse290035(model),
    }
    (OUT / "public_expansion_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
