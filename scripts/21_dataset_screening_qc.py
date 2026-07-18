"""Create the accession screening audit and summarise new QC sensitivities."""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "tables"
CONFIG = ROOT / "config"


def exact_permutation(a: np.ndarray, b: np.ndarray) -> float:
    pooled = np.concatenate([a, b]); observed = abs(a.mean() - b.mean())
    values = []
    for idx in itertools.combinations(range(len(pooled)), len(a)):
        mask = np.zeros(len(pooled), bool); mask[list(idx)] = True
        values.append(abs(pooled[mask].mean() - pooled[~mask].mean()))
    return float(np.mean(np.asarray(values) >= observed - 1e-12))


def screening_table() -> pd.DataFrame:
    rows = [
        ("GSE7890", "included", "measured hydrocortisone response", "five paired keloid fibroblast strains"),
        ("GSE113619", "included", "paired wound induction", "participant-level longitudinal design"),
        ("GSE130973", "included_context", "healthy reference transfer", "healthy skin only; no disease contrast"),
        ("GSE151177", "included_boundary", "protocol transfer boundary", "no cells passed frozen fibroblast transfer"),
        ("GSE156326", "included", "hypertrophic-scar frozen projection", "three donors per group"),
        ("GSE158395", "eligible_not_run", "population/reference sensitivity", "held for a future ancestry-focused analysis; not required for core claims"),
        ("GSE163973", "included_development", "state development", "donor mapping sensitivity retained"),
        ("GSE165816", "excluded_adjacent_disease", "diabetic-foot-ulcer context", "impaired healing rather than pathological scar"),
        ("GSE166120", "excluded_adjacent_disease", "diabetic-foot-ulcer spatial context", "GeoMx impaired-healing panel outside core scar question"),
        ("GSE178411", "included", "clinical injury-to-scar trajectory", "patient-condition aggregation and composition sensitivity"),
        ("GSE181297", "included_development", "state and spatial development", "one participant in communication inference"),
        ("GSE181316", "included", "keloid frozen projection", "anatomical samples collapsed to donor"),
        ("GSE181318", "excluded_adjacent_resource", "neural context", "overlapping source study; not an independent disease cohort"),
        ("GSE188952", "included_exploratory", "scar-type specificity", "three to five samples per scar type"),
        ("GSE191067", "included", "four-tissue reference sensitivity", "12 matrix samples; donor pairing unresolved"),
        ("GSE210434", "eligible_not_run", "culture-domain scar-type contrast", "cultured fibroblasts cannot support tissue-level transfer"),
        ("GSE220300", "excluded_adjacent_resource", "Schwann-cell context", "neural component rather than independent fibroblast-state validation"),
        ("GSE227391", "excluded_feature_coverage", "targeted screen", "33-gene panel has zero overlap with 18 FSSI genes"),
        ("GSE237754", "excluded_small_n", "case-level scRNA context", "one keloid and one control library"),
        ("GSE243716", "included_development", "direct keloid-hypertrophic context", "same-participant two-tissue design"),
        ("GSE246562", "included", "mechanical stiffness response", "three paired donors per disease group"),
        ("GSE253664", "excluded_metadata", "culture context", "GEO source metadata internally inconsistent"),
        ("GSE261116", "excluded_small_n", "radiation response", "one library per condition"),
        ("GSE262112", "excluded_cell_line", "technical culture response", "cell lines with technical repeats"),
        ("GSE266334", "excluded_imbalanced", "PIEZO2 scRNA context", "three keloids and one normal; not balanced validation"),
        ("GSE266338", "included_descriptive", "PIEZO2-positive fibroblast context", "one profile per group; descriptive only"),
        ("GSE270438", "excluded_cell_selection", "immune-side context", "CD45-positive enrichment unsuitable for fibroblast projection"),
        ("GSE274709", "excluded_unrecomputable", "PIEZO2 and recurrence context", "public supplement contains chromosome idxstats only; no sample-level gene matrix or recurrence-to-GSM mapping"),
        ("GSE282885", "included_boundary", "post-lock disease transfer", "library-level; donor independence undocumented"),
        ("GSE293677", "included", "measured siRUNX2 response", "three paired profiles"),
        ("GSE293834", "excluded_small_n", "culture shift", "one matched donor"),
        ("GSE303486", "included", "measured drug screen", "25 conditions across three plates"),
        ("GSE303487", "included_descriptive", "focused FR-1 response", "two profiles per condition"),
        ("GSE307210", "included", "measured siPOSTN response", "three profiles per group; pairing undocumented"),
        ("GSE307504", "included_context", "coordinate-aware spatial boundary", "retrospective spatial context"),
        ("GSE320017", "included_context", "longitudinal laser context", "four participants across repeated times"),
        ("GSE335482", "included_boundary", "paired lesion-zone transfer", "four paired participants"),
        ("GSE44270", "included_lower_tier", "legacy bulk/culture context", "experiment-standardised rather than frozen projection"),
        ("GSE92566", "included_context", "whole-tissue deconvolution", "composition sensitivity; not single-cell transfer"),
    ]
    out = pd.DataFrame(rows, columns=["accession", "disposition", "analysis_role", "reason_or_guardrail"])
    out.insert(1, "source", "GEO")
    out["screening_audit_date"] = "2026-07-12"
    out["audit_status"] = "retrospective_complete_candidate_set_not_prospective_PRISMA_registration"
    return out


def qc_count_audit() -> None:
    qc = pd.read_csv(T / "GSE191067_cell_qc_audit.csv.gz")
    steps = []
    mask = np.ones(len(qc), bool)
    for label, criterion in [
        ("source_matrix", np.ones(len(qc), bool)),
        ("UMI_below_8000", qc.nCount_RNA.to_numpy() < 8000),
        ("genes_at_least_500", qc.nFeature_RNA.to_numpy() >= 500),
        ("genes_at_most_4000", qc.nFeature_RNA.to_numpy() <= 4000),
        ("mitochondrial_below_10_percent", qc.percent_mt.to_numpy() < 10),
    ]:
        mask &= criterion
        steps.append({"step": label, "cells_remaining": int(mask.sum()), "cells_removed_at_step": len(qc) - int(mask.sum()) if label == "source_matrix" else np.nan})
    for i in range(1, len(steps)):
        steps[i]["cells_removed_at_step"] = steps[i - 1]["cells_remaining"] - steps[i]["cells_remaining"]
    out = pd.DataFrame(steps)
    out["source_publication_reported_cells"] = 100987
    out["difference_final_vs_reported"] = out.cells_remaining - 100987
    out.to_csv(T / "GSE191067_cell_count_reconciliation.csv", index=False)


def singlet_effects() -> None:
    scores = pd.read_csv(T / "GSE191067_singlet_only_sample_scores.csv")
    rows = []
    for positive, reference in [("keloid", "normal_skin"), ("keloid", "normal_scar"), ("keloid", "perilesional_skin")]:
        a = scores.loc[scores.condition.eq(positive), "fssi_frozen_weighted"].to_numpy(float)
        b = scores.loc[scores.condition.eq(reference), "fssi_frozen_weighted"].to_numpy(float)
        difference = a.mean() - b.mean()
        se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        df = (a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) ** 2 / (
            (a.var(ddof=1) / len(a)) ** 2 / (len(a) - 1) + (b.var(ddof=1) / len(b)) ** 2 / (len(b) - 1))
        crit = stats.t.ppf(.975, df)
        rows.append({"contrast": f"{positive}_minus_{reference}", "n_positive": len(a), "n_reference": len(b),
                     "mean_difference": difference, "ci95_low": difference - crit * se, "ci95_high": difference + crit * se,
                     "exact_permutation_p": exact_permutation(a, b), "scope": "locked_feature_scDblFinder_singlet_sensitivity"})
    pd.DataFrame(rows).to_csv(T / "GSE191067_singlet_only_effects.csv", index=False)


def main() -> None:
    audit = screening_table()
    audit.to_csv(CONFIG / "dataset_screening_audit.csv", index=False)
    flow = audit.assign(included=audit.disposition.str.startswith("included")).groupby("included", as_index=False).size()
    pd.DataFrame([{"candidate_accessions_screened": len(audit),
                   "included_or_contextual": int(audit.disposition.str.startswith("included").sum()),
                   "excluded_or_not_run": int((~audit.disposition.str.startswith("included")).sum()),
                   "development_accessions": int(audit.disposition.eq("included_development").sum()),
                   "measured_response_accessions": int(audit.analysis_role.str.contains("measured|drug|siRUNX|siPOSTN", case=False).sum()),
                   "screening_scope": "GEO candidate accessions accumulated across project search streams",
                   "registration_status": "retrospective audit; not prospective PRISMA registration"}]).to_csv(
                       T / "dataset_screening_flow_summary.csv", index=False)
    pd.DataFrame([{
        "accession": "GSE274709", "raw_archive_bytes": (ROOT / ".stage" / "clinical" / "GSE274709" / "GSE274709_RAW.tar").stat().st_size,
        "raw_archive_content": "30 chromosome-level idxstats files", "gene_level_expression_available": False,
        "sample_level_recurrence_mapping_available": False,
        "decision": "exclude from FSSI calculation; retain as clinical literature context",
    }]).to_csv(T / "GSE274709_reaudit.csv", index=False)
    qc_count_audit(); singlet_effects()
    print(audit.disposition.value_counts().to_string())


if __name__ == "__main__":
    main()
