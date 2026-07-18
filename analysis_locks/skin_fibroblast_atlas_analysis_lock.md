# Skin-fibroblast atlas analysis lock

Date locked: 2026-07-13 (Asia/Shanghai)

## Medical question

Does the frozen Fibroblast Scar-State Index (FSSI) localise to a published consensus skin-fibroblast state, and is the same programme restricted to pathological scars or shared by other skin diseases?

## Public source

- Steele et al., *Nature Immunology* (2025), DOI: 10.1038/s41590-025-02267-8.
- Processed data: https://cellatlas.io/studies/skin-fibroblast
- Analysis code and published marker tables: https://github.com/haniffalab/skin_fibroblast_atlas

The downloaded snapshot, byte counts and cryptographic hashes will be recorded before analysis.

## Independence boundary

The published atlas contains GSE181297, GSE181316 and GSE163973, which already contribute to ScarVCell development or contextual analyses. These accessions, and every other accession already used by ScarVCell, will be flagged before scoring.

- Analyses containing overlapping accessions are ontology crosswalks only.
- Primary cross-disease transport analyses exclude every overlapping accession.
- No result from this atlas will be called independent keloid replication unless a non-overlapping keloid cohort is present and donor identity is auditable.

## Frozen quantities

- The 18 FSSI genes, directions and weights remain unchanged.
- No atlas outcome may alter the FSSI model.
- The pathological and repair components remain the pre-existing fixed gene partitions.
- F1-F8 or other consensus subtype labels remain exactly as published.

## Analysis units

The preferred biological unit is donor within accession and disease/site context. If only sample identifiers are released, the analysis unit will be called a sample rather than a donor. Cells are never treated as independent replicates.

## Primary analyses

1. **Consensus-state localisation.** Calculate donor/sample pseudobulk FSSI, pathological-component and repair-component scores within each published fibroblast subtype. Estimate within-donor subtype contrasts when at least two subtypes are present.
2. **Non-overlapping cross-disease transport.** Remove all ScarVCell-overlapping accessions before ranking diseases by donor/sample-level FSSI. Report distributions and uncertainty, not cell-level P values.
3. **Scar specificity.** Quantify whether high FSSI is confined to scar/fibrotic contexts or also occurs in inflammatory, malignant and other skin diseases. A shared signal will be retained as evidence against scar specificity.
4. **Published-marker crosswalk.** Test fixed pathological and repair genes against the authors' published subtype marker tables using the table's tested-gene universe. This is an ontology comparison, not expression replication.

## Scoring hierarchy

- If released expression is compatible with the frozen training representation, the frozen transform will be applied without refitting.
- If representation compatibility cannot be established, the primary atlas result will use within-dataset rank-standardised fixed-weight balance and will be labelled a transport sensitivity, not a frozen numerical projection.
- Missing genes are reported. Weights may be renormalised only over observed genes, with full-gene and coverage sensitivities shown.

## Statistical rules

- Effects are calculated within accession before any summary.
- Donor/sample bootstrap intervals preserve accession structure.
- Dataset-level effects are displayed separately; no cross-dataset cell-level test is permitted.
- Multiple subtype and disease comparisons use Benjamini-Hochberg correction.
- Leave-one-accession-out and overlap-exclusion analyses are mandatory.

## Claim ceiling

The atlas can establish consensus-state localisation and cross-disease specificity or non-specificity. It cannot establish prognosis, treatment efficacy, causal regulation or an independent keloid replication when the keloid accessions overlap with ScarVCell development data.

