# Genetics and ancestry analysis lock

Lock date: 12 July 2026. This specification was written before downloading or inspecting
candidate-level outcomes from the keloid GWAS summary statistics or calculating GSE158395 FSSI.

## Medical question

Does germline keloid susceptibility preferentially anchor the frozen matrix-rich fibroblast state,
and does that state retain direction in an African-American paired lesional/non-lesional cohort?

## Immutable transcriptomic objects

- Frozen FSSI genes, signs and weights remain exactly those in `results/tables/fssi_frozen_model.csv`.
- No FSSI gene or weight will be added, removed or refitted after genetic outcome inspection.
- The extended candidate panel is fixed as RUNX2, POSTN, TGFB1, CSF1, CXCL12, MIF, APOD, PI16,
  PTGDS and PDGFB.
- The pathological atlas anchors are clusters 1, 5 and 11; the repair anchor is cluster 2.

## Public genetic inputs

- GCST90652487: 2025 multi-ancestry keloid GWAS meta-analysis.
- GCST90652488: European-ancestry meta-analysis.
- GCST90652489: African-ancestry meta-analysis.
- Published supplementary gene-tissue and colocalisation results are an audit source, not a
  substitute for the locked candidate-level calculations.

## Primary genetic analyses

1. Audit genome build, alleles, effect scale, sample size, variant count and genomic inflation fields.
2. Keep ancestry-specific analyses separate. Multi-ancestry results will not be assigned a single
   ancestry LD reference.
3. Calculate ancestry-appropriate gene-level association for European and African summary
   statistics with MAGMA v1.10, NCBI build-37 protein-coding gene locations, a +/-10 kb gene window
   and the matching 1000 Genomes Phase 3 European or African LD panel. The extended MHC region
   (chromosome 6: 25-34 Mb) is excluded. Sample sizes are fixed from the GWAS metadata.
4. Test the fixed 18-gene FSSI, its pathological and repair components, the eight-gene bootstrap-stable
   core and the fixed extended candidate panel using MAGMA competitive gene-set analysis and 10,000
   size/expression-bin-matched random gene sets. Bonferroni-adjusted gene-level and fixed-set results
   are reported alongside unadjusted values.
5. Apply genetic gene scores to the development single-cell atlas only as a cell-state enrichment
   analysis. Inference is aggregated to donor by cluster; cells are not biological replicates.
6. Report European, African and cross-ancestry direction separately. Discordance is retained.

## Genetic evidence classes

- `colocalised`: candidate has a prespecified relevant skin/fibroblast cis-eQTL signal sharing a
  causal variant with keloid susceptibility under the source study's stated posterior threshold.
- `gene_level`: candidate passes the locked ancestry-specific gene-level threshold but lacks
  relevant-tissue colocalisation.
- `set_level_only`: the locked programme is enriched but the individual candidate is not supported.
- `not_supported`: tested without support.
- `not_evaluable`: absent or technically incompatible.

Variant proximity, pathway membership and literature mention alone cannot be called genetic support.
No broad two-sample Mendelian-randomisation screen will be used.

## GSE158395 analysis

- Processed RNA-seq expression matrix and GEO sample metadata are the public inputs.
- JKR782-001, JKR782-003 and JKR782-004 are the biological keloid participants.
- Two lesional profiles from JKR782-003 are averaged before inference so that the participant,
  rather than the lesion profile, remains the biological unit.
- Primary contrast: paired lesional minus non-lesional skin in three participants.
- Secondary contrast: three participant-level lesional values versus six healthy skin samples.
- Fixed FSSI genes and weights are retained; genes are standardised within experiment because the
  bulk expression scale is not calibrated to the single-cell training representation.
- Composition sensitivity uses a broad-cell reference excluding every FSSI gene. It is reported
  separately and cannot create a cell-intrinsic validation claim.
- Exact sign-flip or label-permutation tests and t-based confidence intervals are both reported.

Metadata amendment, 13 July 2026: GEO releases six normal-skin matrices, whereas the source paper
describes five healthy controls. The six-matrix GEO comparison remains the locked secondary analysis;
a publication-aligned five-control sensitivity excludes GSM4798880, the additional normal matrix not
counted in the paper. Neither comparison affects the primary three-pair lesion/non-lesion analysis.

## Claim ceiling

This layer can provide independent germline anchoring and ancestry-aware transcriptomic consistency.
It cannot establish that a target is causal without relevant-tissue colocalisation, and it cannot
establish recurrence prediction, therapeutic response or ancestry-specific clinical efficacy.

## Locked amendment: externally defined susceptibility-state convergence

Amendment date: 13 July 2026. This amendment was written after auditing the published Greene et al.
supplementary-data schema and before calculating expression of the resulting genes in the ScarVCell
atlas or GSE158395. It does not alter any FSSI object or candidate-level MAGMA test.

- The susceptibility-anchor set is derived programmatically from the source study's multi-ancestry
  S-PrediXcan/colocalisation table. A gene is included only when its tissue is sun-exposed skin,
  non-sun-exposed skin or cultured fibroblasts, its GPGE P value is below the source study's
  relevant-tissue threshold of 1.5e-6, and PP.H4 exceeds 0.90.
- If a gene has multiple qualifying tissues, the result with the largest absolute GPGE Z score
  defines its fixed direction and weight. Genes with conflicting qualifying directions are reported
  and excluded from the signed score. Absolute weights are normalised to sum to one.
- This externally defined signed susceptibility score is evaluated in donor-by-fibroblast-cluster
  pseudobulks. The prespecified contrast is the mean of pathological anchors 1, 5 and 11 minus the
  repair anchor 2 within each donor; cells are never inferential units. Exact donor sign-flip and
  t-based intervals are reported when the required clusters are present.
- A 10,000-iteration expression-matched random-gene null evaluates whether the observed anchor
  contrast exceeds that expected from genes with similar atlas abundance. This is a convergence
  analysis, not genetic validation of FSSI.
- The same frozen susceptibility score is projected to the three GSE158395 lesional/non-lesional
  participant pairs after duplicate-library collapse and experiment-wide gene standardisation.
  Exact paired inference is retained. It cannot be interpreted as ancestry-specific prediction.
- FSSI-susceptibility score correlation is descriptive and is reported at donor or participant level,
  never at cell level.
