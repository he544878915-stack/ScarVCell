# Proteomics and ZNF469 analysis lock

Lock date: 2026-07-14 (Asia/Shanghai)

This document was written after repository-level feasibility inspection but before any ScarVCell fixed-gene protein effect or GSE290035 frozen-score outcome was calculated. Published source-study conclusions are known; no FSSI-specific result from the added records has been inspected.

## Scientific purpose

The primary aim is to determine whether public raw or search-output proteomics can upgrade the existing source-level protein crosswalk to sample-level, within-study evidence. The secondary aim is to test whether an independently published ZNF469 knockdown produces a broad FSSI shift or a narrower extracellular-matrix response. Neither analysis may be used to refit the 18-gene FSSI, alter its weights, or discover a favourable endpoint after outcome inspection.

## Locked resources and roles

| Resource | Locked role | Primary analysis unit | Claim ceiling |
|---|---|---|---|
| PXD066833 | Sample-level clinical-scar proteomics, conditional on auditable sample mapping | One released biological tissue specimen | Within-study protein abundance direction; no clinical prediction |
| PXD039729 / MSV000091168 | Independent keloid-versus-normal proteomic context, conditional on auditable tissue and fraction mapping | One biological specimen within a comparable fraction | Within-study protein abundance direction; no pooling across fractions or diagnoses |
| GSE290035 | Technical perturbation boundary after shZNF469 | One cell-line experiment; three released profiles per arm are technical replicates | Technical expression response only; not donor-level validation |

## Fixed features and directions

The frozen FSSI genes, pathological/repair assignment, signs and weights are read from the existing locked model files. No gene may be added, removed or reweighted based on these resources. Protein aliases may be mapped only through auditable gene/protein identifiers; ambiguous mappings are reported as not evaluable.

Expected favourable protein direction is increased abundance for pathological genes and decreased abundance for repair genes in scar relative to the study-specific reference. This expectation is used only to orient effects. Both concordant and discordant proteins are retained.

For GSE290035, positive control-minus-shZNF469 FSSI denotes movement away from the matrix-rich state. Pathological and repair components, individual fixed genes and prespecified ECM/collagen pathways are reported together. A local ECM response without a broad FSSI response is retained as a valid boundary result.

## Inclusion and failure rules

### PXD066833

1. A public search-output or quantitative table must map protein identifiers to individual raw acquisitions.
2. Scar type and reference type must be recoverable from repository metadata or the source-study sample table without guessing from file order.
3. Technical injections are collapsed before inference.
4. Keloid, hypertrophic scar and normal skin are not silently pooled. The primary contrast is keloid versus normal skin if both are identifiable; hypertrophic scar is separate.
5. If only group-level differential results are public, the resource remains a source-level crosswalk and receives no sample-level confidence interval or P value.

### PXD039729

1. Keloid and normal-skin biological specimens must be identifiable.
2. Cytosolic, membrane, nuclear and cytoskeletal fractions are analysed separately or combined only after within-sample normalisation with fraction represented in the model.
3. Folliculitis keloidalis nuchae, dermatofibrosarcoma protuberans and fibrosarcoma are excluded from the primary keloid contrast and may appear only as labelled specificity contexts.
4. If biological specimen identifiers cannot be distinguished from injections, no sample-level test is performed.

### GSE290035

1. All six released profiles are retained.
2. The three profiles per arm are treated as technical replication of one experiment, exactly as described by GEO.
3. T-based intervals describe technical variability only and are never labelled biological confidence intervals.
4. The experiment cannot promote ZNF469 into an independently validated therapeutic tier.

## Statistical reporting

Protein abundance is log2 transformed when required and normalised within study using the source-compatible workflow. Fixed-gene effects are calculated separately within each dataset. For auditable independent biological specimens, Welch mean differences and 95% confidence intervals are reported; permutation P values are used when the design permits. Multiple fixed proteins are adjusted within study by Benjamini-Hochberg, with raw effect direction retained regardless of significance.

A protein-level programme summary is secondary and is calculated only when at least four pathological and three repair proteins are quantified in the same sample matrix. Available fixed weights are renormalised over measured proteins. It is labelled a partial protein projection, not the RNA-trained FSSI itself. No raw abundance is pooled across proteomic studies.

Cross-study synthesis, if feasible, uses per-protein standardised effects and displays each study separately. A random-effects estimate is secondary and is omitted when fewer than two comparable study effects exist.

