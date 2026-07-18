# Multi-repository expansion analysis lock

**Lock date:** 2026-07-13  
**Purpose:** extend the frozen ScarVCell analysis with patient-paired expression, measured pharmacological response and orthogonal protein evidence without changing the 18-gene model.

## Non-negotiable rules

1. `results/tables/fssi_frozen_model.csv` defines the genes, signs, weights, training means and training standard deviations. No feature, sign, weight or scaling parameter will be modified after inspecting the datasets below.
2. The biological unit is the patient or donor. Technical replicates and multiple regions from one person are collapsed or modelled within person before inference.
3. Statistical contrasts are calculated within each dataset. Raw samples from different platforms are never pooled into one test.
4. All eligible datasets are retained regardless of direction. Missing genes, ambiguous pairing and failed processing are reported as auditable outcomes rather than silently excluded.
5. Whole-tissue expression and protein analyses test transport of a fixed programme, not fibroblast-specific activity, diagnosis, prognosis, treatment efficacy or causality.

## Locked datasets and endpoints

| Dataset | Locked role | Primary unit | Prespecified endpoint |
|---|---|---|---|
| GSE90051 | Primary paired whole-tissue transport | Seven patients | Fixed-weight keloid-minus-adjacent-skin response from the within-array log ratio; exact sign test and t-based CI of patient effects |
| GSE212954 | Primary paired lesion-zone geometry | Patient identifier encoded in sample title, subject to metadata confirmation | Centre-minus-normal and margin-minus-normal fixed-weight responses; centre-minus-margin is secondary |
| GSE151464 | Small paired directional evidence | Three reported pairs | Fixed-weight keloid-minus-corresponding-skin response; descriptive CI and exact sign-flip test |
| GSE83286 | Small paired directional evidence | Three reported pairs | Fixed-weight earlobe-keloid-minus-normal response; descriptive CI and exact sign-flip test |
| GSE282479 | Measured response calibration | Three keloid and four normal fibroblast donors | Paired paricalcitol-minus-untreated fixed-weight response within disease context; interaction described by donor-level response difference, not efficacy |
| GSE145725 | Secondary culture-domain transport | Five keloid and five normal fibroblast lines, with technical replicates collapsed | Fixed-weight keloid-minus-normal response at cell-line level; no tissue-level interpretation |
| PXD015057 | Orthogonal ECM protein crosswalk | Seven keloids and five paired normal-skin/normal-scar patients, as released | Coverage and direction of detected frozen-programme proteins; normal-scar-minus-normal-skin paired contrast and keloid contrasts where sample mapping permits |
| PXD029631 | Orthogonal whole-tissue protein crosswalk | Six reported keloid-adjacent-skin pairs | Detected fixed-programme protein direction and paired effect; no imputation of undetected proteins |
| GSE274709 | Conditional clinical-context analysis | Ten keloid tissues | Raw RNA-seq reprocessing only if GSM-to-patient recurrence mapping is public and unambiguous; otherwise exclusion remains final |

## Scoring and inference

- For separate expression matrices with non-negative expression, genes are transformed according to the source-normalised scale and standardised within experiment before fixed weights are applied. This is a fixed-weight response score, not a calibrated single-cell FSSI projection.
- For GSE90051, each array already represents `log10(keloid/adjacent normal)`. Probe-level ratios are collapsed to one value per gene by the median, and the patient response is the weighted sum of available gene ratios. A coverage sensitivity score rescales absolute weights over observed genes.
- For protein datasets, gene symbols are mapped from released protein identifiers. Only observed proteins contribute. Programme-specific sign concordance and abundance-matched random-set enrichment are reported when the released table supports them.
- Confidence intervals describe variation among independent patients or donors. Exact sign or sign-flip tests are used for small paired samples. With fewer than four independent units, P values are descriptive and no multiplicity-adjusted discovery claim is made.
- A random-effects summary is permitted only across at least three independent, biologically comparable whole-tissue contrasts. Dataset-specific effects remain visible and are the primary evidence.
- GSE173900, GSE218007, GSE232079 and other newly identified accessions enter the expanded screening table. They are analysed only if their public metadata and processed matrices support an unambiguous prespecified patient- or donor-level contrast; otherwise the reason for non-analysis is retained.

## Display and claim guardrails

- - The main display will distinguish `supports direction`, `discordant`, `not tested`, and `not evaluable`.
- Measured paricalcitol response will be labelled `measured response calibration`; it will not be described as treatment validation.
- Protein findings will be labelled `protein-level crosswalk`; incomplete coverage cannot be called FSSI replication.
- GSE274709 cannot be linked to recurrence by inferred sample order, PIEZO2 level or visual matching.
