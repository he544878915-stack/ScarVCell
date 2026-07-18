# ScarVCell

ScarVCell is the final reproducibility release accompanying the manuscript **An auditable fibroblast repair-state coordinate separates established lesion state from inherited keloid susceptibility and measured transcriptional response**.

The release contains the frozen 18-gene Fibroblast Scar-State Index (FSSI), final public-data analysis scripts, accession and unit audits, prespecified analysis locks, machine-readable result tables, the supplementary data workbook and publication figures. It does not contain manuscript drafts, language-editing records, submission utilities or superseded analytical outputs.

## Clinical scope

The analysis tests whether a matrix-rich fibroblast state can be distinguished from comparator-tissue effects, dermal composition, inherited susceptibility and measured molecular response in public pathological-scar resources. FSSI is a research stratification coordinate, not a diagnostic test, prognostic model or treatment recommendation.

## Repository structure

- `analysis_manifest.csv`: ordered index of released analysis modules
- `analysis_locks/`: prespecified rules for added genetic, atlas, transcriptomic and proteomic evidence
- `config/`: frozen FSSI specification, dataset audit, provenance and UniProt mapping
- `scripts/`: final result-generating analysis scripts only
- `results/tables/`: final reported tables and unit-level supporting values
- `results/figures/`: nine main and eleven supplementary figures
- `environment/`: Python, R and external command-line dependencies
- `SHA256SUMS.csv`: integrity record for every released file

## Reproduction boundaries

Original GEO, PRIDE/ProteomeXchange, GWAS and atlas source files are referenced but not redistributed. Raw inputs should be retrieved from the accessions in the dataset audit. Patient, donor, specimen, cell-source, matrix-library and accession units remain distinct. Expression values are not pooled across accessions for inferential testing, and cells are not treated as biological replicates.

## Study team

Junwei He led the study and made the principal contribution, including study conception, methodology, software, public-data curation, primary computational analysis and visualisation. Zezhao Ding contributed independent output checks, data and literature curation, consistency review and interpretation. Wei Xu contributed dataset, unit and clinical-phenotype audit, numerical verification and clinical interpretation. Lei Yi supervised the clinical framing and interpretation. The manuscript contains the complete author-contribution statement.

## Citation

Use `CITATION.cff` and the version-specific Zenodo record supplied with this release.

## Funding

This work was supported by the National Natural Science Foundation of China (Grant No. 82402902). The funder had no role in study design, data analysis or interpretation, manuscript preparation or the decision to submit.
