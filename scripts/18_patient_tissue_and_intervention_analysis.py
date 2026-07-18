"""Analyse the final release public clinical, mechanical and pharmacological expansions.

The FSSI genes, signs and weights are read from the frozen FSSI model. Assay
units are standardised only within each external experiment; outputs are
therefore fixed-weight response coordinates, not refitted single-cell scores.
"""

from __future__ import annotations

import gzip
import itertools
import json
import math
import os
import re
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls


ROOT = Path(__file__).resolve().parents[1]
STAGE = Path(os.environ.get("SCARVCELL_STAGE", str(ROOT / ".stage"))) / "public_expansion"
TABLES = ROOT / "results" / "tables"
SEED = 20260712
RNG = np.random.default_rng(SEED)

URLS = {
    "gene_info": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    "GSE178411": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178411/suppl/GSE178411_counts.txt.gz",
    "GSE113619": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113619/suppl/GSE113619_RNA-seq_keloids_raw.csv.gz",
    "GSE303486_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303486/suppl/GSE303486_Drug-seq2_count.csv.gz",
    "GSE303486_map": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303486/suppl/GSE303486_SampleName_barcode.csv.gz",
    "GSE303487": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303487/suppl/GSE303487_RAW.tar",
    "GSE246562": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246562/suppl/GSE246562_All_Probes_Raw_counts.csv.gz",
    "GSE188952": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188952/suppl/GSE188952_Processed_FPKM.tsv.gz",
}


def download(name: str, suffix: str) -> Path:
    STAGE.mkdir(parents=True, exist_ok=True)
    path = STAGE / f"{name}{suffix}"
    if not path.exists():
        print(f"Downloading {name}", flush=True)
        urllib.request.urlretrieve(URLS[name], path)
    return path


def gene_id_map() -> dict[str, str]:
    path = download("gene_info", ".gene_info.gz")
    info = pd.read_csv(path, sep="\t", compression="gzip", dtype=str, low_memory=False)
    return dict(zip(info["GeneID"], info["Symbol"]))


def collapse_symbol_rows(frame: pd.DataFrame, mapping: dict[str, str] | None = None) -> pd.DataFrame:
    symbols = frame.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if mapping is not None:
        symbols = symbols.map(mapping)
    frame = frame.copy()
    frame.index = symbols
    frame = frame[frame.index.notna() & (frame.index != "-")]
    return frame.groupby(level=0).sum()


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).clip(lower=0)
    library = counts.sum(axis=0).replace(0, np.nan)
    return np.log2(counts.div(library, axis=1) * 1e6 + 0.5)


def score_matrix(expr: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    genes = [g for g in model.gene if g in expr.index]
    x = expr.loc[genes].T
    z = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1).replace(0, 1).fillna(1)
    m = model.set_index("gene").loc[genes]
    weight = m.weight.to_numpy(float).copy(); weight /= np.abs(weight).sum()
    equal = m.equal_weight.to_numpy(float).copy(); equal /= np.abs(equal).sum()
    pathological = [g for g in genes if m.loc[g, "program"] == "pathological"]
    repair = [g for g in genes if m.loc[g, "program"] == "repair"]
    out = pd.DataFrame(index=x.index)
    out["fixed_weight_response"] = z.to_numpy() @ weight
    out["equal_weight_response"] = z.to_numpy() @ equal
    out["pathological_component"] = z[pathological].mean(axis=1)
    out["repair_component_reversed"] = -z[repair].mean(axis=1)
    fib5 = [g for g in ["ADAM12", "COMP", "POSTN"] if g in z]
    out["Fib_5_module"] = z[fib5].mean(axis=1) if fib5 else np.nan
    out["genes_present"] = len(genes)
    out["genes_total"] = len(model)
    out.index.name = "sample_id"
    return out.reset_index()


def welch_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    effect = float(a.mean() - b.mean())
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    se2 = va / len(a) + vb / len(b)
    if se2 <= 0:
        return effect, effect, effect, 1.0
    df = se2**2 / ((va / len(a))**2 / (len(a) - 1) + (vb / len(b))**2 / (len(b) - 1))
    q = stats.t.ppf(.975, df)
    p = 2 * stats.t.sf(abs(effect / math.sqrt(se2)), df)
    return effect, effect - q * math.sqrt(se2), effect + q * math.sqrt(se2), float(p)


def paired_ci(d: np.ndarray) -> tuple[float, float, float, float]:
    d = np.asarray(d, float)
    effect = float(d.mean())
    if len(d) < 2 or np.std(d, ddof=1) == 0:
        return effect, effect, effect, 1.0
    se = np.std(d, ddof=1) / math.sqrt(len(d))
    q = stats.t.ppf(.975, len(d) - 1)
    p = 2 * stats.t.sf(abs(effect / se), len(d) - 1)
    return effect, effect - q * se, effect + q * se, float(p)


def independent_permutation(a: np.ndarray, b: np.ndarray, iterations: int = 20000) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    pooled = np.concatenate([a, b]); n = len(a)
    observed = abs(a.mean() - b.mean())
    combinations = math.comb(len(pooled), n)
    if combinations <= 100000:
        values = []
        for idx in itertools.combinations(range(len(pooled)), n):
            mask = np.zeros(len(pooled), bool); mask[list(idx)] = True
            values.append(abs(pooled[mask].mean() - pooled[~mask].mean()))
        return float(np.mean(np.asarray(values) >= observed - 1e-12))
    exceed = 0
    for _ in range(iterations):
        perm = RNG.permutation(pooled)
        exceed += abs(perm[:n].mean() - perm[n:].mean()) >= observed - 1e-12
    return float((exceed + 1) / (iterations + 1))


def paired_signflip(d: np.ndarray) -> float:
    d = np.asarray(d, float); observed = abs(d.mean())
    if len(d) <= 18:
        null = [abs(np.mean(d * np.asarray(s))) for s in itertools.product([-1, 1], repeat=len(d))]
        return float(np.mean(np.asarray(null) >= observed - 1e-12))
    exceed = 0
    for _ in range(20000):
        exceed += abs(np.mean(d * RNG.choice([-1, 1], len(d)))) >= observed - 1e-12
    return float((exceed + 1) / 20001)


def gds_samples(accession: str) -> pd.DataFrame:
    query = urllib.parse.urlencode({"db": "gds", "retmode": "json", "term": f"{accession}[ACCN]"})
    with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{query}") as response:
        uid = json.load(response)["esearchresult"]["idlist"][0]
    with urllib.request.urlopen(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&id={uid}") as response:
        samples = json.load(response)["result"][uid]["samples"]
    return pd.DataFrame(samples).rename(columns={"accession": "gsm", "title": "sample_title"})


def deconvolve(expr: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    donor = pd.read_csv(TABLES / "scrna_broad_celltype_donor_logcpm_reference.csv", index_col=0)
    meta = pd.read_csv(TABLES / "scrna_broad_celltype_donor_reference_metadata.csv")
    markers = pd.read_csv(TABLES / "scrna_broad_celltype_specificity_markers.csv")
    selected = markers[(markers.rank_within_type <= 50) & ~markers.gene.isin(set(model.gene))]
    genes = sorted(set(selected.gene) & set(donor.index) & set(expr.index))
    types = sorted(meta.broad_cell_type.unique())
    signature = np.column_stack([
        donor.loc[genes, meta.loc[meta.broad_cell_type.eq(ct), "reference_unit"]].mean(axis=1).to_numpy()
        for ct in types
    ])
    signature = np.column_stack([stats.rankdata(signature[:, j]) / len(genes) for j in range(len(types))])
    rows = []
    for sample in expr.columns:
        vector = stats.rankdata(expr.loc[genes, sample].to_numpy()) / len(genes)
        coef, residual = nnls(signature, vector)
        coef = coef / coef.sum() if coef.sum() else coef
        for ct, value in zip(types, coef):
            rows.append({"sample_id": sample, "broad_cell_type": ct, "proportion": value,
                         "fit_rmse": residual / math.sqrt(len(genes)), "marker_genes": len(genes)})
    return pd.DataFrame(rows)


def analyse_gse178411(model: pd.DataFrame, entrez: dict[str, str]) -> None:
    path = download("GSE178411", "_counts.txt.gz")
    with gzip.open(path, "rt") as handle:
        samples = handle.readline().rstrip("\n\r").split("\t")
    raw = pd.read_csv(path, sep="\t", compression="gzip", header=None, skiprows=1)
    counts = raw.copy()
    counts.index = raw.iloc[:, 0].astype(str)
    counts = counts.iloc[:, 1:]
    counts.columns = samples
    counts = collapse_symbol_rows(counts, entrez)
    expr = log_cpm(counts)
    score = score_matrix(expr, model)

    metadata = gds_samples("GSE178411")
    metadata["sample_id"] = metadata.sample_title.str.extract(r"^([^:]+)")[0]
    metadata["patient_id"] = "P" + metadata.sample_title.str.extract(r"Patient\s+(\d+)")[0]
    labels = {
        "Normal skin": "normal_skin", "Normal scar": "normal_scar", "HTS": "hypertrophic_scar",
        "Early Wound": "early_wound", "Late wound": "late_wound", "Chronic wound": "chronic_wound",
    }
    metadata["condition"] = metadata.sample_title.map(
        lambda x: next((value for key, value in labels.items() if key in x), "unresolved"))
    metadata.to_csv(TABLES / "GSE178411_patient_metadata.csv", index=False)
    score = score.merge(metadata[["sample_id", "patient_id", "condition"]], on="sample_id", how="left")
    score.to_csv(TABLES / "GSE178411_sample_scores.csv", index=False)

    patient = score.groupby(["patient_id", "condition"], as_index=False).agg(
        fixed_weight_response=("fixed_weight_response", "mean"),
        equal_weight_response=("equal_weight_response", "mean"),
        pathological_component=("pathological_component", "mean"),
        repair_component_reversed=("repair_component_reversed", "mean"),
        tissue_samples=("sample_id", "size"), genes_present=("genes_present", "min"))
    patient.to_csv(TABLES / "GSE178411_patient_condition_scores.csv", index=False)

    effects = []
    for endpoint, positive, reference in [
        ("hypertrophic_scar_vs_normal_skin", "hypertrophic_scar", "normal_skin"),
        ("hypertrophic_scar_vs_normal_scar", "hypertrophic_scar", "normal_scar"),
        ("late_or_chronic_wound_vs_normal_skin", "late_or_chronic_wound", "normal_skin"),
    ]:
        use = patient.copy()
        if positive == "late_or_chronic_wound":
            use["condition"] = use.condition.replace({"late_wound": positive, "chronic_wound": positive})
            use = use.groupby(["patient_id", "condition"], as_index=False).fixed_weight_response.mean()
        a = use.loc[use.condition.eq(positive), "fixed_weight_response"].to_numpy(float)
        b = use.loc[use.condition.eq(reference), "fixed_weight_response"].to_numpy(float)
        effect, lo, hi, p_t = welch_ci(a, b)
        effects.append({"dataset_id": "GSE178411", "endpoint": endpoint, "unit": "patient_condition",
                        "n_positive": len(a), "n_reference": len(b), "mean_difference": effect,
                        "ci95_low": lo, "ci95_high": hi, "welch_p": p_t,
                        "permutation_p_two_sided": independent_permutation(a, b)})
    pd.DataFrame(effects).to_csv(TABLES / "GSE178411_patient_effects.csv", index=False)

    proportions = deconvolve(expr, model)
    proportions.to_csv(TABLES / "GSE178411_reference_deconvolution.csv", index=False)
    fibro = proportions.query("broad_cell_type == 'fibroblast'")[["sample_id", "proportion"]]
    adjusted = score.merge(fibro, on="sample_id")
    design = np.column_stack([np.ones(len(adjusted)), adjusted.proportion])
    adjusted["composition_residual"] = adjusted.fixed_weight_response - design @ np.linalg.lstsq(
        design, adjusted.fixed_weight_response, rcond=None)[0]
    adjusted_patient = adjusted.groupby(["patient_id", "condition"], as_index=False).agg(
        composition_residual=("composition_residual", "mean"), fibroblast_proportion=("proportion", "mean"),
        tissue_samples=("sample_id", "size"))
    adjusted_patient.to_csv(TABLES / "GSE178411_composition_adjusted_patient_scores.csv", index=False)
    a = adjusted_patient.loc[adjusted_patient.condition.eq("hypertrophic_scar"), "composition_residual"].to_numpy(float)
    b = adjusted_patient.loc[adjusted_patient.condition.eq("normal_skin"), "composition_residual"].to_numpy(float)
    effect, lo, hi, p_t = welch_ci(a, b)
    rho, rho_p = stats.spearmanr(adjusted.fixed_weight_response, adjusted.proportion)
    pd.DataFrame([{"dataset_id": "GSE178411", "endpoint": "hypertrophic_scar_vs_normal_skin",
                   "samples": len(adjusted), "patient_units_positive": len(a), "patient_units_reference": len(b),
                   "spearman_fssi_fibroblast_fraction": rho, "spearman_p": rho_p,
                   "composition_residual_mean_difference": effect, "ci95_low": lo, "ci95_high": hi,
                   "welch_p": p_t, "permutation_p_two_sided": independent_permutation(a, b)}]).to_csv(
                       TABLES / "GSE178411_composition_sensitivity_summary.csv", index=False)


def analyse_gse113619(model: pd.DataFrame, entrez: dict[str, str]) -> None:
    path = download("GSE113619", "_raw.csv.gz")
    raw = pd.read_csv(path, compression="gzip", index_col=0)
    raw.index = raw.index.astype(str)
    counts = collapse_symbol_rows(raw, entrez)
    collapsed = {}
    metadata = []
    for participant in [f"K{i}" for i in range(1, 9)] + [f"N{i}" for i in range(1, 7)]:
        for time_label, token in [("baseline", "1st"), ("six_weeks", "2nd")]:
            columns = [c for c in counts.columns if c.startswith(f"{participant}-{token}_")]
            if columns:
                unit = f"{participant}_{time_label}"
                collapsed[unit] = counts[columns].sum(axis=1)
                metadata.append({"sample_id": unit, "participant_id": participant,
                                 "susceptibility": "keloid_prone" if participant.startswith("K") else "healthy",
                                 "time": time_label, "sequencing_runs": len(columns)})
    expr = log_cpm(pd.DataFrame(collapsed))
    score = score_matrix(expr, model).merge(pd.DataFrame(metadata), on="sample_id")
    score.to_csv(TABLES / "GSE113619_participant_scores.csv", index=False)
    wide = score.pivot(index="participant_id", columns="time", values="fixed_weight_response").dropna()
    delta = (wide.six_weeks - wide.baseline).rename("paired_change").reset_index()
    delta["susceptibility"] = np.where(delta.participant_id.str.startswith("K"), "keloid_prone", "healthy")
    delta.to_csv(TABLES / "GSE113619_paired_changes.csv", index=False)
    a = delta.loc[delta.susceptibility.eq("keloid_prone"), "paired_change"].to_numpy(float)
    b = delta.loc[delta.susceptibility.eq("healthy"), "paired_change"].to_numpy(float)
    effect, lo, hi, p_t = welch_ci(a, b)
    pd.DataFrame([{"dataset_id": "GSE113619", "endpoint": "susceptibility_by_time_interaction",
                   "keloid_prone_n": len(a), "healthy_n": len(b), "difference_in_paired_change": effect,
                   "ci95_low": lo, "ci95_high": hi, "welch_p": p_t,
                   "exact_label_permutation_p": independent_permutation(a, b)}]).to_csv(
                       TABLES / "GSE113619_interaction_summary.csv", index=False)


def parse_drug_profile(name: str) -> tuple[str, str]:
    match = re.match(r"(.+)_([A-H]\d+)_rep(\d+)$", name)
    if not match:
        raise ValueError(f"Unrecognised DRUG-seq profile: {name}")
    return match.group(1), f"rep{match.group(3)}"


def exact_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    observed = float(stats.spearmanr(x, y).statistic)
    null = [abs(stats.spearmanr(x, perm).statistic) for perm in itertools.permutations(y)]
    return observed, float(np.mean(np.asarray(null) >= abs(observed) - 1e-12))


def analyse_gse303486(model: pd.DataFrame) -> None:
    path = download("GSE303486_counts", "_counts.csv.gz")
    counts = pd.read_csv(path, compression="gzip", index_col=0, encoding="latin1")
    counts.index = counts.index.astype(str).str.upper()
    counts = counts.groupby(level=0).sum()
    expr = log_cpm(counts)
    score = score_matrix(expr, model)
    parsed = score.sample_id.map(parse_drug_profile)
    score["compound_raw"] = parsed.map(lambda x: x[0])
    score["replicate"] = parsed.map(lambda x: x[1])
    aliases = {"ICG_001": "ICG-001", "5_FU": "5-FU", "L_Ascorbic acid": "L-Ascorbic acid",
               "SAR_100842": "SAR-100842", "TGFÎ²1": "TGFbeta1", "TGF�1": "TGFbeta1"}
    score["compound"] = score.compound_raw.replace(aliases)
    score.to_csv(TABLES / "GSE303486_profile_scores.csv", index=False)

    replicate = score.groupby(["replicate", "compound"], as_index=False).agg(
        fixed_weight_response=("fixed_weight_response", "mean"),
        pathological_component=("pathological_component", "mean"),
        repair_component_reversed=("repair_component_reversed", "mean"),
        profiles=("sample_id", "size"))
    dmso = replicate[replicate.compound.eq("DMSO")].set_index("replicate")
    rows = []
    for compound, data in replicate[~replicate.compound.eq("DMSO")].groupby("compound"):
        d = data.set_index("replicate").join(dmso[["fixed_weight_response"]], rsuffix="_dmso")
        reduction = d.fixed_weight_response_dmso - d.fixed_weight_response
        effect, lo, hi, p_t = paired_ci(reduction.to_numpy())
        rows.append({"compound": compound, "replicate_n": len(reduction),
                     "mean_measured_fssi_reduction": effect, "ci95_low": lo, "ci95_high": hi,
                     "paired_t_p": p_t, "exact_signflip_p_two_sided": paired_signflip(reduction.to_numpy()),
                     "direction_support_fraction": float(np.mean(reduction > 0))})
    effects = pd.DataFrame(rows).sort_values("mean_measured_fssi_reduction", ascending=False)
    effects.to_csv(TABLES / "GSE303486_compound_effects.csv", index=False)

    lincs = pd.read_csv(TABLES / "drug_reversal_enrichr_lincs_combined_ranked_curated.csv")
    lincs["normalised_name"] = lincs.perturbagen.str.lower().str.replace(r"[^a-z0-9]", "", regex=True)
    effects["normalised_name"] = effects.compound.str.lower().str.replace(r"[^a-z0-9]", "", regex=True)
    matched = effects.merge(lincs.drop_duplicates("normalised_name"), on="normalised_name", how="inner",
                            suffixes=("_measured", "_lincs"))
    matched["lincs_rank"] = matched.normalised_name.map(
        {name: i + 1 for i, name in enumerate(lincs.normalised_name)})
    rho, p = exact_spearman(matched.drug_reversal_score.to_numpy(float),
                            matched.mean_measured_fssi_reduction.to_numpy(float))
    matched.to_csv(TABLES / "GSE303486_LINCS_measured_overlap.csv", index=False)
    pd.DataFrame([{"matched_compounds": len(matched), "spearman_rho": rho,
                   "exact_permutation_p_two_sided": p,
                   "calibration_scope": "compound_level_virtual_to_measured_expression_response"}]).to_csv(
                       TABLES / "GSE303486_LINCS_calibration_summary.csv", index=False)


def analyse_gse303487(model: pd.DataFrame) -> None:
    path = STAGE / "GSM9128017_gene_count.csv.gz"
    if not path.exists():
        urllib.request.urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9128nnn/GSM9128017/suppl/"
            "GSM9128017_gene_count.csv.gz", path)
    data = pd.read_csv(path, compression="gzip")
    counts = data.set_index("Geneid")
    counts.index = counts.index.astype(str).str.upper()
    counts = counts.groupby(level=0).sum()
    score = score_matrix(log_cpm(counts), model)
    score["condition"] = score.sample_id.str.replace(r"_\d+$", "", regex=True).replace({"FR_1": "FR-1"})
    score.to_csv(TABLES / "GSE303487_focused_drug_scores.csv", index=False)
    rows = []
    for compound in ["Rottlerin", "FR-1"]:
        a = score.loc[score.condition.eq("DMSO"), "fixed_weight_response"].to_numpy(float)
        b = score.loc[score.condition.eq(compound), "fixed_weight_response"].to_numpy(float)
        effect, lo, hi, p_t = welch_ci(a, b)
        rows.append({"compound": compound, "control_n": len(a), "treatment_n": len(b),
                     "mean_control_minus_treatment": effect, "ci95_low": lo, "ci95_high": hi,
                     "welch_p": p_t, "exact_label_permutation_p": independent_permutation(a, b)})
    pd.DataFrame(rows).to_csv(TABLES / "GSE303487_focused_drug_effects.csv", index=False)


def analyse_gse246562(model: pd.DataFrame) -> None:
    path = download("GSE246562", "_counts.csv.gz")
    raw = pd.read_csv(path, compression="gzip")
    sample_columns = [c for c in raw.columns if re.match(r"^[NK]\d\s+(8|214)kPa", c.strip())]
    counts = raw.set_index("Probe")[sample_columns]
    counts.index = counts.index.astype(str).str.upper()
    counts = counts.groupby(level=0).sum()
    score = score_matrix(log_cpm(counts), model)
    score["cell_source"] = np.where(score.sample_id.str.strip().str.startswith("K"), "keloid", "normal")
    score["donor_id"] = score.sample_id.str.extract(r"^([NK]\d)")[0]
    score["stiffness_kpa"] = score.sample_id.str.extract(r"(8|214)kPa")[0].astype(int)
    score.to_csv(TABLES / "GSE246562_donor_stiffness_scores.csv", index=False)
    wide = score.pivot(index=["cell_source", "donor_id"], columns="stiffness_kpa",
                       values="fixed_weight_response").reset_index()
    wide["paired_stiff_minus_soft"] = wide[214] - wide[8]
    wide.to_csv(TABLES / "GSE246562_paired_stiffness_changes.csv", index=False)
    a = wide.loc[wide.cell_source.eq("keloid"), "paired_stiff_minus_soft"].to_numpy(float)
    b = wide.loc[wide.cell_source.eq("normal"), "paired_stiff_minus_soft"].to_numpy(float)
    effect, lo, hi, p_t = welch_ci(a, b)
    pd.DataFrame([{"dataset_id": "GSE246562", "endpoint": "disease_by_stiffness_interaction",
                   "keloid_donors": len(a), "normal_donors": len(b), "difference_in_paired_response": effect,
                   "ci95_low": lo, "ci95_high": hi, "welch_p": p_t,
                   "exact_label_permutation_p": independent_permutation(a, b),
                   "keloid_direction_support": float(np.mean(a > 0)),
                   "normal_direction_support": float(np.mean(b > 0))}]).to_csv(
                       TABLES / "GSE246562_interaction_summary.csv", index=False)


def analyse_gse188952(model: pd.DataFrame) -> None:
    path = download("GSE188952", "_FPKM.tsv.gz")
    raw = pd.read_csv(path, sep="\t", compression="gzip")
    gene_col = next(c for c in raw.columns if "gene" in c.lower() or "symbol" in c.lower())
    expr = raw.set_index(gene_col)
    annotation = [c for c in expr.columns if not pd.api.types.is_numeric_dtype(expr[c])]
    expr = expr.drop(columns=annotation, errors="ignore").apply(pd.to_numeric, errors="coerce").fillna(0)
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.groupby(level=0).median()
    score = score_matrix(np.log2(expr + 0.1), model)
    condition = {}
    for sample in score.sample_id:
        key = sample.lower()
        if "hypertrophic" in key:
            condition[sample] = "hypertrophic_scar"
        elif "keloid" in key:
            condition[sample] = "keloid"
        elif "normotrophic" in key:
            condition[sample] = "normal_scar"
        else:
            condition[sample] = "unresolved"
    score["condition"] = score.sample_id.map(condition)
    score.to_csv(TABLES / "GSE188952_scar_tissue_scores.csv", index=False)
    rows = []
    for positive, reference in [("keloid", "normal_scar"), ("hypertrophic_scar", "normal_scar"),
                                ("keloid", "hypertrophic_scar")]:
        a = score.loc[score.condition.eq(positive), "fixed_weight_response"].to_numpy(float)
        b = score.loc[score.condition.eq(reference), "fixed_weight_response"].to_numpy(float)
        effect, lo, hi, p_t = welch_ci(a, b)
        rows.append({"contrast": f"{positive}_vs_{reference}", "n_positive": len(a), "n_reference": len(b),
                     "mean_difference": effect, "ci95_low": lo, "ci95_high": hi, "welch_p": p_t,
                     "exact_permutation_p": independent_permutation(a, b)})
    pd.DataFrame(rows).to_csv(TABLES / "GSE188952_scar_tissue_effects.csv", index=False)


def analyse_gse191067_index_pairing_sensitivity() -> None:
    """Secondary sensitivity only: sample suffixes suggest, but do not prove, donor pairing."""
    scores = pd.read_csv(TABLES / "GSE191067_frozen_sample_scores_sensitivity.csv")
    scores = scores.query("threshold == 0.5 and annotation_sensitivity == 'reference_transfer'").copy()
    subset = scores[scores.condition.isin(["keloid", "perilesional_skin"])].copy()
    subset["index_id"] = subset.sample_id.str.extract(r"(\d+)$")[0]
    wide = subset.pivot(index="index_id", columns="condition", values="fssi_frozen_weighted").dropna()
    wide["index_paired_difference"] = wide["keloid"] - wide["perilesional_skin"]
    values = wide["index_paired_difference"].to_numpy(float)
    mean = float(values.mean())
    if len(values) > 1:
        sem = stats.sem(values)
        lo, hi = stats.t.interval(0.95, len(values) - 1, loc=mean, scale=sem)
    else:
        lo = hi = np.nan
    signs = np.array(list(itertools.product([-1, 1], repeat=len(values))), dtype=float)
    null = np.mean(signs * values, axis=1)
    p_exact = float(np.mean(np.abs(null) >= abs(mean)))
    wide.reset_index().to_csv(TABLES / "GSE191067_index_paired_sensitivity.csv", index=False)
    pd.DataFrame([{
        "dataset_id": "GSE191067", "contrast": "keloid_minus_perilesional_skin",
        "pairing_status": "sample_suffix_assumption_not_confirmed_by_metadata",
        "pair_n": len(values), "mean_index_paired_difference": mean,
        "ci95_low": lo, "ci95_high": hi, "exact_signflip_p": p_exact,
    }]).to_csv(TABLES / "GSE191067_index_paired_sensitivity_summary.csv", index=False)


def main() -> None:
    model = pd.read_csv(TABLES / "fssi_frozen_model.csv")
    entrez = gene_id_map()
    analyse_gse178411(model, entrez)
    analyse_gse113619(model, entrez)
    analyse_gse303486(model)
    analyse_gse303487(model)
    analyse_gse246562(model)
    analyse_gse188952(model)
    analyse_gse191067_index_pairing_sensitivity()
    print("Completed final release public expansion analyses")


if __name__ == "__main__":
    main()
