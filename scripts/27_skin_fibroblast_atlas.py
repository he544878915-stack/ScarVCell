from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path

import fsspec
import numpy as np
import pandas as pd
import zarr
from scipy import stats


ROOT = Path(os.environ.get("SCARVCELL_ROOT", Path(__file__).resolve().parents[1]))
SOURCE = ROOT / "data" / "external" / "steele2025_skin_fibroblast_atlas"
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

ZARR_URL = "https://storage.googleapis.com/haniffalab/skin-fibroblast/zarr/adata_webportal.zarr"
MODEL_PATH = ROOT / "results" / "tables" / "fssi_frozen_model.csv"
AUDIT_PATH = ROOT / "config" / "dataset_screening_audit.csv"
GROUP_COLS = [
    "accession",
    "Patient_status",
    "disease_category_orig",
    "lesional_vs_nonlesional",
    "celltype",
]


def decode_categorical(root: zarr.Group, key: str) -> np.ndarray:
    categories = np.array(list(root[f"obs/{key}/categories"][:]), dtype=object).astype(str)
    codes = np.asarray(root[f"obs/{key}/codes"][:], dtype=int)
    result = np.full(codes.shape, "missing", dtype=object)
    valid = codes >= 0
    result[valid] = categories[codes[valid]]
    return result.astype(str)


def t_interval(values: pd.Series) -> tuple[float, float]:
    x = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if len(x) < 2:
        return math.nan, math.nan
    mean = float(np.mean(x))
    se = float(stats.sem(x))
    margin = float(stats.t.ppf(0.975, len(x) - 1) * se)
    return mean - margin, mean + margin


def bh_adjust(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if not len(pv):
        return pd.Series(out, index=p_values.index)
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[np.where(valid)[0]] = restored
    return pd.Series(out, index=p_values.index)


def accession_name(raw: str) -> str:
    aliases = {
        "Sole-Boldo": "GSE130973",
        "Ganier": "GANIER_SPATIAL_SKIN_ATLAS",
        "Reynolds": "REYNOLDS_HCA_SKIN_PORTAL",
        "GSE": "UNRESOLVED_GSE_LABEL",
    }
    return aliases.get(raw, raw)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    model = pd.read_csv(MODEL_PATH)
    model["gene"] = model.gene.str.upper()
    fixed_genes = model.gene.tolist()

    root = zarr.open_consolidated(fsspec.get_mapper(ZARR_URL), mode="r")
    gene_symbols = np.array(list(root["var/gene_symbol"][:]), dtype=object).astype(str)
    gene_lookup = {gene.upper(): idx for idx, gene in enumerate(gene_symbols)}
    missing = [gene for gene in fixed_genes if gene not in gene_lookup]
    if missing:
        raise RuntimeError(f"Missing frozen genes in atlas: {missing}")

    indices = [gene_lookup[gene] for gene in fixed_genes]
    expression = np.asarray(root["X"].oindex[:, indices], dtype=np.float32)
    if expression.shape != (357276, len(fixed_genes)):
        raise RuntimeError(f"Unexpected selected matrix shape: {expression.shape}")

    means = model.set_index("gene").loc[fixed_genes, "training_mean"].to_numpy(float)
    sds = model.set_index("gene").loc[fixed_genes, "training_sd"].to_numpy(float)
    weights = model.set_index("gene").loc[fixed_genes, "weight"].to_numpy(float)
    equal_weights = model.set_index("gene").loc[fixed_genes, "equal_weight"].to_numpy(float)
    z_expression = (expression - means) / sds
    fssi = z_expression @ weights
    equal_weight_fssi = z_expression @ equal_weights
    pathological_mask = model.set_index("gene").loc[fixed_genes, "program"].eq("pathological").to_numpy()
    repair_mask = ~pathological_mask
    pathological_component = (
        z_expression[:, pathological_mask] @ weights[pathological_mask]
    ) / np.abs(weights[pathological_mask]).sum()
    repair_component = -(
        z_expression[:, repair_mask] @ weights[repair_mask]
    ) / np.abs(weights[repair_mask]).sum()

    obs = pd.DataFrame(
        {
            "GSE_raw": decode_categorical(root, "GSE"),
            "Patient_status": decode_categorical(root, "Patient_status"),
            "disease_category_orig": decode_categorical(root, "disease_category_orig"),
            "lesional_vs_nonlesional": decode_categorical(root, "lesional_vs_nonlesional"),
            "celltype": decode_categorical(root, "celltype"),
            "FSSI": fssi,
            "equal_weight_FSSI": equal_weight_fssi,
            "pathological_component": pathological_component,
            "repair_component": repair_component,
        }
    )
    obs["accession"] = obs.GSE_raw.map(accession_name)
    audit_accessions = set(pd.read_csv(AUDIT_PATH).accession.astype(str))
    obs["scarvcell_overlap"] = obs.accession.isin(audit_accessions)
    obs["accession_auditable"] = obs.accession.ne("UNRESOLVED_GSE_LABEL")

    group = (
        obs.groupby(GROUP_COLS, observed=True)
        .agg(
            n_cells=("FSSI", "size"),
            FSSI=("FSSI", "mean"),
            equal_weight_FSSI=("equal_weight_FSSI", "mean"),
            pathological_component=("pathological_component", "mean"),
            repair_component=("repair_component", "mean"),
            scarvcell_overlap=("scarvcell_overlap", "first"),
            accession_auditable=("accession_auditable", "first"),
        )
        .reset_index()
    )
    group.to_csv(TABLES / "skin_atlas_context_subtype_scores.csv", index=False)

    accession_subtype = (
        group.groupby(["accession", "celltype", "scarvcell_overlap", "accession_auditable"], observed=True)
        .apply(
            lambda x: pd.Series(
                {
                    "n_contexts": len(x),
                    "n_cells": int(x.n_cells.sum()),
                    "FSSI": np.average(x.FSSI, weights=x.n_cells),
                    "equal_weight_FSSI": np.average(x.equal_weight_FSSI, weights=x.n_cells),
                    "pathological_component": np.average(x.pathological_component, weights=x.n_cells),
                    "repair_component": np.average(x.repair_component, weights=x.n_cells),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )
    accession_subtype.to_csv(TABLES / "skin_atlas_accession_subtype_scores.csv", index=False)

    z_columns = [f"z_{gene}" for gene in fixed_genes]
    z_frame = pd.DataFrame(z_expression, columns=z_columns)
    z_frame["accession"] = obs.accession.to_numpy()
    z_frame["celltype"] = obs.celltype.to_numpy()
    z_frame["scarvcell_overlap"] = obs.scarvcell_overlap.to_numpy()
    z_frame["accession_auditable"] = obs.accession_auditable.to_numpy()
    accession_subtype_z = (
        z_frame.groupby(
            ["accession", "celltype", "scarvcell_overlap", "accession_auditable"],
            observed=True,
        )[z_columns]
        .mean()
        .reset_index()
    )
    loo_rows = []
    for gene_index, gene in enumerate(fixed_genes):
        denominator = float(np.abs(weights).sum() - abs(weights[gene_index]))
        loo = accession_subtype_z[
            ["accession", "celltype", "scarvcell_overlap", "accession_auditable"]
        ].copy()
        full_from_z = accession_subtype_z[z_columns].to_numpy(float) @ weights
        loo["score"] = (
            full_from_z
            - weights[gene_index] * accession_subtype_z[f"z_{gene}"].to_numpy(float)
        ) / denominator
        loo = loo[loo.accession_auditable & ~loo.scarvcell_overlap]
        wide = loo.pivot_table(index="accession", columns="celltype", values="score")
        for comparator in ["F2: Universal", "F6: Myofibroblast"]:
            required = {"F7: Fascia-like myofibroblast", comparator}
            if not required.issubset(wide.columns):
                continue
            difference = (wide["F7: Fascia-like myofibroblast"] - wide[comparator]).dropna()
            low, high = t_interval(difference)
            loo_rows.append(
                {
                    "excluded_gene": gene,
                    "contrast": f"F7 minus {comparator}",
                    "n_accessions": len(difference),
                    "mean_difference": difference.mean(),
                    "ci_low": low,
                    "ci_high": high,
                    "positive_accessions": int((difference > 0).sum()),
                }
            )
    pd.DataFrame(loo_rows).to_csv(
        TABLES / "skin_atlas_leave_one_gene_out_subtype_contrasts.csv", index=False
    )

    subtype_rows = []
    for scope, frame in {
        "all_atlas_accessions": accession_subtype[accession_subtype.accession_auditable],
        "nonoverlap_accessions": accession_subtype[
            accession_subtype.accession_auditable & ~accession_subtype.scarvcell_overlap
        ],
    }.items():
        for subtype, sub in frame.groupby("celltype", observed=True):
            low, high = t_interval(sub.FSSI)
            subtype_rows.append(
                {
                    "scope": scope,
                    "celltype": subtype,
                    "n_accessions": sub.accession.nunique(),
                    "mean_FSSI": sub.FSSI.mean(),
                    "mean_equal_weight_FSSI": sub.equal_weight_FSSI.mean(),
                    "ci_low": low,
                    "ci_high": high,
                    "mean_pathological_component": sub.pathological_component.mean(),
                    "mean_repair_component": sub.repair_component.mean(),
                }
            )
    subtype_summary = pd.DataFrame(subtype_rows)
    subtype_summary.to_csv(TABLES / "skin_atlas_consensus_subtype_summary.csv", index=False)

    primary_accession_subtype = accession_subtype[
        accession_subtype.accession_auditable & ~accession_subtype.scarvcell_overlap
    ].copy()
    contrast_rows = []
    for comparator in [
        "F2: Universal",
        "F6: Myofibroblast",
        "F6: Inflammatory myofibroblast",
    ]:
        wide = primary_accession_subtype.pivot_table(
            index="accession", columns="celltype", values="FSSI", aggfunc="mean"
        )
        if {"F7: Fascia-like myofibroblast", comparator}.issubset(wide.columns):
            differences = (
                wide["F7: Fascia-like myofibroblast"] - wide[comparator]
            ).dropna()
            low, high = t_interval(differences)
            positive = int((differences > 0).sum())
            sign_p = float(stats.binomtest(positive, len(differences), 0.5).pvalue)
            contrast_rows.append(
                {
                    "contrast": f"F7 minus {comparator}",
                    "n_accessions": len(differences),
                    "mean_difference": differences.mean(),
                    "ci_low": low,
                    "ci_high": high,
                    "positive_accessions": positive,
                    "exact_sign_p": sign_p,
                }
            )
    pd.DataFrame(contrast_rows).to_csv(
        TABLES / "skin_atlas_consensus_subtype_contrasts.csv", index=False
    )

    context = (
        obs.groupby(
            [
                "accession",
                "Patient_status",
                "disease_category_orig",
                "lesional_vs_nonlesional",
                "scarvcell_overlap",
                "accession_auditable",
            ],
            observed=True,
        )
        .agg(n_cells=("FSSI", "size"), cell_weighted_FSSI=("FSSI", "mean"))
        .reset_index()
    )
    balanced = (
        group.groupby(
            [
                "accession",
                "Patient_status",
                "disease_category_orig",
                "lesional_vs_nonlesional",
                "scarvcell_overlap",
                "accession_auditable",
            ],
            observed=True,
        )
        .agg(n_subtypes=("celltype", "nunique"), subtype_balanced_FSSI=("FSSI", "mean"))
        .reset_index()
    )
    context = context.merge(
        balanced,
        on=[
            "accession",
            "Patient_status",
            "disease_category_orig",
            "lesional_vs_nonlesional",
            "scarvcell_overlap",
            "accession_auditable",
        ],
        how="left",
    )
    context.to_csv(TABLES / "skin_atlas_disease_context_scores.csv", index=False)

    primary_context = context[
        context.accession_auditable & ~context.scarvcell_overlap
    ].copy()
    disease_rows = []
    for disease, sub in primary_context.groupby("Patient_status", observed=True):
        low, high = t_interval(sub.subtype_balanced_FSSI)
        disease_rows.append(
            {
                "Patient_status": disease,
                "n_accessions": sub.accession.nunique(),
                "n_contexts": len(sub),
                "mean_subtype_balanced_FSSI": sub.subtype_balanced_FSSI.mean(),
                "ci_low": low,
                "ci_high": high,
                "mean_cell_weighted_FSSI": sub.cell_weighted_FSSI.mean(),
            }
        )
    disease_summary = pd.DataFrame(disease_rows).sort_values(
        "mean_subtype_balanced_FSSI", ascending=False
    )
    disease_summary.to_csv(TABLES / "skin_atlas_nonoverlap_disease_summary.csv", index=False)

    category_rows = []
    for category, sub in primary_context.groupby("disease_category_orig", observed=True):
        accession_values = (
            sub.groupby("accession", observed=True).subtype_balanced_FSSI.mean()
        )
        low, high = t_interval(accession_values)
        category_rows.append(
            {
                "disease_category_orig": category,
                "n_accessions": len(accession_values),
                "mean_subtype_balanced_FSSI": accession_values.mean(),
                "ci_low": low,
                "ci_high": high,
            }
        )
    pd.DataFrame(category_rows).sort_values(
        "mean_subtype_balanced_FSSI", ascending=False
    ).to_csv(TABLES / "skin_atlas_nonoverlap_category_summary.csv", index=False)

    paired_subtype = group.pivot_table(
        index=["accession", "celltype", "scarvcell_overlap", "accession_auditable"],
        columns="lesional_vs_nonlesional",
        values="FSSI",
        aggfunc="mean",
    ).reset_index()
    if {"Lesional", "Nonlesional"}.issubset(paired_subtype.columns):
        paired_subtype = paired_subtype.dropna(subset=["Lesional", "Nonlesional"]).copy()
        paired_subtype["lesional_minus_nonlesional"] = (
            paired_subtype.Lesional - paired_subtype.Nonlesional
        )
        paired = (
            paired_subtype.groupby(
                ["accession", "scarvcell_overlap", "accession_auditable"], observed=True
            )
            .agg(
                matched_subtypes=("celltype", "nunique"),
                lesional_minus_nonlesional=("lesional_minus_nonlesional", "mean"),
            )
            .reset_index()
        )
    else:
        paired = paired_subtype.iloc[0:0].copy()
        paired["lesional_minus_nonlesional"] = []
    paired_subtype.to_csv(
        TABLES / "skin_atlas_lesional_within_accession_subtype_effects.csv", index=False
    )
    paired.to_csv(TABLES / "skin_atlas_lesional_within_accession_effects.csv", index=False)

    published_ranks = pd.read_csv(SOURCE / "degs_fbs.csv", index_col=0)
    published_ranks = published_ranks.apply(
        lambda column: column.astype("string").str.upper()
    )
    fixed_sets = {
        "pathological": set(model.loc[model.program.eq("pathological"), "gene"]),
        "repair": set(model.loc[model.program.eq("repair"), "gene"]),
    }
    enrichment_rows = []
    rng = np.random.default_rng(150713)
    for subtype in published_ranks.columns:
        ranked_genes = published_ranks[subtype].dropna().drop_duplicates().tolist()
        rank_lookup = {gene: rank for rank, gene in enumerate(ranked_genes, start=1)}
        universe = set(ranked_genes)
        markers = set(ranked_genes[:100])
        background_values = np.arange(1, len(ranked_genes) + 1, dtype=float)
        for programme, genes in fixed_sets.items():
            genes = genes & universe
            overlap = markers & genes
            observed_genes = sorted(genes)
            observed_mean_rank = float(np.mean([rank_lookup[gene] for gene in observed_genes]))
            observed_mean_percentile = float(
                1 - (observed_mean_rank - 1) / max(len(ranked_genes) - 1, 1)
            )
            null = np.empty(10000)
            for iteration in range(len(null)):
                random_mean_rank = rng.choice(
                    background_values, size=len(observed_genes), replace=False
                ).mean()
                null[iteration] = 1 - (random_mean_rank - 1) / max(
                    len(ranked_genes) - 1, 1
                )
            p_value = float(
                (1 + np.sum(null >= observed_mean_percentile)) / (len(null) + 1)
            )
            enrichment_rows.append(
                {
                    "celltype": subtype,
                    "programme": programme,
                    "marker_genes": len(markers),
                    "fixed_genes_in_universe": len(observed_genes),
                    "overlap_n": len(overlap),
                    "overlap_genes": ";".join(sorted(overlap)),
                    "mean_published_rank": observed_mean_rank,
                    "mean_rank_percentile": observed_mean_percentile,
                    "random_percentile_mean": float(null.mean()),
                    "random_percentile_sd": float(null.std(ddof=1)),
                    "empirical_p": p_value,
                }
            )
    enrichment = pd.DataFrame(enrichment_rows)
    enrichment["fdr"] = bh_adjust(enrichment.empirical_p)
    enrichment.to_csv(TABLES / "skin_atlas_published_marker_crosswalk.csv", index=False)

    file_rows = []
    for path in sorted(SOURCE.iterdir()):
        if path.is_file():
            file_rows.append(
                {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
            )
    pd.DataFrame(file_rows).to_csv(
        TABLES / "skin_atlas_public_input_manifest.csv", index=False
    )

    nonoverlap_subtypes = subtype_summary[subtype_summary.scope.eq("nonoverlap_accessions")]
    top_subtype = nonoverlap_subtypes.sort_values("mean_FSSI", ascending=False).iloc[0]
    top_disease = disease_summary.iloc[0] if len(disease_summary) else None
    primary_paired = paired[paired.accession_auditable & ~paired.scarvcell_overlap]
    paired_low, paired_high = t_interval(primary_paired.lesional_minus_nonlesional)
    paired_positive = int((primary_paired.lesional_minus_nonlesional > 0).sum())
    paired_sign_p = (
        float(stats.binomtest(paired_positive, len(primary_paired), 0.5).pvalue)
        if len(primary_paired)
        else None
    )
    summary = {
        "source_cells": int(len(obs)),
        "source_genes": int(len(gene_symbols)),
        "frozen_genes_detected": int(len(fixed_genes) - len(missing)),
        "atlas_accessions": int(obs.accession.nunique()),
        "overlapping_accessions": sorted(obs.loc[obs.scarvcell_overlap, "accession"].unique().tolist()),
        "nonoverlap_auditable_accessions": int(
            obs.loc[obs.accession_auditable & ~obs.scarvcell_overlap, "accession"].nunique()
        ),
        "top_nonoverlap_consensus_subtype": str(top_subtype.celltype),
        "top_nonoverlap_consensus_subtype_mean_fssi": float(top_subtype.mean_FSSI),
        "top_nonoverlap_disease": None if top_disease is None else str(top_disease.Patient_status),
        "top_nonoverlap_disease_mean_fssi": None
        if top_disease is None
        else float(top_disease.mean_subtype_balanced_FSSI),
        "nonoverlap_lesional_effects": int(len(primary_paired)),
        "nonoverlap_lesional_positive": paired_positive,
        "nonoverlap_lesional_mean_difference": None
        if not len(primary_paired)
        else float(primary_paired.lesional_minus_nonlesional.mean()),
        "nonoverlap_lesional_ci_low": None if not len(primary_paired) else paired_low,
        "nonoverlap_lesional_ci_high": None if not len(primary_paired) else paired_high,
        "nonoverlap_lesional_exact_sign_p": paired_sign_p,
        "expression_min": float(np.min(expression)),
        "expression_max": float(np.max(expression)),
        "expression_p99": float(np.quantile(expression, 0.99)),
    }
    (TABLES / "skin_atlas_main_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
