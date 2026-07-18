"""Prepare ancestry-specific GWAS inputs, run MAGMA and summarise locked gene sets."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(os.environ.get("SCARVCELL_ROOT", Path(__file__).resolve().parents[1]))
STAGE = ROOT / ".stage" / "genetics" / "MAGMA"
GWAS = ROOT / ".stage" / "genetics" / "GWAS"
ANALYSIS = STAGE / "analysis"
TABLES = ROOT / "results" / "tables"
MAGMA = STAGE / "program" / "magma.exe"
GENE_LOC = STAGE / "gene_locations" / "NCBI37.3.gene.loc"
CONFIG = {
    "EUR": {"accession": "GCST90652488", "n": 1282582,
            "bfile": STAGE / "reference_eur" / "g1000_eur"},
    "AFR": {"accession": "GCST90652489", "n": 139538,
            "bfile": STAGE / "reference_afr" / "g1000_afr",
            "overlap_bfile": STAGE / "reference_afr" / "g1000_afr_gwas_overlap"},
}
SEED = 20260712
N_NULL = 10000


def prepare(ancestry: str) -> None:
    cfg = CONFIG[ancestry]
    source = GWAS / f"{cfg['accession']}.tsv.gz"
    pval = ANALYSIS / f"{ancestry.lower()}_gwas_pval.txt"
    loc = ANALYSIS / f"{ancestry.lower()}_gwas_snp_loc.txt"
    counts = {"source_rows": 0, "retained_rows": 0, "missing_rsid_or_p": 0,
              "invalid_p": 0, "mhc_excluded": 0}
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rt", newline="") as handle, pval.open("w", newline="") as p_out, loc.open("w", newline="") as l_out:
        reader = csv.DictReader(handle, delimiter="\t")
        p_writer = csv.writer(p_out, delimiter="\t", lineterminator="\n")
        l_writer = csv.writer(l_out, delimiter="\t", lineterminator="\n")
        p_writer.writerow(["SNP", "P"]); l_writer.writerow(["SNP", "CHR", "BP"])
        for row in reader:
            counts["source_rows"] += 1
            rsid, p = row.get("rsid", ""), row.get("p_value", "")
            if not rsid or rsid == "NA" or not p or p == "NA":
                counts["missing_rsid_or_p"] += 1; continue
            try:
                chrom = int(row["chromosome"]); bp = int(row["base_pair_location"]); pv = float(p)
            except (ValueError, TypeError):
                counts["invalid_p"] += 1; continue
            if not 0 < pv <= 1:
                counts["invalid_p"] += 1; continue
            if chrom == 6 and 25_000_000 <= bp <= 34_000_000:
                counts["mhc_excluded"] += 1; continue
            p_writer.writerow([rsid, f"{pv:.16g}"])
            l_writer.writerow([rsid, chrom, bp])
            counts["retained_rows"] += 1
    counts |= {"ancestry": ancestry, "accession": cfg["accession"], "sample_size": cfg["n"],
               "genome_build": "GRCh37", "gene_window_kb": 10}
    (ANALYSIS / f"{ancestry.lower()}_input_audit.json").write_text(json.dumps(counts, indent=2), encoding="utf-8")
    print(json.dumps(counts, indent=2))


def run_command(command: list[str], log: Path) -> None:
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    log.write_text("COMMAND\n" + subprocess.list2cmdline(command) + "\n\nSTDOUT\n" + result.stdout +
                   "\nSTDERR\n" + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"MAGMA failed ({result.returncode}); see {log}")


def run_magma(ancestry: str) -> None:
    cfg = CONFIG[ancestry]; low = ancestry.lower()
    bfile = cfg.get("overlap_bfile", cfg["bfile"])
    if not Path(str(bfile) + ".bed").exists():
        bfile = cfg["bfile"]
    annotation = ANALYSIS / f"{low}_10kb"
    result = ANALYSIS / f"{low}_keloid"
    run_command([str(MAGMA), "--annotate", "window=10", "--snp-loc", str(ANALYSIS / f"{low}_gwas_snp_loc.txt"),
                 "--gene-loc", str(GENE_LOC), "--out", str(annotation)], ANALYSIS / f"{low}_annotate_command.log")
    run_command([str(MAGMA), "--bfile", str(bfile), "--pval", str(ANALYSIS / f"{low}_gwas_pval.txt"),
                 f"N={cfg['n']}", "duplicate=error", "--gene-annot", str(annotation) + ".genes.annot",
                 "--out", str(result)], ANALYSIS / f"{low}_gene_command.log")
    run_command([str(MAGMA), "--gene-results", str(result) + ".genes.raw", "--set-annot",
                 str(ANALYSIS / "locked_gene_sets_magma.txt"), "--out", str(result) + "_sets"],
                ANALYSIS / f"{low}_set_command.log")
    print(f"Completed MAGMA for {ancestry}")


def read_gene_results(ancestry: str) -> pd.DataFrame:
    path = ANALYSIS / f"{ancestry.lower()}_keloid.genes.out.txt"
    data = pd.read_csv(path, sep=r"\s+")
    data["GENE"] = data["GENE"].astype(str)
    locations = pd.read_csv(GENE_LOC, sep=r"\s+", header=None,
                            names=["GENE", "CHR_LOC", "START", "STOP", "STRAND", "SYMBOL"], dtype={"GENE": str})
    data = data.merge(locations[["GENE", "SYMBOL"]], on="GENE", how="left")
    data["ancestry"] = ancestry
    data["bonferroni_p"] = np.minimum(data["P"] * len(data), 1)
    return data


def matched_null(data: pd.DataFrame, membership: pd.DataFrame, ancestry: str) -> pd.DataFrame:
    expression = pd.read_csv(TABLES / "scrna_broad_celltype_mean_logcpm_signature.csv", index_col=0)
    expression.index = expression.index.astype(str).str.upper()
    universe = data.dropna(subset=["SYMBOL", "NSNPS", "ZSTAT"]).copy()
    universe["SYMBOL"] = universe["SYMBOL"].str.upper()
    universe = universe.drop_duplicates("SYMBOL").set_index("SYMBOL")
    universe["fibroblast_expression"] = expression.reindex(universe.index)["fibroblast"]
    universe = universe.dropna(subset=["fibroblast_expression"])
    universe["snp_bin"] = pd.qcut(np.log1p(universe.NSNPS), 5, labels=False, duplicates="drop")
    universe["expression_bin"] = pd.qcut(universe.fibroblast_expression.rank(method="first"), 5, labels=False)
    rng = np.random.default_rng(SEED + (0 if ancestry == "EUR" else 1))
    all_symbols = universe.index.to_numpy()
    rows = []
    for set_name, block in membership.groupby("set"):
        genes = [g for g in block.symbol.str.upper() if g in universe.index]
        observed = float(universe.loc[genes, "ZSTAT"].mean()) if genes else np.nan
        pools = []
        excluded = set(genes)
        for gene in genes:
            row = universe.loc[gene]
            pool = universe.index[(universe.snp_bin == row.snp_bin) &
                                  (universe.expression_bin == row.expression_bin) &
                                  ~universe.index.isin(excluded)].to_numpy()
            if not len(pool):
                pool = np.asarray([symbol for symbol in all_symbols if symbol not in excluded])
            pools.append(pool)
        null = np.empty(N_NULL, dtype=float)
        for _ in range(N_NULL):
            selected = []
            used = set(excluded)
            for pool in pools:
                choice = str(rng.choice(pool))
                while choice in used:
                    choice = str(rng.choice(pool))
                selected.append(choice)
                used.add(choice)
            null[_] = float(universe.loc[selected, "ZSTAT"].mean()) if selected else np.nan
        p = (1 + np.sum(null >= observed)) / (1 + np.sum(np.isfinite(null))) if np.isfinite(observed) else np.nan
        rows.append({"ancestry": ancestry, "set": set_name, "genes_requested": len(block),
                     "genes_evaluable": len(genes), "observed_mean_z": observed,
                     "matched_null_mean": float(np.nanmean(null)), "matched_null_sd": float(np.nanstd(null, ddof=1)),
                     "matched_empirical_p_one_sided": p, "iterations": N_NULL})
    return pd.DataFrame(rows)


def summarise() -> None:
    mapping = pd.read_csv(ANALYSIS / "locked_symbol_entrez_mapping.csv", dtype={"entrez_id": str})
    membership = pd.read_csv(TABLES / "locked_genetic_gene_sets.csv", dtype={"entrez_id": str})
    genes, nulls, sets = [], [], []
    for ancestry in CONFIG:
        data = read_gene_results(ancestry); genes.append(data)
        candidates = data[data.SYMBOL.str.upper().isin(mapping.symbol.str.upper())].copy()
        candidates = mapping.merge(candidates, left_on="entrez_id", right_on="GENE", how="left")
        candidates = candidates.drop(columns=["SYMBOL"])
        candidates["ancestry"] = ancestry
        candidates.to_csv(TABLES / f"{ancestry.lower()}_locked_candidate_magma.csv", index=False)
        nulls.append(matched_null(data, membership, ancestry))
        gsa = pd.read_csv(ANALYSIS / f"{ancestry.lower()}_keloid_sets.gsa.out.txt", sep=r"\s+", comment="#")
        gsa["ancestry"] = ancestry; sets.append(gsa)
    all_genes = pd.concat(genes, ignore_index=True)
    all_genes.to_csv(TABLES / "magma_all_gene_results.csv.gz", index=False, compression="gzip")
    pd.concat(nulls, ignore_index=True).to_csv(TABLES / "locked_gene_set_matched_null.csv", index=False)
    pd.concat(sets, ignore_index=True).to_csv(TABLES / "locked_gene_set_magma_competitive.csv", index=False)
    pivot = all_genes.pivot_table(index="GENE", columns="ancestry", values="ZSTAT").dropna()
    rho, p = spearmanr(pivot["EUR"], pivot["AFR"])
    pd.DataFrame([{"genes_shared": len(pivot), "spearman_rho_gene_z": rho, "p_value": p,
                   "scope": "cross_ancestry_gene_level_concordance"}]).to_csv(
                       TABLES / "cross_ancestry_magma_concordance.csv", index=False)
    print("Summarised ancestry-specific MAGMA results")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("step", choices=["prepare", "magma", "summarise", "all"])
    parser.add_argument("--ancestry", choices=list(CONFIG))
    args = parser.parse_args()
    ancestries = [args.ancestry] if args.ancestry else list(CONFIG)
    if args.step in {"prepare", "all"}:
        for ancestry in ancestries: prepare(ancestry)
    if args.step in {"magma", "all"}:
        for ancestry in ancestries: run_magma(ancestry)
    if args.step in {"summarise", "all"}: summarise()


if __name__ == "__main__":
    main()
