"""Stream the dense GEO CSV, reproduce source QC and retain locked transfer features."""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path(__file__).resolve().parents[1]
STAGE = Path(os.environ.get("SCARVCELL_STAGE", str(ROOT / ".stage"))) / "external" / "GSE191067"
TABLES = ROOT / "results" / "tables"
SOURCE = STAGE / "GSE191067_all.UMI.matrix.csv.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    selected = pd.read_csv(STAGE / "selected_ensembl_symbol_map.csv", dtype=str)
    all_map = pd.read_csv(STAGE / "all_ensembl_symbol_map.csv", dtype=str)
    desired = set((STAGE / "desired_symbols.txt").read_text(encoding="utf-8").splitlines())
    ensg_to_symbols = selected.groupby("ENSEMBL").SYMBOL.apply(list).to_dict()
    mt_ensg = set(all_map.loc[all_map.SYMBOL.str.startswith("MT-", na=False), "ENSEMBL"])

    with gzip.open(SOURCE, "rt", newline="") as handle:
        header = next(csv.reader(handle))
        barcodes = np.asarray(header[1:], dtype=object)
        n = len(barcodes)
        qc_path = TABLES / "GSE191067_cell_qc_audit.csv.gz"
        reuse_qc = qc_path.exists()
        if reuse_qc:
            metadata = pd.read_csv(qc_path)
            if len(metadata) != n or not np.array_equal(metadata.cell_id.astype(str).to_numpy(), barcodes.astype(str)):
                raise ValueError("Existing QC audit does not match the source-matrix header")
            totals = metadata.nCount_RNA.to_numpy(np.int64)
            detected = metadata.nFeature_RNA.to_numpy(np.int32)
            mt_counts = np.rint(metadata.percent_mt.to_numpy(float) * totals / 100.0).astype(np.int64)
        else:
            totals = np.zeros(n, dtype=np.int64)
            detected = np.zeros(n, dtype=np.int32)
            mt_counts = np.zeros(n, dtype=np.int64)
        rows: dict[str, sparse.csr_matrix] = {}
        genes_seen = 0
        selected_rows = 0
        for line in handle:
            split = line.find(",")
            if split < 1:
                continue
            ensg = line[:split].strip('"').split(".")[0]
            if reuse_qc and ensg not in ensg_to_symbols:
                genes_seen += 1
                if genes_seen % 5000 == 0:
                    print(f"streamed_genes={genes_seen} selected_rows={selected_rows}", flush=True)
                continue
            values = np.fromstring(line[split + 1 :], sep=",", dtype=np.int32)
            if values.size != n:
                raise ValueError(f"Malformed row {ensg}: {values.size} values for {n} cells")
            if not reuse_qc:
                totals += values
                detected += values > 0
                if ensg in mt_ensg:
                    mt_counts += values
            if ensg in ensg_to_symbols:
                row = sparse.csr_matrix(values.reshape(1, -1))
                for symbol in ensg_to_symbols[ensg]:
                    if symbol in desired:
                        rows[symbol] = rows.get(symbol, sparse.csr_matrix((1, n), dtype=np.int32)) + row
                selected_rows += 1
            genes_seen += 1
            if genes_seen % 5000 == 0:
                print(f"streamed_genes={genes_seen} selected_rows={selected_rows}", flush=True)

    sample = np.asarray([x.split("_", 1)[0] for x in barcodes], dtype=object)
    condition_map = {
        "HK1": "keloid", "HK2": "keloid", "HK3": "keloid",
        "HK-NS1": "perilesional_skin", "HK-NS2": "perilesional_skin", "HK-NS3": "perilesional_skin",
        "HNS1": "normal_skin", "HNS2": "normal_skin", "HNS3": "normal_skin",
        "HNSR1": "normal_scar", "HNSR2": "normal_scar", "HNSR3": "normal_scar",
    }
    unknown = sorted(set(sample) - set(condition_map))
    if unknown:
        raise ValueError(f"Unexpected sample prefixes: {unknown}")
    mt_pct = np.divide(mt_counts * 100.0, totals, out=np.zeros(n, float), where=totals > 0)
    keep = metadata.pass_source_qc.to_numpy(bool) if reuse_qc else (
        (totals < 8000) & (detected >= 500) & (detected <= 4000) & (mt_pct < 10)
    )
    if not reuse_qc:
        metadata = pd.DataFrame({
            "cell_id": barcodes, "sample_id": sample,
            "condition": [condition_map[x] for x in sample],
            "nCount_RNA": totals, "nFeature_RNA": detected, "percent_mt": mt_pct,
            "pass_source_qc": keep,
        })
        metadata.to_csv(qc_path, index=False, compression="gzip")
    else:
        metadata["sample_id"] = sample
        metadata["condition"] = [condition_map[x] for x in sample]

    symbols = sorted(rows)
    matrix = sparse.vstack([rows[g][:, keep] for g in symbols], format="csr")
    out10x = STAGE / "locked_feature_matrix"
    out10x.mkdir(parents=True, exist_ok=True)
    mtx = out10x / "matrix.mtx"
    with mtx.open("wb") as handle:
        mmwrite(handle, matrix)
    with mtx.open("rb") as source, gzip.open(str(mtx) + ".gz", "wb", compresslevel=6) as target:
        shutil.copyfileobj(source, target)
    mtx.unlink()
    with gzip.open(out10x / "features.tsv.gz", "wt", newline="") as handle:
        for gene in symbols:
            handle.write(f"{gene}\t{gene}\tGene Expression\n")
    with gzip.open(out10x / "barcodes.tsv.gz", "wt", newline="") as handle:
        for barcode in barcodes[keep]:
            handle.write(f"{barcode}\n")

    per_sample = metadata.groupby(["sample_id", "condition"], as_index=False).agg(
        barcodes_before_qc=("cell_id", "size"), cells_after_source_qc=("pass_source_qc", "sum"),
        median_umi=("nCount_RNA", "median"), median_detected_genes=("nFeature_RNA", "median"),
        median_percent_mt=("percent_mt", "median"),
    )
    per_sample["geo_sample_mapping_status"] = np.where(
        per_sample.sample_id.isin(["HNS1", "HNS2", "HNS3"]),
        "three_normal_skin_matrix_prefixes_recoverable_but_not_one_to_one_mappable_to_GEO_titles_HNS3_HNS4_HNS6",
        "matrix_prefix_matches_GEO_sample_title",
    )
    per_sample.to_csv(TABLES / "GSE191067_source_qc_summary.csv", index=False)
    pd.DataFrame([{
        "accession": "GSE191067", "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE191nnn/GSE191067/suppl/GSE191067_all.UMI.matrix.csv.gz",
        "local_file": SOURCE.relative_to(STAGE.parent).as_posix(), "bytes": SOURCE.stat().st_size,
        "sha256": sha256(SOURCE), "barcodes_before_qc": n, "cells_after_source_qc": int(keep.sum()),
        "genes_streamed": genes_seen, "locked_features_retained": len(symbols),
    }]).to_csv(TABLES / "GSE191067_download_and_conversion_manifest.csv", index=False)
    print(per_sample.to_string(index=False))
    print(f"retained_cells={keep.sum()} retained_features={len(symbols)}")


if __name__ == "__main__":
    main()
