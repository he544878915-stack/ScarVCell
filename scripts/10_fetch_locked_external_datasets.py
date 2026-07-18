"""Fetch the post-lock datasets locked in the locked dataset register before matrix inspection."""

from __future__ import annotations

import hashlib
import os
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE = Path(os.environ.get("SCARVCELL_STAGE", str(ROOT / ".stage"))) / "postlock_locked"
TABLES = ROOT / "results" / "tables"

FILES = {
    "GSE282885": [
        ("GSM8652165", "N1-B", "187867"),
        ("GSM8652166", "N1-Z", "189239"),
        ("GSM8652167", "N4", "194909"),
        ("GSM8652168", "K1", "194893"),
        ("GSM8652169", "K2", "194894"),
        ("GSM8652170", "K3", "194907"),
    ],
    "GSE335482": [
        ("GSM9814425", "H2", "234893"),
        ("GSM9814426", "I", "234900"),
        ("GSM9814427", "H1-2", "267150"),
        ("GSM9814428", "I-2", "267157"),
        ("GSM9814429", "H1_2-20240815", "268545"),
        ("GSM9814430", "I1-20240815", "268542"),
        ("GSM9814431", "H1_2-20241010", "279108"),
        ("GSM9814432", "I1-20241010", "279107"),
    ],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(url) as response, partial.open("wb") as handle:
        while block := response.read(1024 * 1024):
            handle.write(block)
    partial.replace(path)


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        base = destination.resolve()
        for member in tar.getmembers():
            target = (destination / member.name).resolve()
            if base not in target.parents and target != base:
                raise ValueError(f"Unsafe archive member: {member.name}")
        tar.extractall(destination, filter="data")


def main() -> None:
    rows: list[dict[str, object]] = []
    for accession, items in FILES.items():
        sample_dir = STAGE / accession
        for gsm, label, run in items:
            filename = f"{gsm}_{label}_EmptyDrops_CR_{run}_matrix.tar.gz"
            prefix = gsm[:7] + "nnn"
            url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix}/{gsm}/suppl/{filename}"
            archive = sample_dir / "archives" / filename
            download(url, archive)
            extracted = sample_dir / "matrices" / gsm
            safe_extract(archive, extracted)
            rows.append({
                "accession": accession,
                "sample_id": gsm,
                "url": url,
                "local_file": archive.relative_to(STAGE).as_posix(),
                "bytes": archive.stat().st_size,
                "sha256": sha256(archive),
                "status": "downloaded_and_extracted",
            })

    accession = "GSE307210"
    filename = "GSE307210_gene_fpkm.txt.gz"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/{accession}/suppl/{filename}"
    path = STAGE / accession / filename
    download(url, path)
    rows.append({
        "accession": accession,
        "sample_id": "series_matrix",
        "url": url,
        "local_file": path.relative_to(STAGE).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "status": "downloaded",
    })

    manifest = pd.DataFrame(rows)
    TABLES.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(TABLES / "postlock_locked_download_manifest.csv", index=False)
    print(manifest.groupby("accession").agg(files=("sample_id", "size"), bytes=("bytes", "sum")))


if __name__ == "__main__":
    main()
