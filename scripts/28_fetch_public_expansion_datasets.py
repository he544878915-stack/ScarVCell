from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from pathlib import Path


ROOT = Path(os.environ.get("SCARVCELL_ROOT", Path(__file__).resolve().parents[1]))
STAGE = ROOT / ".stage" / "public_expansion"
STAGE.mkdir(parents=True, exist_ok=True)

FILES = {
    "GSE90051_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE90nnn/GSE90051/matrix/GSE90051_series_matrix.txt.gz",
    "GPL6480.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",
    "GSE212954_Expression_Gene.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/GSE212954/suppl/GSE212954_Expression_Gene.xlsx",
    "GSE151464_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151464/suppl/GSE151464_RAW.tar",
    "GSE83286_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83286/matrix/GSE83286_series_matrix.txt.gz",
    "GPL19612.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL19nnn/GPL19612/annot/GPL19612.annot.gz",
    "GSE282479_VitD_counts.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282479/suppl/GSE282479_VitD_counts.csv.gz",
    "GSE145725_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145725/matrix/GSE145725_series_matrix.txt.gz",
    "GPL16043.annot.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL16nnn/GPL16043/annot/GPL16043.annot.gz",
    "GSE173900_gene_count_matrix_9samples.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173900/suppl/GSE173900_gene_count_matrix_9samples.csv.gz",
    "GSE218007_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE218nnn/GSE218007/matrix/GSE218007_series_matrix.txt.gz",
    "GSE232079_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232079/suppl/GSE232079_counts.txt.gz",
    "GSE232079_varianceStabilizingTransformation.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232079/suppl/GSE232079_varianceStabilizingTransformation.txt.gz",
}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(name: str, url: str) -> dict[str, object]:
    destination = STAGE / name
    status = "existing"
    if not destination.exists() or destination.stat().st_size == 0:
        temporary = destination.with_suffix(destination.suffix + ".part")
        request = urllib.request.Request(url, headers={"User-Agent": "ScarVCell/16 public-data audit"})
        try:
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
                while block := response.read(1024 * 1024):
                    output.write(block)
            temporary.replace(destination)
            status = "downloaded"
        except Exception as exc:
            temporary.unlink(missing_ok=True)
            return {"file": name, "url": url, "status": "failed", "error": str(exc)}
    return {
        "file": name,
        "url": url,
        "status": status,
        "bytes": destination.stat().st_size,
        "md5": md5(destination),
    }


def main() -> None:
    records = [download(name, url) for name, url in FILES.items()]
    manifest = STAGE / "download_manifest.json"
    manifest.write_text(json.dumps(records, indent=2), encoding="utf-8")
    print(json.dumps(records, indent=2))
    failed = [record for record in records if record["status"] == "failed"]
    if failed:
        raise SystemExit(f"{len(failed)} public files failed; see {manifest}")


if __name__ == "__main__":
    main()
