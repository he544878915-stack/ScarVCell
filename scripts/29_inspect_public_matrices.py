from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def inspect_series(path: Path) -> None:
    metadata: list[str] = []
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                metadata.append(line.rstrip())
            if line.startswith("!series_matrix_table_begin"):
                header = next(handle).rstrip()
                first = next(handle).rstrip()
                break
        else:
            header = "missing"
            first = "missing"
    print(f"\n## {path.name}")
    selected = [
        line for line in metadata
        if line.startswith(("!Sample_title", "!Sample_source_name", "!Sample_characteristics", "!Sample_description"))
    ]
    for line in selected + metadata[-5:]:
        print(line[:2000])
    print("HEADER", header[:2000])
    print("FIRST", first[:2000])


def main() -> None:
    raw_83286 = STAGE / "GSM2198196_N1.txt.gz"
    if raw_83286.exists():
        print(f"\n## {raw_83286.name}")
        with gzip.open(raw_83286, "rt", errors="replace") as handle:
            for _ in range(20):
                print(next(handle).rstrip()[:3000])

    for name in [
        "GSE90051_series_matrix.txt.gz",
        "GSE83286_series_matrix.txt.gz",
        "GSE145725_series_matrix.txt.gz",
        "GSE218007_series_matrix.txt.gz",
        "GSE282479_series_matrix.txt.gz",
    ]:
        path = STAGE / name
        if path.exists():
            inspect_series(path)

    for name in [
        "GSE282479_VitD_counts.csv.gz",
        "GSE173900_gene_count_matrix_9samples.csv.gz",
    ]:
        path = STAGE / name
        if path.exists():
            frame = pd.read_csv(path, nrows=5)
            print(f"\n## {name}")
            print(frame.columns.tolist())
            print(frame.to_string(index=False))

    for name in ["GPL6480.annot.gz", "GPL19612.annot.gz", "GPL16043.annot.gz"]:
        path = STAGE / name
        if path.exists():
            with gzip.open(path, "rt", errors="replace") as handle:
                lines = handle.readlines()
            start = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_begin")) + 1
            from io import StringIO
            frame = pd.read_csv(StringIO("".join(lines[start:start + 6])), sep="\t", dtype=str)
            print(f"\n## {name}")
            print(frame.columns.tolist())
            print(frame.to_string(index=False))

    path = STAGE / "GSE212954_Expression_Gene.xlsx"
    if path.exists():
        workbook = pd.ExcelFile(path)
        print(f"\n## {path.name}: {workbook.sheet_names}")
        for sheet in workbook.sheet_names:
            frame = pd.read_excel(workbook, sheet_name=sheet, header=None, nrows=15)
            print(f"\n### {sheet}")
            print(frame.columns.tolist())
            print(frame.to_string(index=False))
            parsed = pd.read_excel(workbook, sheet_name=sheet, header=9, nrows=2)
            print("HEADER9", [repr(column) for column in parsed.columns])



if __name__ == "__main__":
    main()
