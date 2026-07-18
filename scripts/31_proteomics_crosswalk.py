from __future__ import annotations

import hashlib
import io
import json
import math
import re
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion" / "proteomics"
STAGE.mkdir(parents=True, exist_ok=True)
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)
MODEL = pd.read_csv(TABLES / "fssi_frozen_model.csv")
FROZEN_GENES = MODEL.gene.astype(str).str.upper().tolist()


def supplementary_links(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    links = re.findall(r'(?:xlink:)?href="([^"]+)"', text)
    return sorted(
        {
            link
            for link in links
            if re.search(r"supp|mmc|MOESM|xlsx?|csv|media", link, flags=re.IGNORECASE)
        }
    )


def link_contexts(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    contexts = {}
    for match in re.finditer(r'(?:xlink:)?href="([^"]+)"', text):
        link = match.group(1)
        if re.search(r"mmc|xlsx?", link, flags=re.IGNORECASE):
            window = text[max(0, match.start() - 350):min(len(text), match.end() + 350)]
            contexts[link] = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", window)).strip()
    return contexts


def pdf_text(path: Path) -> str:
    with pdfplumber.open(path) as document:
        return "\n".join(page.extract_text() or "" for page in document.pages)


def extract_pxd015057_rows(path: Path, extract: str) -> pd.DataFrame:
    rows = []
    accession_pattern = re.compile(r"^(.+?)\s+([A-Z][0-9][A-Z0-9]{3}[0-9])\s+(\d+(?:\.\d+)?)\s+(.+)$")
    for line in pdf_text(path).splitlines():
        match = accession_pattern.match(line.strip())
        if not match:
            continue
        description, accession, mass, remainder = match.groups()
        values = remainder.split()
        if len(values) != 12 or any(value != "-" and not re.fullmatch(r"\d+(?:\.\d+)?", value) for value in values):
            continue
        row = {"extract": extract, "description": description, "uniprot": accession, "mass_kda": float(mass)}
        columns = []
        for tissue in ["normal_skin", "normal_scar", "keloid"]:
            columns.extend([f"{tissue}_sequence_coverage", f"{tissue}_peptides", f"{tissue}_unique_peptides", f"{tissue}_psm"])
        row.update({column: (math.nan if value == "-" else float(value)) for column, value in zip(columns, values)})
        rows.append(row)
    return pd.DataFrame(rows)


def uniprot_mapping(accessions: list[str]) -> pd.DataFrame:
    frames = []
    for start in range(0, len(accessions), 40):
        batch = accessions[start:start + 40]
        query = " OR ".join(f"accession:{accession}" for accession in batch)
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
            {"query": f"({query})", "format": "tsv", "fields": "accession,gene_primary,protein_name", "size": 500}
        )
        request = urllib.request.Request(url, headers={"User-Agent": "ScarVCell/16 protein crosswalk"})
        with urllib.request.urlopen(request, timeout=120) as response:
            frames.append(pd.read_csv(io.BytesIO(response.read()), sep="\t"))
    output = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return output.rename(columns={"Entry": "uniprot", "Gene Names (primary)": "gene", "Protein names": "protein_name"})


def parse_pxd029631_workbook(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=0, header=None)
    title = str(raw.iloc[0, 0]) if len(raw) else path.stem
    lower_title = title.lower()
    if "up-reg" in lower_title:
        evidence = "source_reported_upregulated"
    elif "down-reg" in lower_title:
        evidence = "source_reported_downregulated"
    elif "uniquely" in lower_title:
        evidence = "source_reported_keloid_unique"
    else:
        evidence = "supporting_table"
    header_index = None
    gene_column = None
    for index, row in raw.iterrows():
        for column, value in row.items():
            if str(value).strip().lower() in {"gene name", "gene_name", "gene"}:
                header_index = index
                gene_column = column
                break
        if header_index is not None:
            break
    if header_index is None or gene_column is None:
        return pd.DataFrame()
    header = [str(value).strip() for value in raw.iloc[header_index].tolist()]
    body = raw.iloc[header_index + 1:].copy()
    body.columns = header
    gene_header = header[gene_column]
    body = body.loc[body[gene_header].notna()].copy()
    body["gene"] = body[gene_header].astype(str).str.upper().str.strip()
    normalised_columns = {str(column).strip().lower().replace("_", " "): column for column in body.columns}
    protein_column = next((column for key, column in normalised_columns.items() if key in {"protein name", "protein_name"}), None)
    fold_column = next((column for key, column in normalised_columns.items() if "fold change" in key), None)
    p_column = next((column for key, column in normalised_columns.items() if "p-value" in key or "p value" in key), None)
    output = pd.DataFrame(
        {
            "gene": body.gene,
            "protein_name": body[protein_column] if protein_column is not None else "",
            "fold_change": pd.to_numeric(body[fold_column], errors="coerce") if fold_column is not None else math.nan,
            "p_value": pd.to_numeric(body[p_column], errors="coerce") if p_column is not None else math.nan,
            "evidence": evidence,
            "source_file": path.name,
        }
    )
    return output


def analyse_proteomics() -> dict[str, object]:
    pxd015 = pd.concat(
        [
            extract_pxd015057_rows(STAGE / "PMC7852214_supplementary" / "mmc2.pdf", "NaCl"),
            extract_pxd015057_rows(STAGE / "PMC7852214_supplementary" / "mmc3.pdf", "GuHCl"),
        ],
        ignore_index=True,
    )
    mapping = uniprot_mapping(sorted(pxd015.uniprot.unique().tolist()))
    pxd015 = pxd015.merge(mapping, on="uniprot", how="left")
    pxd015["gene"] = pxd015.gene.fillna("").astype(str).str.upper().str.split().str[0]
    s5_text = pdf_text(STAGE / "PMC7852214_supplementary" / "mmc5.pdf")
    s5_symbols = {
        match.group(1)
        for line in s5_text.splitlines()
        if (match := re.match(r"^([A-Z][A-Z0-9]+)\s+", line.strip()))
    }
    pxd015["source_keloid_unique_or_increased"] = pxd015.gene.isin(s5_symbols)
    pxd015.to_csv(TABLES / "PXD015057_aggregated_ecm_protein_crosswalk.csv", index=False)

    pxd029_frames = [
        parse_pxd029631_workbook(path)
        for path in sorted((STAGE / "PMC9541363_supplementary").glob("*.xlsx"))
    ]
    pxd029 = pd.concat([frame for frame in pxd029_frames if not frame.empty], ignore_index=True)
    pxd029.to_csv(TABLES / "PXD029631_source_reported_protein_tables.csv", index=False)

    rows = []
    for gene in FROZEN_GENES:
        observed_015 = pxd015.loc[pxd015.gene.eq(gene)]
        observed_029 = pxd029.loc[pxd029.gene.eq(gene)]
        evidence_015 = sorted(observed_015.extract.unique().tolist())
        evidence_029 = sorted(observed_029.evidence.unique().tolist())
        rows.append(
            {
                "gene": gene,
                "programme": MODEL.set_index("gene").loc[gene, "program"],
                "PXD015057_detected_extracts": ";".join(evidence_015),
                "PXD015057_source_keloid_unique_or_increased": bool(observed_015.source_keloid_unique_or_increased.any()) if len(observed_015) else False,
                "PXD029631_source_evidence": ";".join(evidence_029),
                "orthogonal_protein_evidence": bool(evidence_015 or evidence_029),
            }
        )
    crosswalk = pd.DataFrame(rows)
    crosswalk.to_csv(TABLES / "fixed_programme_protein_crosswalk.csv", index=False)
    return {
        "PXD015057_parsed_proteins": int(pxd015.gene.ne("").sum()),
        "PXD015057_fixed_genes_detected": int(crosswalk.PXD015057_detected_extracts.ne("").sum()),
        "PXD015057_fixed_genes_source_classified_keloid_unique_or_increased": int(crosswalk.PXD015057_source_keloid_unique_or_increased.sum()),
        "PXD029631_source_table_rows": int(len(pxd029)),
        "PXD029631_fixed_genes_with_source_evidence": int(crosswalk.PXD029631_source_evidence.ne("").sum()),
        "fixed_genes_with_any_orthogonal_protein_evidence": int(crosswalk.orthogonal_protein_evidence.sum()),
        "scope": "source-level protein crosswalk; released supplements do not expose patient-level abundance matrices",
    }


def main() -> None:
    records = {}
    for pmc in ["PMC7852214", "PMC9541363"]:
        path = STAGE / f"{pmc}.xml"
        records[pmc] = {"links": supplementary_links(path), "contexts": link_contexts(path)}
    (STAGE / "supplementary_link_audit.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    result = analyse_proteomics()
    (TABLES / "proteomics_crosswalk_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
