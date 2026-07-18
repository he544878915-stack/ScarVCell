from __future__ import annotations

import itertools
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion"
TABLES = ROOT / "results" / "tables"
MODEL = pd.read_csv(TABLES / "fssi_frozen_model.csv")
MODEL["gene"] = MODEL.gene.str.upper()
GENES = MODEL.gene.tolist()
WEIGHTS = MODEL.set_index("gene").weight

SOURCE_CONTEXT = {
    "PCS-201-021": "normal",
    "CRL-1762": "keloid",
    "PK1": "keloid",
    "PK2": "keloid",
    "PK3": "keloid",
    "PK4": "keloid",
    "PNF4": "normal",
}


def fixed_score(expression: pd.DataFrame) -> pd.Series:
    available = [gene for gene in GENES if gene in expression.index]
    matrix = expression.loc[available].astype(float)
    gene_sd = matrix.std(axis=1, ddof=1).replace(0, np.nan)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(gene_sd, axis=0).fillna(0)
    weights = WEIGHTS.loc[available]
    return z.mul(weights, axis=0).sum(axis=0) / weights.abs().sum()


def exact_sign_flip(values: np.ndarray) -> float:
    observed = abs(float(np.mean(values)))
    null = [abs(float(np.mean(values * signs))) for signs in itertools.product([-1, 1], repeat=len(values))]
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def t_interval(values: np.ndarray) -> tuple[float, float]:
    if len(values) < 2:
        return math.nan, math.nan
    margin = float(stats.t.ppf(0.975, len(values) - 1) * stats.sem(values))
    return float(np.mean(values) - margin), float(np.mean(values) + margin)


def summarise(context: str, values: pd.Series) -> dict[str, object]:
    vector = values.to_numpy(float)
    low, high = t_interval(vector)
    return {
        "dataset": "GSE232079",
        "context": context,
        "contrast": "DMSO_minus_corin",
        "n_cell_sources": len(vector),
        "mean_response": float(np.mean(vector)),
        "ci95_low": low,
        "ci95_high": high,
        "positive_sources": int(np.sum(vector > 0)),
        "exact_signflip_p_two_sided": exact_sign_flip(vector),
    }


def sample_metadata(columns: list[str]) -> pd.DataFrame:
    rows = []
    for column in columns:
        source = next((candidate for candidate in SOURCE_CONTEXT if column.startswith(candidate)), None)
        if source is None:
            raise ValueError(f"Unrecognised GSE232079 sample name: {column}")
        treatment = "corin" if re.search(r"\bcorin\b", column, flags=re.IGNORECASE) else "DMSO"
        replicate = int(re.search(r"(\d+)$", column).group(1))
        passage_match = re.search(r"\bP(\d+)\b", column)
        rows.append(
            {
                "sample_id": column,
                "cell_source": source,
                "context": SOURCE_CONTEXT[source],
                "treatment": treatment,
                "replicate": replicate,
                "passage": int(passage_match.group(1)) if passage_match else math.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    path = STAGE / "GSE232079_varianceStabilizingTransformation.txt.gz"
    frame = pd.read_csv(path, sep="\t")
    frame = frame.rename(columns={frame.columns[0]: "ENSEMBL"})
    frame["ENSEMBL"] = frame.ENSEMBL.astype(str).str.split(".").str[0]
    mapping = pd.read_csv(STAGE.parent / "external" / "GSE191067" / "all_ensembl_symbol_map.csv")
    mapping["ENSEMBL"] = mapping.ENSEMBL.astype(str).str.split(".").str[0]
    merged = frame.merge(mapping[["ENSEMBL", "SYMBOL"]].drop_duplicates(), on="ENSEMBL", how="left")
    merged["SYMBOL"] = merged.SYMBOL.astype(str).str.upper()
    sample_columns = frame.columns[1:].tolist()
    expression = merged.groupby("SYMBOL")[sample_columns].mean()
    score = fixed_score(expression)

    metadata = sample_metadata(sample_columns)
    metadata["frozen_score"] = metadata.sample_id.map(score)
    metadata["genes_detected"] = len(expression.index.intersection(GENES))
    metadata["genes_total"] = len(GENES)
    metadata.to_csv(TABLES / "GSE232079_corin_sample_scores.csv", index=False)

    source_scores = metadata.groupby(["cell_source", "context", "treatment"], as_index=False).frozen_score.mean()
    wide = source_scores.pivot(index=["cell_source", "context"], columns="treatment", values="frozen_score").reset_index()
    wide["DMSO_minus_corin"] = wide.DMSO - wide.corin
    wide.to_csv(TABLES / "GSE232079_corin_source_responses.csv", index=False)
    summaries = [summarise(context, subset.DMSO_minus_corin) for context, subset in wide.groupby("context")]

    summary = {
        "dataset": "GSE232079",
        "source_file": path.name,
        "source_file_sha256": "065cd789d80f613f5aa3c709693150829e903f0126f9de05f89dc125c9d8789c",
        "processed_matrix_dimensions": [int(frame.shape[0]), int(len(sample_columns))],
        "biological_unit": "cell source; two expression replicates were averaged within treatment",
        "cell_source_context": SOURCE_CONTEXT,
        "genes_detected": int(len(expression.index.intersection(GENES))),
        "genes_total": len(GENES),
        "summaries": summaries,
    }
    (TABLES / "GSE232079_corin_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True), encoding="utf-8"
    )
    print(wide.to_string(index=False))
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
