"""Project frozen FSSI into the measured siRUNX2 keloid-fibroblast experiment."""

from __future__ import annotations

import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = Path(os.environ.get("SCARVCELL_STAGE", str(ROOT / ".stage"))) / "GSE293677" / "GSE293677_TPM.tab.gz"
T = ROOT / "results" / "tables"


def exact_sign_flip(d: np.ndarray) -> tuple[float, float]:
    observed = float(np.mean(d))
    null = np.asarray([np.mean(d * np.asarray(s)) for s in itertools.product([-1, 1], repeat=len(d))])
    return float(np.mean(null >= observed - 1e-12)), float(np.mean(np.abs(null) >= abs(observed) - 1e-12))


def main() -> None:
    model = pd.read_csv(T / "fssi_frozen_model.csv")
    raw = pd.read_csv(INPUT, sep="\t", compression="gzip")
    raw["gene"] = raw.FEATURE_NAME.astype(str).str.split("|").str[-1].str.upper()
    expression = raw.groupby("gene")[raw.columns[1:7]].median()
    genes = [g for g in model.gene if g in expression.index]
    x = np.log1p(expression.loc[genes].T)
    z = (x - x.mean()) / x.std(ddof=1).replace(0, 1.0).fillna(1.0)
    m = model.set_index("gene").loc[genes]
    w = m.weight.to_numpy(float).copy(); w /= np.abs(w).sum()
    ew = m.equal_weight.to_numpy(float).copy(); ew /= np.abs(ew).sum()
    scores = pd.DataFrame({
        "sample_id": x.index,
        "pair_id": [1, 2, 3, 1, 2, 3],
        "treatment": ["siControl"] * 3 + ["siRUNX2"] * 3,
        "fssi_weighted": z.to_numpy() @ w,
        "fssi_equal_weight": z.to_numpy() @ ew,
        "genes_present": len(genes),
    })
    rows = []
    for score_col in ("fssi_weighted", "fssi_equal_weight"):
        wide = scores.pivot(index="pair_id", columns="treatment", values=score_col)
        reversal = (wide.siControl - wide.siRUNX2).to_numpy(float)
        p1, p2 = exact_sign_flip(reversal)
        rows.append({
            "dataset_id": "GSE293677", "target": "RUNX2", "perturbation": "siRUNX2",
            "score_variant": score_col, "paired_donor_n": len(reversal),
            "mean_sri_reduction": float(reversal.mean()), "median_sri_reduction": float(np.median(reversal)),
            "support_fraction": float(np.mean(reversal > 0)),
            "exact_p_one_sided_expected": p1, "exact_p_two_sided": p2,
            "inference": "measured_target_specific_expression_response_small_n",
        })

    contributions = z.mul(w, axis=1)
    contribution_rows = []
    for pair in (1, 2, 3):
        ctrl, kd = f"CTRL{pair}.TPM", f"KD{pair}.TPM"
        for gene in genes:
            contribution_rows.append({"pair_id": pair, "gene": gene,
                                      "weighted_contribution_change_kd_minus_control": float(contributions.loc[kd, gene] - contributions.loc[ctrl, gene])})

    runx2 = np.log1p(expression.loc["RUNX2"]) if "RUNX2" in expression.index else pd.Series(dtype=float)
    target_qc = pd.DataFrame({
        "pair_id": [1, 2, 3],
        "runx2_log1p_tpm_control": [runx2.get(f"CTRL{i}.TPM", np.nan) for i in (1, 2, 3)],
        "runx2_log1p_tpm_knockdown": [runx2.get(f"KD{i}.TPM", np.nan) for i in (1, 2, 3)],
    })
    target_qc["runx2_change_kd_minus_control"] = target_qc.runx2_log1p_tpm_knockdown - target_qc.runx2_log1p_tpm_control

    scores.to_csv(T / "runx2_measured_perturbation_sri_scores.csv", index=False)
    pd.DataFrame(rows).to_csv(T / "runx2_measured_perturbation_summary.csv", index=False)
    pd.DataFrame(contribution_rows).to_csv(T / "runx2_measured_perturbation_gene_contributions.csv", index=False)
    target_qc.to_csv(T / "runx2_knockdown_expression_qc.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print(target_qc.to_string(index=False))


if __name__ == "__main__":
    main()
