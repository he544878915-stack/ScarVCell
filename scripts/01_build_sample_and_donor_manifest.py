"""Lock donor independence and dataset roles before the v2 analyses."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "tables"


def add(rows, dataset, sample, donor, condition, role, independence, note=""):
    rows.append({
        "dataset_id": dataset,
        "sample_id": sample,
        "donor_id": donor,
        "condition": condition,
        "dataset_role": role,
        "independence_status": independence,
        "analysis_unit": "donor",
        "note": note,
    })


def main() -> None:
    rows: list[dict[str, str]] = []

    for i in range(1, 4):
        add(rows, "GSE163973", f"KF{i}_matrix", f"GSE163973_KF{i}", "keloid",
            "single_cell_development", "reported_independent_donor")
    for i in range(1, 4):
        add(rows, "GSE163973", f"NF{i}_matrix", f"GSE163973_NF{i}_UNRESOLVED", "normal_scar",
            "single_cell_development", "donor_mapping_uncertain",
            "GEO exposes three matrices but the protocol text reports two normal-scar donors; strict inference must collapse one pair in sensitivity analyses.")

    for sample, condition in [
        ("GSM5494438_Ke01", "keloid"),
        ("GSM5494439_Ke02", "keloid"),
        ("GSM5610149_NS02", "normal_scar"),
    ]:
        add(rows, "GSE181297", sample, f"GSE181297_{sample.split('_')[-1]}", condition,
            "single_cell_development", "reported_independent_donor")

    for sample, condition in [("GSM7794710_K", "keloid"), ("GSM7794711_H", "hypertrophic_scar")]:
        add(rows, "GSE243716", sample, "GSE243716_P1", condition,
            "single_cell_development", "paired_same_donor",
            "Keloid and hypertrophic scar tissues were collected from the same 34-year-old participant.")

    for cohort in ("GR1", "PR1", "GR2", "PR2"):
        for timepoint in ("D0", "W3", "M6", "Y1"):
            add(rows, "GSE320017", f"{cohort}{timepoint}", f"GSE320017_{cohort}",
                "burn_hypertrophic_scar_laser_series", "retrospective_context_validation",
                "repeated_longitudinal_same_donor",
                "Four longitudinal matrices belong to one participant and must not be treated as independent samples.")

    for dataset in ("GSE92566", "GSE44270", "GSE7890"):
        add(rows, dataset, "TO_BE_PARSED", "TO_BE_PARSED", "keloid_or_control",
            "locked_external_validation", "pending_metadata_audit",
            "Locked before expression retrieval; no feature or weight tuning is permitted after projection.")

    donor = pd.DataFrame(rows)
    donor.to_csv(OUT / "sample_and_donor_manifest.csv", index=False)

    design = pd.DataFrame([
        {"dataset_id": "GSE163973", "role": "development", "allowed_use": "feature_selection;weight_estimation;internal_effect"},
        {"dataset_id": "GSE181297", "role": "development", "allowed_use": "state_mapping;weight_estimation;communication"},
        {"dataset_id": "GSE243716", "role": "development", "allowed_use": "paired_state_mapping;weight_estimation"},
        {"dataset_id": "GSE320017", "role": "retrospective_context", "allowed_use": "four-donor longitudinal context only"},
        {"dataset_id": "GSE307504", "role": "retrospective_context", "allowed_use": "coordinate-aware spatial falsification"},
        {"dataset_id": "GSE92566", "role": "locked_external_validation", "allowed_use": "frozen-score projection only"},
        {"dataset_id": "GSE44270", "role": "locked_external_validation", "allowed_use": "frozen-score projection;culture-stratified"},
        {"dataset_id": "GSE7890", "role": "locked_external_validation", "allowed_use": "frozen-score projection only"},
    ])
    design["lock_date"] = "2026-07-11"
    design["prohibited_use"] = "post-projection feature selection or weight tuning"
    design.to_csv(OUT / "validation_design_manifest.csv", index=False)

    print(f"Wrote {len(donor)} donor/sample rows and {len(design)} dataset-role rows")


if __name__ == "__main__":
    main()
