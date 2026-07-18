from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / ".stage" / "public_expansion"
OUTPUT = ROOT / "results" / "tables" / "Lou2026_source_data_inventory.json"


def main() -> None:
    inventory = {}
    for path in sorted(STAGE.glob("41467_2026_72823_MOESM*_ESM.xlsx")):
        workbook = pd.ExcelFile(path)
        sheets = []
        for sheet in workbook.sheet_names:
            frame = pd.read_excel(workbook, sheet_name=sheet, header=None, nrows=8)
            sheets.append(
                {
                    "sheet": sheet,
                    "preview_rows": int(frame.shape[0]),
                    "preview_columns": int(frame.shape[1]),
                    "preview": frame.fillna("").astype(str).values.tolist(),
                }
            )
        inventory[path.name] = sheets
    OUTPUT.write_text(json.dumps(inventory, indent=2), encoding="utf-8")
    print(json.dumps(inventory, indent=2))


if __name__ == "__main__":
    main()
