import pandas as pd
from pathlib import Path

IN_DIR = Path("StringDB_results")
IN_DIR.mkdir(parents=True, exist_ok=True)  # just in case

def combine_stringdb(pattern: str, out_csv: str, label: str):
    files = sorted(IN_DIR.glob(pattern))
    if not files:
        print(f"[{label}] No files matched pattern: {pattern}")
        return

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)

            # Extract ChemicalName from filename before the last "__<LABEL>"
            # e.g., "Zinc__DISEASES.tsv" -> "Zinc"
            stem = f.stem  # e.g. "Zinc__DISEASES"
            chemical_name = stem.rsplit("__", 1)[0]

            # Clean up expected cols if present
            if "term" in df.columns:
                df["term"] = df["term"].astype(str).str.strip()
            if "description" in df.columns:
                df["description"] = df["description"].astype(str).str.strip()

            # Add ChemicalName as the first column
            df.insert(0, "ChemicalName", chemical_name)

            dfs.append(df)
        except Exception as e:
            print(f"[{label}] Skipping {f.name} due to error: {e}")

    if not dfs:
        print(f"[{label}] No valid dataframes loaded; nothing to save.")
        return

    combined = pd.concat(dfs, ignore_index=True)

    # Drop exact duplicate rows (optional but usually helpful)
    before = len(combined)
    combined = combined.drop_duplicates()
    after = len(combined)

    out_path = IN_DIR / out_csv
    combined.to_csv(out_path, index=False)

    print(
        f"[{label}] Combined {len(files)} files → {out_path}\n"
        f"         Rows (raw): {before:,} | Rows (unique): {after:,}"
    )

# Build ALL_DISEASES.csv
combine_stringdb(pattern="*__DISEASES.tsv", out_csv="ALL_DISEASES.csv", label="DISEASES")

# Build ALL_HPO.csv
combine_stringdb(pattern="*__HPO.tsv", out_csv="ALL_HPO.csv", label="HPO")
