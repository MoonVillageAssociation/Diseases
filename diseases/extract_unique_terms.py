import pandas as pd
from pathlib import Path

IN_DIR = Path("StringDB_results")
OUT_HPO = Path("HPO_analysis/HPO_terms.csv")
OUT_DIS = Path("DISEASES_analysis/DISEASES_terms.csv")

# Ensure output directories exist
OUT_HPO.parent.mkdir(parents=True, exist_ok=True)
OUT_DIS.parent.mkdir(parents=True, exist_ok=True)

def collect_unique(pattern: str) -> pd.DataFrame:
    files = sorted(IN_DIR.glob(pattern))
    rows = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t", dtype=str)
        except Exception:
            continue
        if not {"term", "description"} <= set(df.columns):
            continue
        df["term"] = df["term"].astype(str).str.strip()
        df["description"] = df["description"].astype(str).str.strip()
        rows.append(df[["term", "description"]])
    if not rows:
        return pd.DataFrame(columns=["term", "description"])
    out = pd.concat(rows, ignore_index=True)
    out = out.drop_duplicates(subset=["term", "description"]).sort_values(["term", "description"])
    out = out.reset_index(drop=True)
    return out

# Collect and save
hpo_df = collect_unique("*__HPO.tsv")
dis_df = collect_unique("*__DISEASES.tsv")

hpo_df.to_csv(OUT_HPO, index=False)
dis_df.to_csv(OUT_DIS, index=False)

print(f"Saved {len(hpo_df)} unique HPO rows -> {OUT_HPO}")
print(f"Saved {len(dis_df)} unique DISEASES rows -> {OUT_DIS}")
