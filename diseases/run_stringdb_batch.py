import os
import time
import re
import sys
import json
from pathlib import Path

import pandas as pd
import requests

# ---------------- CONFIG ----------------
INPUT_XLSX = "chem_gen_dis_astdegprot_pval_log2fc_filt.xlsx"
OUT_DIR = "StringDB_results"
SPECIES = 9606  # Human
CALLER_ID = "anna_batch_stringdb_script"
API_DELAY = 1.0
MIN_MAPPED_FOR_ENRICH = 3  # STRING needs ≥3 IDs for enrichment
# ---------------------------------------

STRING_BASE = "https://string-db.org/api"

def sanitize_filename(name: str) -> str:
    name = re.sub(r"[^\w\-]+", "_", str(name).strip())
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "unnamed"

def string_post(method: str, params: dict, out_format: str = "json") -> requests.Response:
    url = f"{STRING_BASE}/{out_format}/{method}"
    p = dict(params or {})
    p.setdefault("caller_identity", CALLER_ID)
    r = requests.post(url, data=p, timeout=120)
    r.raise_for_status()
    return r

def map_gene_symbols_to_string_ids(genes, species: int):
    """
    Map a list of HGNC symbols to STRING IDs using newline separators (CR).
    We also set limit=1 so each input yields at most 1 best mapping.
    """
    if not genes:
        return [], {}

    method = "get_string_ids"
    # IMPORTANT: use real newlines for POST, not "%0d"
    identifiers = "\r".join(genes)

    params = {
        "identifiers": identifiers,
        "species": species,
        "limit": 1,         # best hit per query
        "echo_query": 1     # returns the original queryItem so we can track
    }

    resp = string_post(method, params, out_format="tsv")
    lines = resp.text.strip().splitlines()
    if not lines:
        return [], {}

    header = lines[0].split("\t")
    # Robustly locate columns
    cols = {name: i for i, name in enumerate(header)}
    def col(name, default=None):
        return cols.get(name, default)

    idx_stringId = col("stringId", 2)
    idx_pref     = col("preferredName", 5)
    idx_query    = col("queryItem", 0)

    string_ids = []
    mapping = {}

    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) <= max(idx_stringId, idx_pref, idx_query):
            continue
        string_id  = parts[idx_stringId].strip()
        pref       = parts[idx_pref].strip() if idx_pref is not None else ""
        query_item = parts[idx_query].strip() if idx_query is not None else ""
        if string_id:
            string_ids.append(string_id)
            if query_item:
                mapping[query_item] = pref or query_item

    # Deduplicate while preserving order
    seen = set()
    uniq_ids = []
    for sid in string_ids:
        if sid not in seen:
            seen.add(sid)
            uniq_ids.append(sid)

    return uniq_ids, mapping

def run_enrichment(string_ids, species: int):
    if not string_ids:
        return []
    method = "enrichment"
    identifiers = "\r".join(string_ids)  # use real newline separator
    params = {"identifiers": identifiers, "species": species}
    resp = string_post(method, params, out_format="json")
    try:
        data = resp.json()
    except json.JSONDecodeError:
        return []
    if isinstance(data, dict) and data.get("error"):
        return []
    return data if isinstance(data, list) else []

def filter_categories(rows):
    rows_hpo, rows_dis = [], []
    for r in rows:
        cat = str(r.get("category", "")).strip().upper()
        if cat == "HPO":
            rows_hpo.append(r)
        elif cat in {"DISEASES", "DISEASE"}:
            rows_dis.append(r)
    return rows_hpo, rows_dis

def normalize_rows(rows):
    """Convert list fields to comma-joined strings for TSV output."""
    if not rows:
        return rows
    norm = []
    for r in rows:
        rr = dict(r)
        if isinstance(rr.get("inputGenes"), list):
            rr["inputGenes"] = ",".join(rr["inputGenes"])
        if isinstance(rr.get("preferredNames"), list):
            rr["preferredNames"] = ",".join(rr["preferredNames"])
        norm.append(rr)
    return norm

def main():
    in_path = Path(INPUT_XLSX)
    if not in_path.exists():
        sys.exit(f"[ERROR] Input file not found: {in_path}")

    out_dir = Path(OUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(in_path)
    need = {"ChemicalName", "GeneSymbol"}
    if not need <= set(df.columns):
        sys.exit(f"[ERROR] Missing required columns: {need - set(df.columns)}")

    # Prepare groups
    df["GeneSymbol"] = df["GeneSymbol"].astype(str).str.strip()
    df = df[df["GeneSymbol"].str.len() > 0]

    def split_symbols(s):
        # Split on commas/semicolons/pipes/whitespace, but keep hyphens (e.g., HLA-DRA)
        return [p for p in re.split(r"[,\;\|\s]+", s.strip()) if p]

    grouped = {}
    for chem, sub in df.groupby("ChemicalName", dropna=False):
        genes = []
        for v in sub["GeneSymbol"]:
            genes.extend(split_symbols(v))
        # unique, preserve order
        seen = set()
        uniq = []
        for g in genes:
            if g not in seen:
                seen.add(g)
                uniq.append(g)
        grouped[chem if pd.notna(chem) else "NA"] = uniq

    print(f"[INFO] Found {len(grouped)} ChemicalName groups.")
    summary = []

    for i, (chem, genes) in enumerate(grouped.items(), start=1):
        safe = sanitize_filename(chem)
        print(f"[{i}/{len(grouped)}] {chem} -> {len(genes)} genes")

        # Map
        try:
            string_ids, mapping = map_gene_symbols_to_string_ids(genes, SPECIES)
        except Exception as e:
            print(f"    [ERROR] ID mapping failed: {e}")
            time.sleep(API_DELAY)
            continue

        print(f"    Mapped {len(string_ids)}/{len(genes)} genes to STRING IDs")

        if len(string_ids) < MIN_MAPPED_FOR_ENRICH:
            print(f"    Skipping enrichment (<{MIN_MAPPED_FOR_ENRICH} mapped IDs)")
            summary.append({
                "ChemicalName": chem,
                "n_input_genes": len(genes),
                "n_mapped_string_ids": len(string_ids),
                "has_HPO": False,
                "has_DISEASES": False
            })
            continue

        time.sleep(API_DELAY)

        # Enrichment
        try:
            enr = run_enrichment(string_ids, SPECIES)
        except Exception as e:
            print(f"    [ERROR] Enrichment failed: {e}")
            time.sleep(API_DELAY)
            continue

        rows_hpo, rows_dis = filter_categories(enr)
        rows_hpo = normalize_rows(rows_hpo)
        rows_dis = normalize_rows(rows_dis)

        if rows_hpo:
            pd.DataFrame(rows_hpo).to_csv(out_dir / f"{safe}__HPO.tsv", sep="\t", index=False)
            print(f"    Saved HPO ({len(rows_hpo)} rows)")
        else:
            print("    No HPO terms")

        if rows_dis:
            pd.DataFrame(rows_dis).to_csv(out_dir / f"{safe}__DISEASES.tsv", sep="\t", index=False)
            print(f"    Saved DISEASES ({len(rows_dis)} rows)")
        else:
            print("    No DISEASES terms")

        summary.append({
            "ChemicalName": chem,
            "n_input_genes": len(genes),
            "n_mapped_string_ids": len(string_ids),
            "has_HPO": bool(rows_hpo),
            "has_DISEASES": bool(rows_dis)
        })

        time.sleep(API_DELAY)

    if summary:
        pd.DataFrame(summary).to_csv(out_dir / "summary.tsv", sep="\t", index=False)
        print("[INFO] Wrote summary.tsv")

if __name__ == "__main__":
    main()
