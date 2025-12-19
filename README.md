# Disease analysis

### StringDB data

To obtain the DISEASES and HPO data from StrongDB, run the following scripts.

```bash
python -m pip install pandas requests openpyxl
```

```bash
python run_stringdb_batch.py
```

```bash
extract_unique_terms.py
```

```bash
concatenate_stringdb_output.py
```

### DISEASES analysis

To reproduce the figures and percentages presented in the paper, run the R script DISEASES_analysis.R.

### HPO analysis

To reproduce the figures and percentages presented in the paper, run the R script HPO_analysis.R.
