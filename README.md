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

# DISEASES_analysis folder

To reproduce the figures and percentages presented in the paper, run the R script DISEASES_analysis.R.

Below is an explanation of two files generated to facilitate an in-depth analysis of the data. 

---

### `diseases_description_counts.csv`

This file summarizes the frequency of disease annotations in the dataset. Each row corresponds to a unique disease term (`description`), together with the number of times it appears in the dataset.

**Contents:**
- `description` – disease name or term  
- `n` – number of occurrences of the disease in the dataset  

---

### `diseases_unique_inputGenes_and_preferredNames.csv`

This file contains a disease-centric aggregation of associated genes. Each row represents a single disease term and all the associated genes.

**Contents:**
- `description` – disease name or term  
- `inputGenes_unique` – unique gene identifiers associated with the disease  
- `preferredNames_unique` – unique preferred gene names associated with the disease

---

# HPO_analysis folder

To reproduce the figures and percentages presented in the paper, run the R script HPO_analysis.R.

Below is an explanation of two files generated to facilitate an in-depth analysis of the data. 

---

### `hpo_description_counts.csv`

This file summarizes the frequency of Human Phenotype Ontology (HPO) terms in the dataset. Each row corresponds to a unique HPO term (`description`), together with the number of times it occurs in the dataset.

**Contents:**
- `description` – HPO term description  
- `n` – number of occurrences of the term in the filtered dataset  

---

### `hpo_unique_inputGenes_and_preferredNames.csv`

This file provides an HPO-centric aggregation of associated genes. Each row represents a single HPO term and all the associated genes.

**Contents:**
- `description` – HPO term description  
- `inputGenes_unique` – unique gene identifiers associated with the HPO term  
- `preferredNames_unique` – unique preferred gene names associated with the HPO term  

---

