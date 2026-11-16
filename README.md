# AML data jamboree: data processing

Harmonized multi-omics leukemia resource, AML first. Three R scripts turn the public datasets in
`data/` into processed, harmonized data in `results/`. The web portal
([AML_data_visualization](../AML_data_visualization), Python/Streamlit) reads only `results/`.

```
prepare_cell_line_data.R           data/cell_lines/  → results/cell_lines/              [working]
prepare_patient_data.R             data/patients/    → results/patients/<DATASET_ID>/   [planned]
prepare_single_cell_RNAseq_data.R  data/single_cell/ → results/single_cell/<DATASET_ID>/ [planned]
data/                              inputs (curated files in git, raw downloads gitignored; data/README.md)
results/                           processed/harmonized data for the portal (docs/results_contract.md)
catalog/datasets.tsv               inventory of every dataset (planned → processed)
schema/                            harmonized study/sample annotation (OncoTree, NCIt, HGNC)
```

## Run

From the repo root:

```bash
Rscript prepare_cell_line_data.R       # ~ a few minutes; needs the raw downloads in data/cell_lines/
```

Each script also copies `catalog/` and `schema/` into `results/`. To share the results:
`tar -czf aml_results.tar.gz results`.

Then in the portal repo, `streamlit run app.py`. It finds `../AML_data_processing/results`
automatically, or you can set `AML_PORTAL_DATA=/path/to/results`.

## Cell-line data (`prepare_cell_line_data.R`)

- **Cell lines:** all 70 DepMap 26Q3 models with OncotreePrimaryDisease "Acute Myeloid Leukemia".
- **Annotation:** OncoTree code, sex, age/population, mutation status of 20 recurrent AML genes
  (DepMap hotspot + damaging calls), KMT2A rearrangement (curated, else DepMap fusion calls), and
  the MYC/mTOR RNA classification for 56 of the lines (`MYC_high`, `mTOR_high`, `MYC_mTOR_category`).
- **Groupings:** `results/cell_lines/groupings.tsv` lists the annotation columns the portal offers
  for two-group comparisons (volcano plots). To add one, add a column to `samples.tsv` and a row there.
- **Layers:** RNA expression, CRISPR (Chronos) and RNAi (DEMETER2) dependency, drug AUC (PRISM,
  CTRP, GDSC2), RRBS promoter and CpG-cluster methylation, miRNA, global chromatin marks, metabolomics.
- **Filtering:** only features measured in at least 3 lines and not constant are kept.
- **No statistics:** the portal computes group comparisons on the fly.
