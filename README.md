# AML data jamboree: data processing

Harmonized multi-omics leukemia resource, AML first. Three R scripts turn the public datasets in
`data/` into processed, harmonized data in `results/`. The web portal
([AML_data_visualization](../AML_data_visualization), Python/Streamlit) reads only `results/`.

```
prepare_cell_line_data.R           data/cell_lines/  → results/cell_lines/              [working]
prepare_patient_data.R             data/patients/    → results/patients/<DATASET_ID>/   [working]
prepare_single_cell_RNAseq_data.R  data/single_cell/ → results/single_cell/<DATASET_ID>/ [planned]
data/                              inputs (curated files in git, raw downloads gitignored; data/README.md)
results/                           processed/harmonized data for the portal (docs/results_contract.md)
catalog/datasets.tsv               inventory of every dataset (planned → processed)
schema/                            harmonized study/sample annotation (OncoTree, NCIt, HGNC)
```

## Run

From the repo root:

```bash
Rscript prepare_cell_line_data.R       # ~ 5 min; needs the raw downloads in data/cell_lines/
Rscript prepare_patient_data.R         # ~ 1 min; needs data/patients/ (and DepMap PortalCompounds.csv for drug targets)
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

## Patient data (`prepare_patient_data.R`)

Six cohorts, written to `results/patients/<DATASET_ID>/` in the same layout as `results/cell_lines/`
(one column per sample in `values/`):

| DATASET_ID | cohort | samples / patients | layers | survival |
|---|---|---|---|---|
| `BEAT_AML` | Beat AML waves 1–4 (Bottomly et al. 2022) | 942 / 805 | RNA (671), drug AUC, 166 inhibitors (631) | OS |
| `LEE_2018` | Lee et al. 2018 Nat Commun | 30 / 30 | RNA array (29), drug AUC, 159 drugs (30) | none (no event indicator) |
| `WANG_2022_PED` | Wang et al. 2022 Blood Cancer Discov, pediatric | 66 / 57 | RNA (48), drug AUC, 45 drugs (61) | EFS |
| `MALANI_2022` | Malani et al. 2022 Cancer Discov (FPM tumor board) | 252 / 186 | RNA (163), DSS and sDSS, 515 drugs (164) | OS, EFS |
| `PTRC_BEAT_AML` | Beat AML proteomics, OHSU / PNNL PTRC (Pino et al. 2024) | 210 / 210 | protein (210), phosphosites (210), drug AUC from Beat AML (167) | OS (from Beat AML) |
| `CPTAC_AML_ADULT` | CPTAC adult AML, data freeze v8 (Chu et al. 2026) | 162 / 162 | protein (162), phosphosites (162) | none |

- **Annotation** (`schema/sample_annotation.tsv`): population, age, sex, OncoTree code, tissue,
  timepoint (diagnosis / relapse / refractory / remission / residual), risk group with its classification
  (ELN 2017, ELN 2010 or cytogenetic), KMT2A rearrangement, induction response, OS / EFS in days.
- **Mutations:** `mut_<GENE>` for the same 20 AML genes as the cell lines; NA when the gene was not
  assayed for that sample (Lee 2018 has NPM1 and FLT3 only). Beat AML: WES calls plus clinical FLT3-ITD
  and NPM1; Wang 2022: panel variants at VAF ≥ 5% plus clinical FLT3 / NPM1 / CEBPA; Malani 2022: binary
  mutation table plus clinical FLT3-ITD / TKD.
- **MYC / mTOR RNA signatures:** mean z-score of the genes in `data/*_signature_RNA.tsv` within each
  cohort; `MYC_high` / `mTOR_high` = above the cohort median. This reproduces the published cell-line
  classification closely but not exactly (r = 0.90 MYC, 0.82 mTOR on the 56 DepMap lines).
- **Drugs:** `features/drug_*.tsv` adds `compound_key` (normalized name, matches the same drug across
  cohorts), HGNC `targets` from DepMap `PortalCompounds.csv` (name or synonym match), `mechanism`,
  `drug_class`. `lower_is_stronger` is TRUE for AUC and FALSE for DSS / sDSS.
- **Genes:** HGNC symbols, protein-coding only; Ensembl IDs mapped with the Beat AML annotation;
  duplicated probes / IDs averaged.
- **Proteomics:** `protein` (gene-level log2 ratio) and `phospho` (single-site log2 values; `features/phospho.tsv`
  adds `site` and `sequence_window`). CPTAC: data freeze v8 `prot_ratio_gene_*_RBC_corrected` and
  `phospho_ratio_single-site_*_RBC_corrected` (not imputed). Acetyl, glyco and succinyl sites are not used yet.
- **PTRC and Beat AML:** PTRC `sample_id` is the Beat AML labId. The 167 specimens in Beat AML waves 1–4 take
  annotation, survival and drug AUC from the matching Beat AML sample (`beat_aml_sample_id`); RNA is not
  repeated. The PTRC drug layer has `duplicate_of = BEAT_AML/drug_auc` in `datasets.tsv` so analyses can
  avoid counting it twice.
