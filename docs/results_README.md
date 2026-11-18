# Harmonized AML multi-omic data

This folder holds the harmonized data behind the **AML Multi-Omic Data Analysis Portal**. It was built for
the 2026 NCI Data Jamboree project *Building a Harmonized Multi-Omics Data Resource and Web Application for
Biomarker Discovery and Drug Target Nomination in Leukemia*.

Public acute myeloid leukemia (AML) datasets are brought into one layout and one vocabulary, so that patient
cohorts and cell lines can be analysed side by side:

- **70 AML cell lines** with 13 data types: RNA, CRISPR and RNAi dependency, three drug screens, methylation,
  miRNA, chromatin marks, metabolomics, copy number and the cell-surface proteome.
- **15 patient cohorts**, adult and pediatric: 4,064 samples from 3,600 patients, with RNA, protein, surface protein,
  phosphosites, mutations, clinical data and ex vivo drug response.

The data was produced by `prepare_cell_line_data.R` and `prepare_patient_data.R` in
[AML_data_processing](https://github.com/2026-NCI-Data-Jamboree/AML_data_processing). It is read by the
portal in [AML_data_visualization](https://github.com/2026-NCI-Data-Jamboree/AML_data_visualization).
Do not edit files here by hand: rerun the scripts instead. When this README sits in `results/`, it is a copy
of `docs/results_README.md`.

## Folder layout

```
results/
  README.md                 this file
  dataset_sources.tsv       every integrated dataset: reference, source, accession, NCI resource
  catalog/datasets.tsv      dataset inventory (see catalog/README.md)
  schema/                   harmonized annotation fields and allowed values (see schema/README.md)
  cell_lines/               all AML cell lines, one folder
  patients/<DATASET_ID>/    one folder per patient cohort
  single_cell/<REFERENCE_ID>/  the single-cell reference used for deconvolution (see Deconvolution)
```

Each data folder (`cell_lines/`, and each `patients/<DATASET_ID>/`) has the same files:

| file | what it contains |
|---|---|
| `samples.tsv` | One row per sample, with harmonized annotation: age, sex, disease, timepoint, risk group, mutations, KMT2A, survival, and more. The first column (`cell_line` or `sample_id`) is the sample ID used as the column names in `values/`. |
| `datasets.tsv` | One row per data layer, with its `key`, label, modality, unit, source, the number of features and samples, and `lower_is_stronger`. |
| `features/<key>.tsv` | One row per feature (gene, phosphosite, drug, ...), with a display label and HGNC `gene_symbol`. Drug layers add `targets`, `mechanism` and `compound_key`. |
| `values/<key>.tsv.gz` | The data matrix: rows are features (first column `feature`), columns are samples. |
| `groupings.tsv` | The annotation columns that can be used to compare groups of samples, with display labels. |
| `manifest.tsv` | Every file in the folder, with a description. |
| `provenance.tsv` | When and how the folder was made: date, code commit, R version, data version and notes. |

All tables are tab-separated with a header row. Empty cells mean missing. `.tsv.gz` files are gzip-compressed.

## Patient cohorts

| `DATASET_ID` | cohort | population | samples / patients | data layers (samples) | survival |
|---|---|---|---|---|---|
| `BEAT_AML` | Beat AML waves 1–4 (Bottomly et al. 2022) | adult | 942 / 805 | RNA (671); drug AUC, 166 inhibitors (631) | OS |
| `PTRC_BEAT_AML` | Beat AML proteomics, PTRC (Pino et al. 2024) | adult | 210 / 210 | protein (210); phosphosites (210); drug AUC from Beat AML (167) | OS |
| `CPTAC_AML_ADULT` | CPTAC adult AML, data freeze v8 (Chu et al. 2026) | adult | 162 / 162 | protein (162); phosphosites (162) | none |
| `ZENODO_14160480` | AML proteogenomics, Zenodo record 14160480 | adult | 374 / 374 | RNA (361); protein (374); phosphosites (217); drug IC50, 41 drugs (102) | OS, EFS |
| `LEE_2018` | Lee et al. 2018 Nat Commun | adult | 30 / 30 | RNA array (29); drug AUC, 159 drugs (30) | none |
| `MALANI_2022` | FPM tumor board (Malani et al. 2022) | adult | 252 / 186 | RNA (163); DSS and sDSS, 515 drugs (164) | OS, EFS |
| `WANG_2022_PED` | Pediatric AML pharmacogenomics (Wang et al. 2022) | pediatric | 66 / 57 | RNA (48); drug AUC, 45 drugs (61) | EFS |
| `AML_STJUDE_2024` | St. Jude pediatric AML (Umeda et al. 2024) | pediatric | 887 / 887 | RNA (887) | none |
| `TCGA_LAML` | TCGA LAML (from the AMLmap compendium) | adult | 150 / 150 | RNA (150) | OS |
| `TARGET_AML` | TARGET AML (from the AMLmap compendium) | pediatric | 145 / 145 | RNA (145) | OS |
| `LEUCEGENE` | Leucegene (from the AMLmap compendium) | adult | 399 / 399 | RNA (399) | none |
| `LUMC_AML` | LUMC AML (from the AMLmap compendium) | adult | 95 / 95 | RNA (95) | none |
| `JAYAVELU_2022_DISCOVERY` | Jayavelu et al. 2022, discovery cohort | adult | 177 / 177 | protein, label-free (177) | OS |
| `JAYAVELU_2022_VALIDATION` | Jayavelu et al. 2022, validation cohort | adult | 75 / 75 | protein, SILAC, imputed (75) | OS |
| `LEUCEGENE_SURFACEOME_2024` | Leucegene AML surfaceome (Bordeleau et al. 2024 Cell Rep) | adult | 100 / 100 | surface protein (100) | none |

Each cohort's `provenance.tsv` has notes on how it was processed. References, sources and accessions are in
`dataset_sources.tsv`.

## Deconvolution

Cell type and cell state proportions of the bulk RNA samples in some cohorts, estimated against the AML scAtlas
single-cell reference (bone marrow, diagnosis samples). They are added by `prepare_deconvolution.py` as ordinary
layers of the cohort (modality `deconvolution`, unit `proportion`), so they join with every other layer by `sample_id`.

- **Layers:** `values/deconv_<reference>_<method>_<level>.tsv.gz`, e.g. `deconv_aml_scatlas_instaprism_cell_type`.
  Rows are cell types or cell states, columns are samples, proportions sum to 1 per sample. Methods: InstaPrism initial
  fit (cell type and cell state), InstaPrism updated fit (cell type; depends on the other samples in the cohort) and NNLS
  (cell type).
- **`deconvolution_runs.tsv`:** one row per layer: method, level, quantity (`mrna_share`: share of the sample's mRNA, not of
  cells), settings, input source and conversion, run date.
- **`deconvolution_warnings.tsv`:** flags on a run, a sample or a label, with a severity (`info`, `warning`, `caution`)
  and a code: low initial malignant proportion (the updated fit may be unreliable), hemoglobin-dominated samples, low share
  of the reference's genes in the bulk matrix, converted input scale, batch dependence. The `key` column names the layer
  (empty = every layer of the cohort). Warnings flag a result; the values are always kept.
- **`single_cell/<REFERENCE_ID>/`:** `reference.tsv` (cells and reference samples per cell type and state),
  `markers.tsv` (top marker genes per label) and `provenance.tsv` (citation, build, limitations).
- **Cohorts:** `BEAT_AML`, `MALANI_2022`, `WANG_2022_PED`, `AML_STJUDE_2024`, `TCGA_LAML`, `TARGET_AML`, `LEUCEGENE`,
  `LUMC_AML`.
- **Input:** raw counts (Beat AML and the AMLmap cohorts from AMLmap's uncorrected counts; Malani from its raw reads) or
  counts per million converted back from a log2 layer (Wang; St. Jude, whose values were batch-corrected after the CPM
  step and are only approximately CPM). **Transcript-length-normalized data (TPM, RPKM, FPKM) is never used**, so
  cohorts with only such data (Zenodo 14160480) and cell lines are not deconvolved.
- **Read with care:** the scAtlas malignant profile under-represents monocytic and more committed AML, and the proportions are
  estimates, not measurements (in Beat AML the malignant proportion tracks the clinical blast percentage, Spearman 0.4-0.5,
  but is not calibrated).

## Cell lines

`cell_lines/` covers 70 AML cell lines from DepMap 26Q3 and related CCLE resources:

| `key` | data | unit | lines |
|---|---|---|---|
| `expression` | RNA expression | log2(TPM+1) | 56 |
| `crispr` | CRISPR dependency (Chronos) | gene effect; lower = more dependent | 32 |
| `rnai` | RNAi dependency (DEMETER2 v6) | gene effect; lower = more dependent | 23 |
| `drug_PRISM`, `drug_CTRP`, `drug_GDSC2` | drug response | AUC; lower = more sensitive | 22–33 |
| `meth_tss`, `meth_cpg` | promoter and CpG-cluster methylation (RRBS) | beta value | 33 |
| `mirna` | miRNA expression | log2(expression+1) | 35 |
| `chromatin` | global chromatin marks | relative abundance | 34 |
| `metabolomics` | metabolomics | log10 abundance | 37 |
| `cn` | copy number, gene level (DepMap 26Q3, WGS) | log2(relative CN + 1); 1 = no change | 50 |
| `surface_george` | cell-surface proteome, surface biotinylation (George et al. 2026, bioRxiv preprint) | log2 mean abundance; empty = not detected on the surface | 19 |

The `cn` features add `chromosome` and `cytoband`, so genes can be grouped into arms or bands. The
`surface_george` features add `spc` (the paper's surface score, 0–4) and `surfaceome` (TRUE when UniProt places
the protein at the cell surface; the capture also picks up some abundant intracellular proteins).

## How the data is harmonized

- **Genes** use current HGNC symbols in every cohort. Older names are mapped to the approved symbol
  (for example HDGFRP3 → HDGFL3), so a gene can be found under one name everywhere. Patient RNA keeps
  protein-coding genes only. When a gene has several measurements (for example array probes), they are
  averaged.
- **Units are kept as the source reports them** (see `unit` in `datasets.tsv`). Values are not
  rescaled across cohorts, so compare effect sizes within a cohort.
- **Direction of drug and dependency scores:** `lower_is_stronger = TRUE` means a lower value is a
  stronger effect. This applies to AUC, IC50 and gene effect. DSS and sDSS are the other way round:
  higher means more sensitive.
- **Drugs** are matched across cohorts with `compound_key`: the upper-case name without salt forms,
  brackets or punctuation. Targets come from DepMap.
- **Mutations:** `mut_<GENE>` columns cover 20 recurrent AML genes: TP53, FLT3, NPM1, NRAS, KRAS, PTPN11,
  KIT, IDH1, IDH2, DNMT3A, TET2, ASXL1, RUNX1, CEBPA, WT1, SRSF2, U2AF1, SF3B1, STAG2 and BCOR.
  - The columns are TRUE (mutated) or FALSE (not mutated).
  - Empty means the gene was **not tested** for that sample, which is different from not mutated.
- **FLT3-ITD** (cell lines): `FLT3_ITD` is TRUE for an internal tandem duplication, called from the
  variant-level mutation table. `mut_FLT3` also counts ITD and the D835 / I836 TKD mutations.
- **KMT2A rearrangement** (`KMT2A_rearranged`, with `KMT2A_partner` when known) comes from fusion calls,
  karyotype or clinical records, depending on the cohort.
- **Risk group** is favorable, intermediate or adverse, as the cohort reports it. `risk_classification`
  says which system: ELN 2017, ELN 2010 or cytogenetic.
- **Age and population:** `age_years` is age at diagnosis. The portal's age groups follow the NCI bands:
  pediatric < 15, AYA 15–39, adult ≥ 40.
- **Survival:**
  - `overall_survival_days` and `event_free_survival_days` are in days.
  - `os_event` and `efs_event`: TRUE means the event was observed, FALSE means the patient was censored.
- **MYC and mTOR RNA signatures** (cell lines only): `MYC_score` and `mTOR_score`, with `MYC_high` and
  `mTOR_high`, from the published classification. Patient cohorts have no MYC / mTOR groups.
- **Same patient in two cohorts:** PTRC samples are Beat AML specimens. Their clinical data and drug AUC
  come from the matching Beat AML sample, and that drug layer is marked
  `duplicate_of = BEAT_AML/drug_auc` so it is not counted twice.

The allowed values for each harmonized field are in `schema/sample_annotation.tsv` and
`schema/vocabularies.tsv`.

## Reading the data

Python (pandas):

```python
import pandas as pd

d = "results/patients/BEAT_AML"
samples = pd.read_csv(f"{d}/samples.tsv", sep="\t")
layers = pd.read_csv(f"{d}/datasets.tsv", sep="\t")
rna = pd.read_csv(f"{d}/values/expression.tsv.gz", sep="\t", index_col="feature")  # genes x samples
rna.loc["MAP1A", samples["sample_id"]]
```

R:

```r
d <- "results/patients/BEAT_AML"
samples <- read.delim(file.path(d, "samples.tsv"), check.names = FALSE)
auc <- read.delim(gzfile(file.path(d, "values/drug_auc.tsv.gz")), row.names = 1, check.names = FALSE)
drugs <- read.delim(file.path(d, "features/drug_auc.tsv"))   # names, targets, compound_key
```

## Things to keep in mind

- Cohorts differ in how drug response was measured, which drugs were tested and how samples were
  collected.
- Drug response is measured on patient cells in the lab (ex vivo), not in the patient.
- Some cohorts include more than one sample per patient, for example diagnosis and relapse. Use
  `participant_id` and `timepoint` to keep one sample per patient.
- Some pediatric samples may appear in more than one cohort:
  - The St. Jude cohort includes samples from TARGET, which is also a cohort of its own here.
  - The St. Jude cohort and Wang 2022 both come from St. Jude, so they may share patients.
  - The patient IDs cannot be matched between these cohorts.
- Please cite the original publications listed in `dataset_sources.tsv` when you use a cohort.

## Versions

`provenance.tsv` records `results_version`, which is 3 at the time of writing. Adding files or columns
keeps the version. Renaming or removing anything the portal reads raises it. The full file
specification is in `docs/results_contract.md` in the processing repository.
