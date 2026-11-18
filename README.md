# AML data jamboree: data processing

Harmonized multi-omics leukemia resource, AML first. Three scripts turn the public datasets in
`data/` into processed, harmonized data in `results/`. The web portal
([AML_data_visualization](../AML_data_visualization), Python/Streamlit) reads only `results/`.

```
prepare_cell_line_data.R           data/cell_lines/  → results/cell_lines/              [working]
prepare_patient_data.R             data/patients/    → results/patients/<DATASET_ID>/   [working]
prepare_deconvolution.py           results/ + data/  → deconvolution layers, results/single_cell/<REFERENCE_ID>/ [working; run last]
data/                              inputs (curated files in git, raw downloads gitignored; data/README.md)
results/                           processed/harmonized data for the portal (results/README.md; spec in docs/results_contract.md)
catalog/datasets.tsv               inventory of the integrated datasets
schema/                            harmonized study/sample annotation (OncoTree, NCIt, HGNC)
```

## Run

Environments: the deconvolution tool in `deconvolution/` has its own renv (R 4.6.2, Bioconductor, InstaPrism): start R
there and run `renv::restore()` once. `prepare_deconvolution.py` is a Python script with inline dependencies (pandas),
run by [uv](https://docs.astral.sh/uv/). The two `prepare_*.R` scripts need `xml2` and have no renv (a repo-root renv
can be added later).

From the repo root:

```bash
Rscript prepare_cell_line_data.R       # ~ 5 min; needs the raw downloads in data/cell_lines/
Rscript prepare_patient_data.R         # ~ 5 min; needs perl; needs data/patients/ (and DepMap PortalCompounds.csv for drug targets)
uv run prepare_deconvolution.py        # run last; ~ 10 min for all cohorts; needs data/reference/ and the tool's renv (deconvolution/)
```

Each script also copies `catalog/`, `schema/`, `dataset_sources.tsv` and the data README (`docs/results_README.md` → `results/README.md`) into `results/`. To share the results:
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

Fifteen cohorts, written to `results/patients/<DATASET_ID>/` in the same layout as `results/cell_lines/`
(one column per sample in `values/`):

| DATASET_ID | cohort | samples / patients | layers | survival |
|---|---|---|---|---|
| `BEAT_AML` | Beat AML waves 1–4 (Bottomly et al. 2022) | 942 / 805 | RNA (671), drug AUC, 166 inhibitors (631) | OS |
| `LEE_2018` | Lee et al. 2018 Nat Commun | 30 / 30 | RNA array (29), drug AUC, 159 drugs (30) | none (no event indicator) |
| `WANG_2022_PED` | Wang et al. 2022 Blood Cancer Discov, pediatric | 66 / 57 | RNA (48), drug AUC, 45 drugs (61) | EFS |
| `MALANI_2022` | Malani et al. 2022 Cancer Discov (FPM tumor board) | 252 / 186 | RNA (163), DSS and sDSS, 515 drugs (164) | OS, EFS |
| `PTRC_BEAT_AML` | Beat AML proteomics, OHSU / PNNL PTRC (Pino et al. 2024) | 210 / 210 | protein (210), phosphosites (210), drug AUC from Beat AML (167) | OS (from Beat AML) |
| `CPTAC_AML_ADULT` | CPTAC adult AML, data freeze v8 (Chu et al. 2026) | 162 / 162 | protein (162), phosphosites (162) | none |
| `AML_STJUDE_2024` | St. Jude pediatric AML (Umeda et al. 2024 Nat Genet), cBioPortal | 887 / 887 | RNA (887) | none |
| `TCGA_LAML` | TCGA LAML, from the AMLmap compendium | 150 / 150 | RNA (150) | OS |
| `TARGET_AML` | TARGET AML (pediatric), from the AMLmap compendium | 145 / 145 | RNA (145) | OS |
| `LEUCEGENE` | Leucegene, from the AMLmap compendium | 399 / 399 | RNA (399) | none |
| `LUMC_AML` | LUMC AML, from the AMLmap compendium | 95 / 95 | RNA (95) | none |
| `JAYAVELU_2022_DISCOVERY` | Jayavelu et al. 2022 Cancer Cell, discovery cohort | 177 / 177 | protein, label-free, unimputed (177) | OS |
| `JAYAVELU_2022_VALIDATION` | Jayavelu et al. 2022 Cancer Cell, validation cohort | 75 / 75 | protein, SILAC, imputed (75) | OS |
| `LEUCEGENE_SURFACEOME_2024` | Leucegene AML surfaceome, Cell Reports 2024 | 100 / 100 | surface protein (100) | none |
| `ZENODO_14160480` | AML proteogenomics, Zenodo record 14160480 (374 newly diagnosed) | 374 / 374 | RNA (361), protein (374), phosphosites (217), drug IC50, 41 drugs (102) | OS, EFS |

- **Annotation** (`schema/sample_annotation.tsv`): population, age, sex, OncoTree code, tissue,
  timepoint (diagnosis / relapse / refractory / remission / residual), risk group with its classification
  (ELN 2017, ELN 2010 or cytogenetic), KMT2A rearrangement, induction response, OS / EFS in days.
- **Mutations:** `mut_<GENE>` for the same 20 AML genes as the cell lines; NA when the gene was not
  assayed for that sample (Lee 2018 has NPM1 and FLT3 only). Beat AML: WES calls plus clinical FLT3-ITD
  and NPM1; Wang 2022: panel variants at VAF ≥ 5% plus clinical FLT3 / NPM1 / CEBPA; Malani 2022: binary
  mutation table plus clinical FLT3-ITD / TKD; St. Jude 2024: somatic WGS / WES calls (665 samples)
  plus FLT3-ITD from the structural variant calls.
- **MYC / mTOR RNA signatures:** mean z-score of the genes in `data/*_signature_RNA.tsv` within each
  cohort; `MYC_high` / `mTOR_high` = above the cohort median. This reproduces the published cell-line
  classification closely but not exactly (r = 0.90 MYC, 0.82 mTOR on the 56 DepMap lines).
- **Drugs:** `features/drug_*.tsv` adds `compound_key` (normalized name, matches the same drug across
  cohorts), HGNC `targets` from DepMap `PortalCompounds.csv` (name or synonym match), `mechanism`,
  `drug_class`. `lower_is_stronger` is TRUE for AUC and FALSE for DSS / sDSS.
- **Genes:** current HGNC symbols. Older names are mapped to the approved symbol (previous symbol, else a
  unique alias; e.g. HDGFRP3 → HDGFL3) with `data/reference/hgnc_complete_set.txt`, so the same gene has the
  same name in every cohort. RNA keeps HGNC protein-coding genes; duplicated probes / IDs are averaged.
- **Proteomics:** `protein` (gene-level log2 ratio) and `phospho` (single-site log2 values; `features/phospho.tsv`
  adds `site` and `sequence_window`). CPTAC: data freeze v8 `prot_ratio_gene_*_RBC_corrected` and
  `phospho_ratio_single-site_*_RBC_corrected` (not imputed). Acetyl, glyco and succinyl sites are not used yet.
- **PTRC and Beat AML:** PTRC `sample_id` is the Beat AML labId. The 167 specimens in Beat AML waves 1–4 take
  annotation, survival and drug AUC from the matching Beat AML sample (`beat_aml_sample_id`); RNA is not
  repeated. The PTRC drug layer has `duplicate_of = BEAT_AML/drug_auc` in `datasets.tsv` so analyses can
  avoid counting it twice.
- **Zenodo 14160480:** protein values equal to a protein's repeated minimum (fill values for missing
  measurements, about 11% of values) are set to NA; phosphosite intensities are log2-transformed; RNA is
  log2(TPM+1). The 217 MB RNA workbook is streamed with `perl` (too large for the xml2 reader). Drug
  response is log10 IC50 (nM; lower = more sensitive); drug-combination CI values are not used. Clinical
  data from Supp. Data 1 (sex, age, FAB, ELN 2017, WHO 2022, fusions, induction, OS / EFS in days, proteomic
  subtype S1-S8); there are no per-gene mutation calls, so `mut_<GENE>` is NA.
- **St. Jude 2024:** cBioPortal datahub study `aml_stjude_2024`; RNA is log2 CPM as provided. KMT2A-r from
  the molecular category or a KMT2A fusion (partner from the SV file). Extra groupings: molecular category
  (KMT2Ar, RUNX1::RUNX1T1, CBFB::MYH11, NUP98r, UBTF, ...) and FAB. No risk group, survival or drug
  response. 157 samples come from TARGET, and some patients may overlap Wang 2022 (IDs cannot be matched).
- **Jayavelu 2022 (discovery, validation):** Cancer Cell supplementary Data S1 / S3. Discovery protein = the unimputed
  label-free (DIA) proteome, 0 = not quantified -> NA; validation protein = the authors' imputed SILAC ratios (no
  unimputed version). Mutations: targeted myeloid panels (pathological variants; discovery panel II only `mutated`) at
  VAF >= 5% plus clinical NPM1 / FLT3-ITD / TKD; genes with no call in a cohort = not assayed. KMT2A-r from balanced
  11q23 translocations in the karyotype. Extra groupings: AML type (de novo / sAML / tAML), FAB, FLT3-ITD. OS in
  days; relapse-free survival kept as `relapse_free_survival_days` / `rfs_event`.
- **Leucegene surfaceome (Cell Reports 2024):** Tables S1 / S2. Layer `surface` (modality `surface_proteomics`, kept
  apart from whole-cell `protein`): vsn-normalised, knn-imputed log2 intensity as provided; features carry `uniprot`,
  `spat_score` and `likely_surface` (SPAT >= 8). Three HLA-DRB1 accessions listed twice are averaged (higher SPAT
  kept). Mutations tested: NPM1, DNMT3A, RUNX1 (mut_*) and FLT3-ITD (`FLT3_ITD`; TKD not tested, so mut_FLT3 is NA).
  Risk = cytogenetic risk. Extra groupings: AML subgroup, WHO 2022, FAB, FLT3-ITD. 85 specimens are in the AMLmap
  Leucegene RNA cohort, whose participant_id drops AMLmap's leading X so the two cohorts share patients.
- **AMLmap (TCGA, TARGET, Leucegene, LUMC):** one compendium (`data/patients/AMLmap`), split into one cohort per
  study; its Beat AML samples are not used. RNA = log2(CPM+1) of the counts AMLmap corrected for study and
  library batch. Mutations = protein-altering calls at VAF >= 5% plus the FLT3-ITD / TKD, NPM1, IDH1/2 and
  DNMT3A R882 hotspot columns; a gene with no call in a study is treated as not assayed (e.g. KRAS, SRSF2,
  BCOR in LUMC; WT1 in Leucegene). KMT2A from the AMLmap KMT2A diagnosis. Extra groupings: AMLmap cluster,
  WHO 2022, ICC 2022, FAB; `*_like_score` columns hold the AMLmap cell-type scores. OS (days) for TCGA and
  TARGET. TARGET patients may also be in the St. Jude cohort.

## Deconvolution (`prepare_deconvolution.py`)

Estimates the cell type and cell state composition of bulk RNA samples with the tool in `deconvolution/`
(R, InstaPrism and NNLS; see `deconvolution/README.md`; a copy of the jamboree repo's tool, plus `--bulk-input`) against the
AML scAtlas reference, and adds the proportions to `results/` as ordinary layers (modality `deconvolution`;
`docs/results_contract.md`, Deconvolution). **Run it last**: the two `prepare_*.R` scripts wipe their `results/` subfolders.

```bash
uv run prepare_deconvolution.py                       # all sources
uv run prepare_deconvolution.py --only BEAT_AML       # one or more sources
```

| source | bulk input | what the script does with it |
|---|---|---|
| `WANG_2022_PED` (48 samples) | `results/patients/WANG_2022_PED/values/expression.tsv.gz`: recorded unit `log2 CPM` (Table S13 "RNA-seq CPM", described as CPM in the paper; minimum is exactly `log2(0.001)`, so `log2(CPM + 0.001)`), protein-coding genes, HGNC symbols | invert with `2^x - 0.001`; map symbols to Ensembl IDs with `data/reference/hgnc_complete_set.txt`; rescale to 1e6; `--bulk-input normalized`; `caution` warning for the converted scale |
| `BEAT_AML` (AMLmap's Beat AML samples) | `data/patients/AMLmap/counts/uncorrected_counts.csv`: raw integer counts, versioned Ensembl IDs | subset to the cohort's RNA samples (`BA2872R` is `BA2872` in AMLmap), strip version suffixes; `--bulk-input counts` |
| `MALANI_2022` (163 AML samples) | `data/patients/Functional_Precision_Medicine_Tumor_Board_AML/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv`: raw read counts, Ensembl IDs (the 4 healthy CD34+ controls are not in the cohort) | subset to the cohort's samples; `--bulk-input counts` |
| `AML_STJUDE_2024` (887) | `results/patients/AML_STJUDE_2024/values/expression.tsv.gz`: recorded unit `log2 CPM` (cBioPortal `data_mrna_seq_cpm.txt`); the paper (Umeda et al. 2024) transformed counts to log2(CPM) with voom, then corrected library strand and the TARGET batch with ComBat, so values reach -13.6 | invert with `2^x`; map symbols to Ensembl IDs with HGNC; rescale to 1e6; `--bulk-input normalized`; `caution` warning that the values are only approximately CPM (the ComBat step cannot be undone) |
| `TCGA_LAML`, `TARGET_AML`, `LEUCEGENE`, `LUMC_AML` | the same AMLmap counts, one cohort per study (the cohorts `prepare_patient_data.R` makes from AMLmap) | as Beat AML; one run per study |

AMLmap's own processing in `prepare_patient_data.R` uses the batch-corrected counts for the RNA layer; deconvolution uses the raw (`uncorrected_counts.csv`) counts of the same samples. AMLmap's counts are re-processed, so its Beat AML counts are not identical to the earlier biodev file (Spearman about 0.90 per sample).

**Input requirement.** Raw or estimated counts, or linear CPM (library-size normalized only). **Transcript-length-normalized data (TPM, RPKM, FPKM) must not be used**: the reference is built from single-cell counts, which have no gene-length bias, and dividing the bulk by transcript length changes the relative weight of genes. Log-scale values must be inverted to linear first, and only from a CPM unit. A source is judged by the unit recorded for its data (`datasets.tsv`, or the file's documentation): `read_log_expression()` in `prepare_deconvolution.py` accepts only the units in `UNIT_TO_LINEAR` (log2 CPM) and refuses TPM, RPKM and FPKM; the tool has `--bulk-input counts|normalized` (normalized = CPM-like) and cannot detect a wrong unit itself.

Not deconvolved: Zenodo 14160480 (its RNA is `log2(TPM+1)`, TPM, and no counts are available), Wang's raw GEO matrix `GSE192638_Gene_count_matrix_updated.txt` (identical to Table S13, which the paper describes as CPM; the file is named "count matrix", so it is not used directly), cell lines (one line is one cell type), CPTAC and PTRC (no usable RNA; PTRC duplicates Beat AML), Lee 2018
(microarray, deferred), 
Beat AML's filtered `beatAML.rpkm.RDS` was tried first and gave implausible compositions; raw counts fixed that, from the biodev file and now from AMLmap.

Per cohort it writes four layers (InstaPrism initial fit at cell type and cell state, InstaPrism updated fit at cell type,
NNLS at cell type), `deconvolution_runs.tsv` and `deconvolution_warnings.tsv`. Warnings flag, never remove, results:
low initial malignant proportion (the updated fit may be unreliable), hemoglobin-dominated samples, low share of the
reference's genes in the bulk matrix, converted input scale (Wang), and warnings from the tool. Thresholds are constants at the top
of the script. `results/single_cell/aml_scatlas/` gets `reference.tsv`, `markers.tsv` and `provenance.tsv`; the tool's
run logs, bulk matrices and outputs stay in `deconvolution/runs/` (gitignored).

- **Setup:** `renv::restore()` in `deconvolution/` (R 4.6.2; system libraries for some packages, e.g. `libpng-dev`,
  `libuv1-dev`, `libharfbuzz-dev`), `data/reference/aml-scatlas.rds` and `data/reference/hgnc_complete_set.txt` (`data/README.md`).
- **Reference:** everything specific to the reference (ID, `.rds` path, malignant label, citation, limitations) is in one
  block at the top of the script, so a new reference means rebuilding the `.rds`, editing that block and rerunning.
- **Adding a cohort:** a `Source` entry and a reader that returns the bulk matrix (Ensembl IDs x the folder's `sample_id`s)
  in `prepare_deconvolution.py`. `read_log_expression()` already handles a harmonized log layer whose unit is in `UNIT_TO_LINEAR`.
- **Reference limitation:** the scAtlas malignant profile under-represents monocytic and more committed AML.
