# Results contract (version 3)

`results/` holds the processed, harmonized data. `AML_data_visualization` reads **only** `results/`.
All tables are tab-separated with a header row and empty cells for NA; `.tsv.gz` files are
gzip-compressed TSVs.

```
results/
  README.md                   overview of the harmonized data (copy of docs/results_README.md)
  catalog/datasets.tsv        dataset inventory (catalog/README.md)
  dataset_sources.tsv         every integrated dataset: source, accession, NCI resource status
  schema/*.tsv                harmonized annotation schema
  cell_lines/                 prepare_cell_line_data.R
  patients/<DATASET_ID>/      prepare_patient_data.R
  single_cell/<REFERENCE_ID>/ prepare_deconvolution.R (planned): deconvolution reference summary
```

Every data folder (`cell_lines/`, `patients/<DATASET_ID>/`) has the same layout:

| file | content |
|---|---|
| `samples.tsv` | one row per sample. First column is the sample key used as column names in `values/` (`cell_line` for cell lines, `sample_id` for patients). Harmonized fields follow `schema/sample_annotation.tsv` (`study_id`, `sample_type`, `participant_id`, `disease_oncotree`, `population`, `age_years`, `sex`, `mutations`, ...) |
| `groupings.tsv` | sample annotation columns offered as two-group comparisons: `column`, `label`, `true_label`, `false_label`, `kind` (`binary` = TRUE/FALSE column, `categorical` = text column). The portal offers exactly these, so a new grouping only needs a column in `samples.tsv` and a row here |
| `datasets.tsv` | one row per omics layer: `key`, `label`, `modality` (`schema/vocabularies.tsv`), `unit`, `source`, `lower_is_stronger` (TRUE = lower value means more dependent / sensitive), `n_features`, `n_lines` (cell lines) or `n_samples` (patients); patients also `duplicate_of` (`<DATASET_ID>/<key>` when the layer repeats another cohort's values, e.g. PTRC drug AUC) |
| `features/<key>.tsv` | one row per feature: `feature` (ID used in `values/`), `label` (display name), `gene_symbol` (HGNC; empty if not a gene), `entrez_id`; drug layers add `compound_name`, `targets` (`;`-separated HGNC), `mechanism` (patients also `compound_key`, the normalized name shared across cohorts, `depmap_compound`, `drug_class`); methylation layers add `region`; phospho layers add `site`, `sequence_window` |
| `values/<key>.tsv.gz` | matrix: column `feature`, then one column per sample |
| `manifest.tsv` | `path`, `kind`, `description` for every file |
| `provenance.tsv` | key/value: `prepared_at`, `processing_commit`, `working_tree_dirty`, `r_version`, `results_version` |

## Deconvolution (additive; `results_version` stays 3)

`prepare_deconvolution.R` runs last (the other scripts wipe their `results/` subfolders) and adds the
following to a data folder (`cell_lines/`, `patients/<DATASET_ID>/`) that has bulk RNA. The portal checks
that a layer exists before showing it.

**Layers.** `values/deconv_<reference>_<method>_<level>.tsv.gz`, e.g. `deconv_aml_scatlas_instaprism_cell_type`.
Rows are cell types or cell states (`feature`), columns are the folder's sample key (as for every other layer),
values are proportions summing to 1 per sample; the exception is the updated fit's cell-state rows, which
cover only the normal cell types. Each method and level a run produces is its own layer (InstaPrism initial
fit at cell type and cell state, InstaPrism updated, NNLS at cell type). The quantity (mRNA share or cell
share) is not in the key; it is in `deconvolution_runs.tsv`. Each layer has `features/<key>.tsv` and a
`datasets.tsv` row with modality `deconvolution`, unit `proportion`, `lower_is_stronger` empty.

| file | content |
|---|---|
| `deconvolution_runs.tsv` | one row per layer: `key`, `reference`, `method`, `level`, `quantity` (`mrna_share`/`cell_share`), batch-dependent settings, tool version, run date, input source and conversion used |
| `deconvolution_warnings.tsv` | `key` (empty = all layers in the folder), `scope` (`run`/`sample`/`label`), `sample_id`, `label` (empty when `scope` is `run`), `severity` (`info`/`warning`/`caution`), `code` (stable short code), `message`, `source` (`tool` copies `run_info.json`, `processing` is added by the script). Warnings flag a result; they never hide or block it |

**`single_cell/<REFERENCE_ID>/`** summarizes the reference the layers were built with:

| file | content |
|---|---|
| `reference.tsv` | one row per cell type and per cell state: `level`, `label`, `cell_type` (parent; equals the label for a cell type row), `n_cells`, `n_reference_samples`, `is_malignant` |
| `markers.tsv` | top 20 genes per label: `level`, `label`, `rank`, `gene_id`, `gene_symbol`, `score` |
| `provenance.tsv` | key/value: reference name, source, citation, build date, tool version, label columns, excluded labels, unlabelled policy, `limitations` |

The reference also gets a dataset ID and a row in `dataset_sources.tsv`.

## `patients/<DATASET_ID>/samples.tsv` survival

The portal's Patients → Survival tab uses, per cohort, one row per participant (the `diagnosis` sample when
`timepoint` is given):

| columns | content |
|---|---|
| `overall_survival_days`, `os_event` | OS time in days; TRUE = death observed, FALSE = censored |
| `event_free_survival_days`, `efs_event` | optional EFS, same coding |

Patients can be grouped by any `groupings.tsv` row, by `risk_group`, `KMT2A_rearranged`,
`treatment_response`, `disease_oncotree`, `population`, `sex`, by each gene in `mutations` with ≥ 5 mutant
patients, or by high vs low values of a gene in any `values/` layer.

## `cell_lines/samples.tsv` groupings

Binary groupings are TRUE/FALSE, and empty when the line is not annotated:

| column | TRUE means | annotated lines |
|---|---|---|
| `mTOR_high` | mTOR RNA signature score above the median | 56 (curated classification) |
| `MYC_high` | MYC RNA signature score above the median | 56 |
| `KMT2A_rearranged` | KMT2A fusion (`KMT2A_partner`; `KMT2A_source` = curated or DepMap fusion calls) | 70 |
| `cn` layer | gene-level relative copy number, log2(relative CN + 1) (1 = ploidy-neutral); features add `chromosome`, `cytoband` (NCBI) so the portal can aggregate genes into arms / bands. Modality `copy_number` | lines with CN |
| `FLT3_ITD` | FLT3 internal tandem duplication (juxtamembrane in-frame insertion, variant-level calls); `mut_FLT3` also counts ITD and TKD (D835/I836/N841/Y842) | lines with mutation calls |
| `mut_<GENE>` | hotspot or damaging mutation in GENE (20 recurrent AML genes) | lines with mutation calls |
| `Mito_high` | mitochondrial High group | lines in `data/mito_classification.tsv` |

Categorical groupings: `MYC_mTOR_category` (`MYC-High / mTOR-Low`, ...), `disease_oncotree`,
`population`, `sex`. The raw scores are in `MYC_score` and `mTOR_score`.

## Rules for changes

- Adding files or columns is backward compatible.
- Renaming or removing anything the portal reads means bumping `results_version` and updating
  the portal at the same time.
