# Results contract (version 3)

`results/` holds the processed, harmonized data. `AML_data_visualization` reads **only** `results/`.
All tables are tab-separated with a header row and empty cells for NA; `.tsv.gz` files are
gzip-compressed TSVs.

```
results/
  catalog/datasets.tsv        dataset inventory (catalog/README.md)
  schema/*.tsv                harmonized annotation schema
  cell_lines/                 prepare_cell_line_data.R
  patients/<DATASET_ID>/      prepare_patient_data.R
  single_cell/<DATASET_ID>/   prepare_single_cell_RNAseq_data.R (planned)
```

Every data folder (`cell_lines/`, `patients/<DATASET_ID>/`) has the same layout:

| file | content |
|---|---|
| `samples.tsv` | one row per sample. First column is the sample key used as column names in `values/` (`cell_line` for cell lines, `sample_id` for patients). Harmonized fields follow `schema/sample_annotation.tsv` (`study_id`, `sample_type`, `participant_id`, `disease_oncotree`, `population`, `age_years`, `sex`, `mutations`, ...) |
| `groupings.tsv` | sample annotation columns offered as two-group comparisons: `column`, `label`, `true_label`, `false_label`, `kind` (`binary` = TRUE/FALSE column, `categorical` = text column). The portal offers exactly these, so a new grouping only needs a column in `samples.tsv` and a row here |
| `datasets.tsv` | one row per omics layer: `key`, `label`, `modality` (`schema/vocabularies.tsv`), `unit`, `source`, `lower_is_stronger` (TRUE = lower value means more dependent / sensitive), `n_features`, `n_lines` (cell lines) or `n_samples` (patients) |
| `features/<key>.tsv` | one row per feature: `feature` (ID used in `values/`), `label` (display name), `gene_symbol` (HGNC; empty if not a gene), `entrez_id`; drug layers add `compound_name`, `targets` (`;`-separated HGNC), `mechanism` (patients also `compound_key`, the normalized name shared across cohorts, `depmap_compound`, `drug_class`); methylation layers add `region` |
| `values/<key>.tsv.gz` | matrix: column `feature`, then one column per sample |
| `manifest.tsv` | `path`, `kind`, `description` for every file |
| `provenance.tsv` | key/value: `prepared_at`, `processing_commit`, `working_tree_dirty`, `r_version`, `results_version` |

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
| `mut_<GENE>` | hotspot or damaging mutation in GENE (20 recurrent AML genes) | lines with mutation calls |
| `Mito_high` | mitochondrial High group | lines in `data/mito_classification.tsv` |

Categorical groupings: `MYC_mTOR_category` (`MYC-High / mTOR-Low`, ...), `disease_oncotree`,
`population`, `sex`. The raw scores are in `MYC_score` and `mTOR_score`.

## Rules for changes

- Adding files or columns is backward compatible.
- Renaming or removing anything the portal reads means bumping `results_version` and updating
  the portal at the same time.
