# Inputs

```
data/
  cell_lines/      read by prepare_cell_line_data.R
  patients/        read by prepare_patient_data.R             (one folder per cohort)
  single_cell/     read by prepare_single_cell_RNAseq_data.R  (planned, one folder per DATASET_ID)
```

## Committed (small, curated)

| file | content |
|---|---|
| `cell_lines/fig2_MYC_mTOR_DepMap_all_cell_line_classification.tsv` | MYC/mTOR RNA scores and groups, curated KMT2A status (56 lines) |
| `mito_classification.tsv` | mitochondrial High/Low group per line |
| `MYC_signature_RNA.tsv`, `mTOR_signature_RNA.tsv` | MYC / mTOR RNA signature genes (from AML_cell_lines) |

## Raw downloads (gitignored, ~4.3 GB). Place under `data/cell_lines/`:

| folder | files used | source |
|---|---|---|
| `DepMap_26Q3/` | `Model.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `CRISPRGeneEffect.csv`, `OmicsSomaticMutationsMatrixHotspot.csv`, `OmicsSomaticMutationsMatrixDamaging.csv`, `OmicsFusionFiltered.csv`, `PortalCompounds.csv` | DepMap portal, release 26Q3 |
| `DEMETER2_v6/` | `D2_combined_gene_dep_scores.csv` | DepMap portal, DEMETER2 Data v6 |
| `DepMap_Portal_Harmonized_Drug_Screens/` | `PRISM/REPURPOSINGAUCMatrix.csv`, `CTD2/CTRPAUCMatrix.csv`, `GDSC/GDSC2AUCMatrix.csv` | DepMap portal, harmonized drug screens |
| `CCLE_2019/` | `CCLE_miRNA_20181103.gct`, `CCLE_GlobalChromatinProfiling_20181130.csv` | DepMap portal, CCLE 2019 |
| `methylation_RRBS/` | `CCLE_RRBS_TSS_1kb_20180614.txt`, `CCLE_RRBS_TSS_CpG_clusters_20180614.txt` | DepMap portal, CCLE 2019 |
| `targeted_metabolomics/` | `Supp_Tables.xlsx` | Li et al. 2019 Nat Med, supplementary tables |

Downloaded but not yet prepared: `GSE270494/` (EPIC methylation) and `mitoAML/` (cell-line proteome).
Symlinks to an existing copy work too (`ln -s /path/to/DepMap_26Q3 data/cell_lines/`).

## Patient cohorts (gitignored). Place under `data/patients/`:

| folder | DATASET_ID | files used | source |
|---|---|---|---|
| `BeatAML/` | `BEAT_AML` | `beataml_wv1to4_clinical.xlsx`, `beataml_waves1to4_norm_exp_dbgap.txt`, `beataml_probit_curve_fits_v4_dbgap.txt`, `beataml_wes_wv1to4_mutations_dbgap.txt`, `beataml_waves1to4_sample_mapping.xlsx`, `beataml.pvc.csv` (Primitive / Committed cell state) | Bottomly et al. 2022 Cancer Cell; biodev.github.io/BeatAML2 |
| `Nat_Comm_2018/` | `LEE_2018` | `Supp_Data_1.xlsx` (expression, AUC), `Supp_Data_2.xlsx` (drug panel), `Supp_Data_3.xlsx` (clinical) | Lee et al. 2018 Nat Commun, supplementary data |
| `Blood_Cancer_Discovery_2022/` | `WANG_2022_PED` | `Supp_Table.xlsx` (S2 clinical, S3 samples, S4 drugs, S5 drug profile, S9–S10 panel, S13 RNA CPM) | Wang et al. 2022 Blood Cancer Discov; GSE192638 |
| `Functional_Precision_Medicine_Tumor_Board_AML/` | `MALANI_2022` | `File_0` (samples), `File_1.1` (clinical), `File_2` (drug library), `File_3.2` (DSS / sDSS), `File_6` (mutations), `File_7` (RNA log2 CPM) | Malani et al. 2022 Cancer Discov, supplementary files |

| `PTRC/` | `PTRC_BEAT_AML` | `PTRC_metadata.tsv`, `PTRC_protein_gene_ratio_MD.tsv`, `PTRC_phospho_single_site_ratio_MD.tsv` | OHSU / PNNL PTRC (Pino et al. 2024) |
| `zenodo_14160480/` | `ZENODO_14160480` | `Protein abundance of 374 AML patients after adjusting batch effects.csv`, `Phosphopeptides of 217 AML patients.csv`, `Gene expression matrix of 361 AML patients ... (TPM) format.xlsx`, `Supp_Data_1.xlsx` (clinical), `Supp_Data_10.xlsx` (drugs), `Supp_Data_11.xlsx` (IC50) | Zenodo record 14160480 |
| `adult_AML_data_freeze/` | `CPTAC_AML_ADULT` | `proteogenomics_metadata_20250112.tsv`, `v8/prot_ratio_gene_MD_tidy_mapFile_proteogenomics_RBC_corrected.tsv`, `v8/phospho_ratio_signle-site_MD_tidy_mapFile_proteogenomics_RBC_corrected.tsv` | CPTAC adult AML data freeze v8 |

`prepare_patient_data.R` also reads `cell_lines/DepMap_26Q3/PortalCompounds.csv` for drug targets.
Not used: PTRC RNA files (the RNA is in Beat AML), data freeze PTM (acetyl, glyco, succinyl) and imputed / `all` files, `beatAML.rpkm.RDS`, Wang `GSE192638_Gene_count_matrix_updated.txt(.gz)`, Malani `File_3.1`,
`File_4`, `File_5`, `File_8`, `File_9`.
