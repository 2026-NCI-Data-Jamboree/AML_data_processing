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
| `DepMap_26Q3/` | `Model.csv`, `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`, `CRISPRGeneEffect.csv`, `OmicsSomaticMutationsMatrixHotspot.csv`, `OmicsSomaticMutationsMatrixDamaging.csv`, `OmicsSomaticMutations.csv` (FLT3-ITD / TKD), `OmicsFusionFiltered.csv`, `PortalCompounds.csv` | DepMap portal, release 26Q3 |
| `DEMETER2_v6/` | `D2_combined_gene_dep_scores.csv` | DepMap portal, DEMETER2 Data v6 |
| `DepMap_Portal_Harmonized_Drug_Screens/` | `PRISM/REPURPOSINGAUCMatrix.csv`, `CTD2/CTRPAUCMatrix.csv`, `GDSC/GDSC2AUCMatrix.csv` | DepMap portal, harmonized drug screens |
| `CCLE_2019/` | `CCLE_miRNA_20181103.gct`, `CCLE_GlobalChromatinProfiling_20181130.csv` | DepMap portal, CCLE 2019 |
| `methylation_RRBS/` | `CCLE_RRBS_TSS_1kb_20180614.txt`, `CCLE_RRBS_TSS_CpG_clusters_20180614.txt` | DepMap portal, CCLE 2019 |
| `targeted_metabolomics/` | `Supp_Tables.xlsx` | Li et al. 2019 Nat Med, supplementary tables |
| `surfaceome_George_2026/` | `media-4_surface_mean_abundance.xlsx` (surface proteome; media-2/3 and the paper PDF also kept here) | George et al. 2026 bioRxiv, https://doi.org/10.64898/2026.09.16.752145, supplementary files media-2/3/4 |
| `annotation/` (optional) | `uniprot_cellmembrane.tsv`, `uniprot_sp_tm.tsv`; `Homo_sapiens.gene_info.gz` | UniProt, reviewed human (surface-annotated proteins); NCBI gene_info (chromosome and cytoband of each gene, for the copy-number layer) |
| `DepMap_26Q3/` (copy number) | `OmicsCNGeneWGS.csv` | DepMap portal, release 26Q3: gene-level relative copy number from WGS, for the `cn` layer (an older `OmicsCNGene.csv`, 26Q3 or `DepMap_24Q4/`, is used only if this file is absent) |

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
| `aml_stjude_2024/` | `AML_STJUDE_2024` | `data_clinical_sample.txt`, `data_clinical_patient.txt`, `data_mutations.txt`, `data_sv.txt`, `data_gene_panel_matrix.txt`, `data_mrna_seq_cpm.txt` | cBioPortal datahub, Umeda et al. 2024 Nat Genet |
| `AML_cancer_cell_2022/` | `JAYAVELU_2022_DISCOVERY`, `JAYAVELU_2022_VALIDATION` | `Data S1 Discovery Cohort.xlsx` (clinical, panels I / II, unimputed proteome), `Data S3 Validation Cohort.xlsx` (clinical, panel, imputed proteome) | Jayavelu et al. 2022 Cancer Cell, supplementary Data S1 / S3 |
| `AML_surfaceome_Cell_Reports_2024/` | `LEUCEGENE_SURFACEOME_2024` | `Table_S1.xlsx` (AML surfaceome cohort sheet), `Table_S2.xlsx` (Mass spect data per specimen) | Leucegene AML surfaceome, Cell Reports 2024, supplementary tables |
| `AMLmap/` | `TCGA_LAML`, `TARGET_AML`, `LEUCEGENE`, `LUMC_AML` | `meta.csv`, `scores.csv`, `counts/corrected_counts.csv`, `aberrations/mutations.csv`, `aberrations/aberrations_oh.csv` (Beat AML samples not used) | AMLmap compendium |
| `adult_AML_data_freeze/` | `CPTAC_AML_ADULT` | `proteogenomics_metadata_20250112.tsv`, `v8/prot_ratio_gene_MD_tidy_mapFile_proteogenomics_RBC_corrected.tsv`, `v8/phospho_ratio_signle-site_MD_tidy_mapFile_proteogenomics_RBC_corrected.tsv` | CPTAC adult AML data freeze v8 |

`prepare_patient_data.R` also reads `cell_lines/DepMap_26Q3/PortalCompounds.csv` for drug targets and
`reference/hgnc_complete_set.txt` for current gene symbols (HGNC complete set, from
https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt; gitignored).
Not used: PTRC RNA files (the RNA is in Beat AML), data freeze PTM (acetyl, glyco, succinyl) and imputed / `all` files, `beatAML.rpkm.RDS` (a filtered 15,814-gene RPKM file; it gave implausible deconvolution results), Malani `File_3.1`,
`File_4`, `File_5`, `File_9` (`File_8` is used for deconvolution, below), St. Jude z-scores, GISTIC and ancestry files (St. Jude's `data_mrna_seq_cpm.txt` is deconvolved from the harmonized layer in `results/`, as approximate CPM: voom log2(CPM), then ComBat-corrected), AMLmap's `normalised_and_corrected_counts.csv` and its drug response (Beat AML only).

## Inputs for `prepare_deconvolution.py`

Bulk RNA for deconvolution, in addition to the harmonized `results/`:

| File | Used for |
|---|---|
| `patients/Blood_Cancer_Discovery_2022/GSE192638_Gene_count_matrix_updated.txt` | Not used. Despite its name it is not a count matrix: its values equal Wang's Table S13 (recorded as CPM; the harmonized layer is `log2(CPM + 0.001)`), which is what `prepare_deconvolution.py` reads, from `results/`. Confirmed CPM by the paper (Wang et al. 2022: "counts per million mapped reads (CPM) values ... Supplementary Table S13"; Partek Flow gene-level counts). TPM, RPKM and FPKM must not be deconvolved |
| `patients/Functional_Precision_Medicine_Tumor_Board_AML/File_8_RNA_seq_Raw_Reads_163S_4Healthy.csv` | Malani 2022: raw integer read counts, 60,619 Ensembl genes, 163 AML samples plus 4 healthy CD34+ controls; used for deconvolution (the harmonized RNA layer is log2 CPM from `File_7`) |
| `patients/AMLmap/counts/uncorrected_counts.csv`, `patients/AMLmap/meta.csv` | AMLmap (added 2026-11-17; zip `AMLmap-...zip`): raw integer counts of 1,224 samples from five studies re-processed with one pipeline (60,660 genes, versioned Ensembl IDs; Beat AML 435, Leucegene 399, TCGA 150, TARGET 145, 100LUMC 95). Used for Beat AML and the four temporary `AMLMAP_*` cohorts. Not used: its batch-corrected and normalised count files (not extracted), `scores.csv`, `aberrations/`, `drug_response/` (kept for the upstream processing) |
| `patients/BeatAML/beataml_waves1to4_counts_dbgap.txt` | Beat AML raw counts from the public biodev `beataml2.0_data` repository (https://github.com/biodev/beataml2.0_data, CC-BY-4.0, Git LFS 122 MB). Used for the first Beat AML run and the blast-percent check; **replaced by the AMLmap counts** so that all cohorts but Wang and Zenodo come from one pipeline |
| `reference/aml-scatlas.rds` | The deconvolution reference (433 MB, gitignored). Built with `deconvolution/references/aml-scatlas/` (three commands, see its README; copied from the jamboree repo's `analyses/aml-scatlas-reference`, where the copy in use was built on 2026-11-17) from the AML scAtlas on CELLxGENE (Whittle et al., eLife 2026, https://doi.org/10.7554/eLife.104978; collection https://cellxgene.cziscience.com/collections/071b706a-7ea7-47a4-bddf-6457725839fc, dataset version `be1f3ced-6b8f-47b0-a3e2-7edc8d466abd`; not the authors' figshare copy) with `deconvolution/build_reference.R`. md5 of the copy used: `3d5f6d49185321319e2bc30ce0a9a245` |
| `reference/hgnc_complete_set.txt` | HGNC complete set, also used to map gene symbols to Ensembl IDs (and for marker gene symbols). Record the download date |

Zenodo 14160480 is not deconvolved: it has no raw file and its harmonized RNA layer is `log2(TPM+1)`, i.e. transcript-length-normalized, which must not be used. TPM, RPKM and FPKM are never inputs to deconvolution; counts and CPM are.
