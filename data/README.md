# Inputs

```
data/
  cell_lines/      read by prepare_cell_line_data.R
  patients/        read by prepare_patient_data.R             (planned, one folder per DATASET_ID)
  single_cell/     read by prepare_single_cell_RNAseq_data.R  (planned, one folder per DATASET_ID)
```

## Committed (small, curated)

| file | content |
|---|---|
| `cell_lines/fig2_MYC_mTOR_DepMap_all_cell_line_classification.tsv` | MYC/mTOR RNA scores and groups, curated KMT2A status (56 lines) |
| `mito_classification.tsv` | mitochondrial High/Low group per line |

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
