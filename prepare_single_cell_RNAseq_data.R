# Prepare AML single-cell RNA-seq reference data for the portal.   STATUS: planned (jamboree)
#
# Used by built-in analysis 3 (deconvolution of bulk transcriptomes into AML differentiation
# stages, Karakaslar et al. 2024), with the Childhood Cancer Data Lab (Alex's Lemonade Stand
# Foundation).
#
#   data/single_cell/<DATASET_ID>/   raw downloads per study
#   ->  results/single_cell/<DATASET_ID>/
#         cells.tsv.gz             one row per cell: cell_id, sample_id, cell_type, differentiation stage
#         samples.tsv              sample annotation, fields from schema/sample_annotation.tsv
#         pseudobulk.tsv.gz        rows = genes (HGNC), columns = sample x cell type
#         signature_matrix.tsv     rows = genes, columns = cell types (deconvolution reference)
#         manifest.tsv, provenance.tsv
#   and the deconvolution estimates for bulk cohorts:
#     results/patients/<DATASET_ID>/deconvolution.tsv, results/cell_lines/deconvolution.tsv
#
# Usage: Rscript prepare_single_cell_RNAseq_data.R     (from the repo root)

message("prepare_single_cell_RNAseq_data.R: not implemented yet; see the plan at the top of this file.")
quit(status = 0)
