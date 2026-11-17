#!/usr/bin/env Rscript
# Summarize a reference (from build_reference.R) for the portal's Reference tab.
#
# Writes to --output-dir:
#   reference.tsv   level, label, cell_type, n_cells, n_reference_samples, is_malignant
#   markers.tsv     level, label, rank, gene_id, score (top --n-markers per label,
#                   by the rule select_markers() uses; gene symbols are added later)
#   build.tsv       key, value: what the reference was built from
#
# Example:
#   Rscript summarize_reference.R --reference ../data/reference/aml-scatlas.rds \
#     --malignant-cell-type Malignant --output-dir /tmp/ref_summary
#
# Added in AML_data_processing; not part of the jamboree copy of the tool.

suppressPackageStartupMessages({
  library(optparse)
  library(SummarizedExperiment)
})

module_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", file_arg)))
})
for (f in list.files(file.path(module_dir, "R"), "\\.R$", full.names = TRUE)) {
  source(f)
}

opts <- parse_args(OptionParser(
  option_list = list(
    make_option("--reference", help = "Reference .rds from build_reference.R [required]"),
    make_option("--output-dir", dest = "output_dir", help = "Output directory [required]"),
    make_option("--malignant-cell-type", dest = "malignant_cell_type", help = "Malignant cell type label"),
    make_option("--min-cells", dest = "min_cells", type = "integer", default = 10L, help = "As in deconvolve.R [default: %default]"),
    make_option("--n-markers", dest = "n_markers", type = "integer", default = 20L, help = "Markers per label [default: %default]"),
    make_option("--exclude-gene-categories", dest = "exclude_gene_categories", default = "mitochondrial,ribosomal,sex", help = "As in deconvolve.R [default: %default]"),
    make_option("--gene-categories", dest = "gene_categories", default = file.path(module_dir, "reference-files", "gene_categories.tsv"), help = "Gene category table")
  )
))
for (required in c("reference", "output_dir")) {
  if (is.null(opts[[required]])) stop("--", gsub("_", "-", required), " is required", call. = FALSE)
}
dir.create(opts$output_dir, showWarnings = FALSE, recursive = TRUE)

reference <- load_reference(opts$reference)
exclude <- strsplit(opts$exclude_gene_categories, ",")[[1]]
categories <- read_gene_categories(opts$gene_categories)
excluded_ids <- unique(categories$gene_id[categories$category %in% exclude])

label_rows <- list()
marker_rows <- list()
for (level in c("cell_type", "cell_state")) {
  if (level == "cell_state" && !has_cell_states(reference)) next
  profiles <- prepare_profiles(reference, level == "cell_state", opts$min_cells)
  labels <- profiles$labels
  if (level == "cell_type") {
    labels <- data.frame(
      label = labels$cell_type, cell_type = labels$cell_type,
      n_cells = labels$n_cells, n_reference_samples = labels$n_samples
    )
    counts <- profiles$counts
    colnames(counts) <- labels$label
  } else {
    labels <- data.frame(
      label = labels$cell_state, cell_type = labels$cell_type,
      n_cells = labels$n_cells, n_reference_samples = labels$n_samples
    )
    counts <- profiles$counts
  }
  # a cell type with a single state of the same name appears once, at both levels
  labels$level <- level
  labels$is_malignant <- if (is.null(opts$malignant_cell_type)) NA else labels$cell_type == opts$malignant_cell_type
  label_rows[[level]] <- labels[, c("level", "label", "cell_type", "n_cells", "n_reference_samples", "is_malignant")]

  counts <- counts[!rownames(counts) %in% excluded_ids, , drop = FALSE]
  ranked <- rank_markers(counts, opts$n_markers)
  marker_rows[[level]] <- data.frame(
    level = rep(level, nrow(ranked)), label = ranked$label,
    rank = ranked$rank, gene_id = ranked$gene, score = round(ranked$score, 3)
  )
}
write_tsv <- function(x, name) {
  write.table(x, file.path(opts$output_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
write_tsv(do.call(rbind, label_rows), "reference.tsv")
write_tsv(do.call(rbind, marker_rows), "markers.tsv")

build <- S4Vectors::metadata(reference)$build
as_text <- function(x) if (length(x) == 0) "" else paste(x, collapse = ";")
write_tsv(data.frame(
  key = c("reference_name", "build_date", "cell_type_column", "cell_state_column",
          "sample_column", "unlabelled", "exclude_labels", "n_input_files", "min_cells"),
  value = c(
    as_text(S4Vectors::metadata(reference)$reference_name), as_text(build$timestamp),
    as_text(build$cell_type_column), as_text(build$cell_state_column),
    as_text(build$sample_column), as_text(build$unlabelled),
    as_text(unlist(build$exclude_labels)), as_text(length(build$input_files)),
    as_text(opts$min_cells)
  )
), "build.tsv")
message("Wrote reference.tsv, markers.tsv and build.tsv to ", opts$output_dir)
