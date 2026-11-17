# Prepare cell type / cell state proportions of bulk RNA samples for the portal.
#
#   results/ + data/  ->  deconvolution layers added to results/patients/<DATASET_ID>/ and results/cell_lines/
#         values/deconv_<reference>_<method>_<level>.tsv.gz   rows = cell types or states, columns = sample key
#         features/deconv_*.tsv, rows in datasets.tsv (modality deconvolution, unit proportion), manifest.tsv
#         deconvolution_runs.tsv       one row per layer
#         deconvolution_warnings.tsv   flags on runs, samples and labels (never hides a result)
#       results/single_cell/<REFERENCE_ID>/{reference,markers,provenance}.tsv
#
# The deconvolution tool is in deconvolution/ (own renv; launched as a separate Rscript process).
# Layout and decisions: docs/results_contract.md (Deconvolution) and the effort's spec.
#
# RUN LAST: prepare_cell_line_data.R and prepare_patient_data.R wipe their results/ subfolders.
#   Rscript prepare_cell_line_data.R
#   Rscript prepare_patient_data.R
#   Rscript prepare_deconvolution.R [--only <DATASET_ID>[,<DATASET_ID>...]]
#
# Run from the repo root. Needs base R + jsonlite, and the tool's renv restored (deconvolution/README.md).
# Inputs: data/reference/aml-scatlas.rds, data/reference/hgnc_complete_set.txt, and the raw files in SOURCES.
#
# Sources done so far: WANG_2022_PED. To do: BEAT_AML, ZENODO_14160480, cell_lines (log inversion and
# symbol -> Ensembl mapping with HGNC), MALANI_2022 (once harmonized).

suppressMessages(library(jsonlite))

# ---- Reference configuration: everything specific to the reference lives here -------------------
# Swapping the reference means rebuilding the .rds and changing this block, then rerunning.
REFERENCE <- list(
  id          = "aml_scatlas",                 # used in layer keys and in results/single_cell/<id>/
  dataset_id  = "AML_SCATLAS",                 # row in dataset_sources.tsv
  rds         = "data/reference/aml-scatlas.rds",
  malignant   = "Malignant",                   # malignant cell type label in the reference
  name        = "AML scAtlas",
  source      = "AML scAtlas AnnData, https://doi.org/10.48420/27269946; built by 2026-nci-data-jamboree/analyses/aml-scatlas-reference",
  citation    = "Whittle et al., eLife 2026, https://doi.org/10.7554/eLife.104978",
  limitations = paste(
    "Bone marrow, diagnosis samples only; Smart-seq2 libraries (Velten, Stetson) dropped because their counts are reads, not UMIs.",
    "The HSPC label is set only for cells whose author label is HSPC, so blasts in monocytic AML are labelled CD14+ Mono or ProMono:",
    "the normal myeloid profiles carry leukaemic signal and the malignant profile under-represents monocytic and more committed AML."
  )
)

# ---- Settings -----------------------------------------------------------------------------------
N_CORES        <- 2
N_MARKERS      <- 20                      # markers per label in markers.tsv
MIN_BULK_SAMPLES        <- 10             # below this InstaPrism's filter lets one sample decide
LOW_MALIGNANT           <- 0.1            # initial malignant proportion below this: updated fit unreliable
HEMOGLOBIN_SHARE        <- 0.10           # HBB, HBA1, HBA2, HBD above this share of a sample's signal
HEMOGLOBIN_GENES        <- c("HBB", "HBA1", "HBA2", "HBD")
MIN_REFERENCE_GENES_PRESENT <- 0.80       # share of the reference's genes present in the bulk matrix
HGNC_PATH      <- "data/reference/hgnc_complete_set.txt"
TOOL_DIR       <- "deconvolution"
RUNS_DIR       <- file.path(TOOL_DIR, "runs")        # run logs, bulk matrices, tool outputs (gitignored)

# One run per cohort (and one for the cell lines). Each run: InstaPrism (initial + updated fit) and NNLS.
METHODS <- list(
  instaprism = c("--method", "instaprism", "--malignant-cell-type", REFERENCE$malignant),
  nnls       = c("--method", "nnls")
)
METHOD_LABEL <- c(
  instaprism                = "InstaPrism initial fit",
  instaprism_update_all     = "InstaPrism updated fit",
  instaprism_update_malignant = "InstaPrism updated fit (malignant profile only)",
  instaprism_update_pooled  = "InstaPrism updated fit (pooled)",
  nnls                      = "NNLS"
)
UPDATE_METHODS <- c("instaprism_update_all", "instaprism_update_malignant", "instaprism_update_pooled")

# ---- Sources ------------------------------------------------------------------------------------
# Each source returns the linear bulk matrix (genes x samples, columns = the folder's sample key,
# rescaled to 1e6 per sample) plus how it was obtained.
read_wang <- function(sample_ids) {
  path <- "data/patients/Blood_Cancer_Discovery_2022/GSE192638_Gene_count_matrix_updated.txt"
  if (!file.exists(path)) stop("Missing ", path, " (data/README.md)")
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  in_samples <- intersect(sample_ids, colnames(x))
  m <- as.matrix(x[, in_samples, drop = FALSE])
  rownames(m) <- x$gene_id
  list(
    matrix     = sweep(m, 2, colSums(m), "/") * 1e6,
    source     = path,
    conversion = "none (linear TPM-like values rescaled to 1e6 per sample)",
    input      = "normalized"
  )
}

SOURCES <- list(
  list(id = "WANG_2022_PED", folder = "results/patients/WANG_2022_PED", key = "sample_id",
       unit = "samples", reader = read_wang)
)

# ---- Helpers ------------------------------------------------------------------------------------
wt <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
rt <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", quote = "")
write_gz <- function(df, path) {
  con <- gzfile(path, "w"); on.exit(close(con))
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
commit <- suppressWarnings(tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
                                    error = function(e) NA))
commit <- if (length(commit)) commit[1] else NA

args <- commandArgs(trailingOnly = TRUE)
only <- NULL
if (length(args)) {
  if (!identical(args[1], "--only") || length(args) != 2) stop("Usage: Rscript prepare_deconvolution.R [--only ID[,ID...]]")
  only <- strsplit(args[2], ",")[[1]]
  unknown <- setdiff(only, vapply(SOURCES, `[[`, "", "id"))
  if (length(unknown)) stop("Unknown or not yet supported source(s): ", paste(unknown, collapse = ", "))
}

if (!dir.exists("results")) stop("results/ is missing; run prepare_cell_line_data.R and prepare_patient_data.R first")
if (!file.exists(REFERENCE$rds)) stop("Reference not found: ", REFERENCE$rds, " (data/README.md)")
if (!file.exists(HGNC_PATH)) stop("HGNC table not found: ", HGNC_PATH, " (data/README.md)")
if (!dir.exists(file.path(TOOL_DIR, "renv", "library"))) stop("Restore the tool's renv first: see ", TOOL_DIR, "/README.md")

hgnc <- read.delim(HGNC_PATH, stringsAsFactors = FALSE, colClasses = "character", quote = "", comment.char = "")
hgnc_date <- max(hgnc$date_modified, na.rm = TRUE)   # newest record in the snapshot
hgnc <- hgnc[hgnc$ensembl_gene_id != "", c("symbol", "ensembl_gene_id")]
hgnc_symbol_of <- setNames(hgnc$symbol, hgnc$ensembl_gene_id)
hgnc_ensembl_of <- setNames(hgnc$ensembl_gene_id, hgnc$symbol)

tool <- function(script, args, log) {
  force(args); force(log)  # resolve relative paths before changing directory
  old <- setwd(TOOL_DIR); on.exit(setwd(old))
  status <- system2("Rscript", c(script, shQuote(args)), stdout = log, stderr = log)
  if (status != 0) stop(script, " failed; see ", log)
}
abs_path <- function(p) {
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  file.path(normalizePath(dirname(p)), basename(p))
}

# ---- Warning rows -------------------------------------------------------------------------------
warning_rows <- function(key = "", scope = "run", sample_id = "", label = "", severity, code, message, source) {
  data.frame(key = key, scope = scope, sample_id = sample_id, label = label, severity = severity,
             code = code, message = message, source = source, stringsAsFactors = FALSE)
}
tool_warning_code <- function(msg) {
  if (grepl("^Initial '.*' proportion is below", msg)) return(c("low_malignant_initial", "caution"))
  if (grepl("not all integers", msg)) return(c("non_integer_counts", "info"))
  if (grepl("appear more than once", msg)) return(c("duplicate_gene_ids", "info"))
  if (grepl("shared by|gene overlap|overlap", msg, ignore.case = TRUE)) return(c("low_gene_overlap", "warning"))
  if (grepl("min-cells|fewer than .* cells|dropp", msg, ignore.case = TRUE)) return(c("labels_dropped", "info"))
  c("tool_warning", "warning")
}

# ---- One source ---------------------------------------------------------------------------------
process_source <- function(src) {
  message("== ", src$id)
  folder <- src$folder
  samples <- rt(file.path(folder, "samples.tsv"))
  sample_ids <- samples[[src$key]]
  bulk <- src$reader(sample_ids)
  message(ncol(bulk$matrix), " of ", length(sample_ids), " samples have RNA; ", nrow(bulk$matrix), " genes")
  run_dir <- file.path(RUNS_DIR, src$id)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  bulk_path <- abs_path(file.path(run_dir, "bulk.tsv"))
  wt(data.frame(gene_id = rownames(bulk$matrix), bulk$matrix, check.names = FALSE), bulk_path)

  # run the tool once per method
  props <- list(); infos <- list(); info_of_method <- list()
  for (m in names(METHODS)) {
    message("running ", m)
    tool("deconvolve.R",
         c("--reference", abs_path(REFERENCE$rds), "--bulk", bulk_path, "--bulk-input", bulk$input,
           METHODS[[m]], "--output-dir", abs_path(run_dir), "--output-prefix", paste0(src$id, "_", m),
           "--n-cores", N_CORES),
         abs_path(file.path(run_dir, paste0(m, ".log"))))
    props[[m]] <- read.delim(file.path(run_dir, paste0(src$id, "_", m, "_proportions.tsv")), stringsAsFactors = FALSE)
    infos[[m]] <- fromJSON(file.path(run_dir, paste0(src$id, "_", m, "_run_info.json")))
    for (mm in unique(props[[m]]$method)) info_of_method[[mm]] <- m
  }
  props <- do.call(rbind, props)

  # drop earlier deconvolution output of this folder
  unlink(list.files(file.path(folder, "values"), "^deconv_", full.names = TRUE))
  unlink(list.files(file.path(folder, "features"), "^deconv_", full.names = TRUE))

  # layers
  combos <- unique(props[, c("method", "level", "quantity")])
  datasets <- rt(file.path(folder, "datasets.tsv"))
  datasets <- datasets[!startsWith(datasets$key, "deconv_"), , drop = FALSE]
  manifest <- rt(file.path(folder, "manifest.tsv"))
  manifest <- manifest[!grepl("deconv", manifest$path), , drop = FALSE]
  n_col <- if (src$unit == "samples") "n_samples" else "n_lines"
  runs <- list(); new_datasets <- list(); new_manifest <- list()
  keys <- character(0); key_method <- character(0)
  for (i in seq_len(nrow(combos))) {
    cm <- combos[i, ]
    key <- paste("deconv", REFERENCE$id, cm$method, cm$level, sep = "_")
    p <- props[props$method == cm$method & props$level == cm$level, ]
    wide <- tapply(p$proportion, list(p$label, p$sample_id), sum)
    wide <- wide[order(rownames(wide)), intersect(sample_ids, colnames(wide)), drop = FALSE]
    write_gz(data.frame(feature = rownames(wide), round(wide, 6), check.names = FALSE),
             file.path(folder, "values", paste0(key, ".tsv.gz")))
    wt(data.frame(feature = rownames(wide), label = rownames(wide), gene_symbol = "", entrez_id = ""),
       file.path(folder, "features", paste0(key, ".tsv")))
    ds <- setNames(as.list(rep("", ncol(datasets))), names(datasets))
    ds$key <- key
    ds$label <- paste0(METHOD_LABEL[[cm$method]], ", ", if (cm$level == "cell_type") "cell types" else "cell states",
                       " (", REFERENCE$name, ")")
    ds$modality <- "deconvolution"; ds$unit <- "proportion"
    ds$source <- paste0(REFERENCE$name, " reference, ", sub(" \\(.*", "", METHOD_LABEL[[cm$method]]))
    ds[[n_col]] <- as.character(ncol(wide))
    ds$n_features <- as.character(nrow(wide))
    new_datasets[[key]] <- as.data.frame(ds, stringsAsFactors = FALSE)
    new_manifest[[key]] <- data.frame(
      path = c(paste0("features/", key, ".tsv"), paste0("values/", key, ".tsv.gz")),
      kind = c("features", "values"),
      description = paste0(ds$label, c(" - feature annotation", " - values per sample")))
    info <- infos[[info_of_method[[cm$method]]]]
    runs[[key]] <- data.frame(
      key = key, reference = REFERENCE$id, method = cm$method, level = cm$level, quantity = cm$quantity,
      instaprism_filter = if (grepl("^instaprism", cm$method)) tolower(info$method$options$instaprism_filter) else "",
      update_mode = if (cm$method %in% UPDATE_METHODS) info$method$options$update_mode else "",
      min_cells = info$parameters$min_cells,
      malignant_cell_type = if (grepl("^instaprism", cm$method)) REFERENCE$malignant else "",
      n_samples = ncol(wide),
      input_source = bulk$source, input_conversion = bulk$conversion,
      tool_version = commit, run_date = format(Sys.Date()), hgnc_snapshot = hgnc_date)
    keys <- c(keys, key); key_method <- c(key_method, cm$method)
  }
  wt(rbind(datasets, do.call(rbind, new_datasets)), file.path(folder, "datasets.tsv"))
  runs <- do.call(rbind, runs)
  wt(runs, file.path(folder, "deconvolution_runs.tsv"))
  names(keys) <- key_method
  update_keys <- keys[key_method %in% UPDATE_METHODS]

  # warnings
  w <- list()
  seen <- list()
  for (m in names(infos)) {
    tw <- infos[[m]]$warnings
    if (!length(tw) || !nrow(tw)) next
    run_keys <- keys[key_method %in% names(info_of_method)[unlist(info_of_method) == m]]
    for (msg in tw$message) {
      cs <- tool_warning_code(msg)
      # the low-malignant warning is about the reference update, so it flags the updated layers only
      k <- if (cs[1] == "low_malignant_initial") intersect(run_keys, update_keys) else run_keys
      if (length(k)) seen[[length(seen) + 1]] <- list(msg = msg, code = cs[1], severity = cs[2], keys = unname(k))
    }
  }
  for (s in seen) {
    w[[length(w) + 1]] <- warning_rows(key = s$keys, severity = s$severity, code = s$code,
                                       message = s$msg, source = "tool")
  }
  if (length(update_keys)) {
    w[[length(w) + 1]] <- warning_rows(
      key = unname(update_keys), severity = "caution", code = "updated_fit_batch_dependent",
      message = paste("The updated fit re-estimates the normal cell type profiles pooled over all samples in this run, so a sample's",
                      "estimate depends on the other samples; it also uses the", REFERENCE$name, "malignant profile, which",
                      "under-represents monocytic and more committed AML. Prefer the initial fit unless the update is needed."),
      source = "processing")
    ini <- props[props$method == "instaprism" & props$level == "cell_type" & props$label == REFERENCE$malignant, ]
    low <- ini[ini$proportion < LOW_MALIGNANT, ]
    if (nrow(low)) {
      w[[length(w) + 1]] <- do.call(rbind, lapply(unname(update_keys), function(k) warning_rows(
        key = k, scope = "sample", sample_id = low$sample_id, severity = "caution", code = "low_malignant_initial",
        message = sprintf("Initial %s proportion %.3f is below %.2f: the updated fit may let the malignant profile absorb normal cell types",
                          REFERENCE$malignant, low$proportion, LOW_MALIGNANT), source = "processing")))
    }
  }
  if (ncol(bulk$matrix) < MIN_BULK_SAMPLES) {
    w[[length(w) + 1]] <- warning_rows(severity = "warning", code = "small_bulk_set", source = "processing",
      message = sprintf("Only %d samples in this run (below %d): InstaPrism's gene filter then lets a single sample decide which genes are dropped",
                        ncol(bulk$matrix), MIN_BULK_SAMPLES))
  }
  hb_ids <- hgnc_ensembl_of[HEMOGLOBIN_GENES]
  hb <- colSums(bulk$matrix[rownames(bulk$matrix) %in% hb_ids, , drop = FALSE]) / colSums(bulk$matrix)
  hb <- hb[hb > HEMOGLOBIN_SHARE]
  if (length(hb)) {
    w[[length(w) + 1]] <- warning_rows(scope = "sample", sample_id = names(hb), severity = "caution",
      code = "hemoglobin_dominated", source = "processing",
      message = sprintf("%s make up %.0f%% of this sample's signal (above %.0f%%); red cell contamination can distort the estimates",
                        paste(HEMOGLOBIN_GENES, collapse = ", "), 100 * hb, 100 * HEMOGLOBIN_SHARE))
  }
  funnel <- infos[[1]]$genes$funnel
  present <- funnel$shared / funnel$reference
  if (present < MIN_REFERENCE_GENES_PRESENT) {
    w[[length(w) + 1]] <- warning_rows(severity = "warning", code = "low_reference_gene_coverage", source = "processing",
      message = sprintf("Only %.0f%% of the reference's genes (%d of %d) are present in the bulk matrix (below %.0f%%)",
                        100 * present, funnel$shared, funnel$reference, 100 * MIN_REFERENCE_GENES_PRESENT))
  }
  w <- if (length(w)) do.call(rbind, w) else warning_rows(severity = character(0), code = character(0), message = character(0), source = character(0))
  w <- unique(w[order(w$key, w$scope, w$sample_id), ])
  wt(w, file.path(folder, "deconvolution_warnings.tsv"))
  message(nrow(w), " warning rows; ", nrow(runs), " layers")

  new_manifest[["files"]] <- data.frame(
    path = c("deconvolution_runs.tsv", "deconvolution_warnings.tsv"), kind = c("deconvolution", "deconvolution"),
    description = c("Deconvolution runs, one row per layer", "Deconvolution warnings on runs, samples and labels"))
  wt(rbind(manifest, do.call(rbind, new_manifest)), file.path(folder, "manifest.tsv"))
  invisible(runs)
}

for (src in SOURCES) if (is.null(only) || src$id %in% only) process_source(src)

# ---- Reference summary --------------------------------------------------------------------------
message("== reference ", REFERENCE$id)
dir.create(RUNS_DIR, recursive = TRUE, showWarnings = FALSE)
ref_tmp <- abs_path(file.path(RUNS_DIR, "reference_summary"))
tool("summarize_reference.R",
     c("--reference", abs_path(REFERENCE$rds), "--malignant-cell-type", REFERENCE$malignant,
       "--n-markers", N_MARKERS, "--output-dir", ref_tmp),
     abs_path(file.path(RUNS_DIR, "reference_summary.log")))
ref_dir <- file.path("results", "single_cell", REFERENCE$id)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(file.path(ref_tmp, "reference.tsv"), ref_dir, overwrite = TRUE))
mk <- read.delim(file.path(ref_tmp, "markers.tsv"), stringsAsFactors = FALSE)
mk$gene_symbol <- unname(hgnc_symbol_of[mk$gene_id])
wt(mk[, c("level", "label", "rank", "gene_id", "gene_symbol", "score")], file.path(ref_dir, "markers.tsv"))
bd <- rt(file.path(ref_tmp, "build.tsv")); build <- setNames(bd$value, bd$key)
wt(data.frame(
  key = c("reference_id", "reference_name", "source", "citation", "build_date", "processing_commit", "hgnc_snapshot",
          "cell_type_column", "cell_state_column", "sample_column", "unlabelled", "excluded_labels", "n_input_files",
          "min_cells", "malignant_cell_type", "limitations"),
  value = c(REFERENCE$id, REFERENCE$name, REFERENCE$source, REFERENCE$citation, build[["build_date"]], commit, hgnc_date,
            build[["cell_type_column"]], build[["cell_state_column"]], build[["sample_column"]], build[["unlabelled"]],
            build[["exclude_labels"]], build[["n_input_files"]], build[["min_cells"]], REFERENCE$malignant,
            REFERENCE$limitations)), file.path(ref_dir, "provenance.tsv"))
message("wrote ", ref_dir)

invisible(file.copy("dataset_sources.tsv", "results", overwrite = TRUE))
message("Done.")
