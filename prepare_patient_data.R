# Prepare harmonized AML patient-cohort data for the portal.
#
#   data/patients/<folder>/   raw downloads per cohort (see data/README.md)
#   ->  results/patients/<DATASET_ID>/      (IDs from catalog/datasets.tsv)
#         samples.tsv              one row per sample: schema/sample_annotation.tsv fields
#                                  (disease, timepoint, risk group, mutations, KMT2A, survival, ...)
#         groupings.tsv            sample annotation columns offered as group comparisons
#         datasets.tsv             one row per layer (RNA expression, drug response, ...)
#         features/<key>.tsv       feature annotation; drug layers add compound_key / targets / mechanism
#         values/<key>.tsv.gz      matrix: rows = features, columns = sample_id
#         manifest.tsv, provenance.tsv
#
# Cohorts with ex vivo drug response (built-in analysis 1):
#   BEAT_AML       Beat AML waves 1-4 (Bottomly et al. 2022)       RNA-seq, WES, 166 inhibitors (AUC)
#   LEE_2018       Lee et al. 2018 Nat Commun                      expression array, 160 drugs (AUC)
#   WANG_2022_PED  Wang et al. 2022 Blood Cancer Discov (pediatric) RNA-seq, panel, 45 drugs (AUC)
#   MALANI_2022    Malani et al. 2022 Cancer Discov (FPM)          RNA-seq, panel, 515 drugs (DSS, sDSS)
# RNA-only cohort (no drug response):
#   AML_STJUDE_2024  Umeda et al. 2024 Nat Genet (pediatric)       RNA-seq, WGS / WES, fusions
#   TCGA_LAML, TARGET_AML, LEUCEGENE, LUMC_AML  AMLmap compendium   RNA-seq (batch-corrected), mutations
#
# Harmonization
#   genes     current HGNC symbols (previous symbols and aliases mapped; data/reference/hgnc_complete_set.txt);
#             RNA keeps protein-coding genes only
#   drugs     compound_key = upper-case name without salts / brackets / punctuation, so the same drug
#             matches across cohorts; HGNC targets from DepMap PortalCompounds (name or synonym match)
#   mutations mut_<GENE> TRUE/FALSE for 20 recurrent AML genes (NA = gene not assayed for that sample)
#   MYC / mTOR RNA signature groups are not derived for patients (cell lines only, prepare_cell_line_data.R)
#
# Usage: Rscript prepare_patient_data.R     (from the repo root; base R + xml2)

suppressMessages(library(xml2))

dd       <- "data/patients"
out_root <- "results/patients"
unlink(out_root, recursive = TRUE)

aml_genes <- c("TP53", "FLT3", "NPM1", "NRAS", "KRAS", "PTPN11", "KIT", "IDH1", "IDH2", "DNMT3A",
               "TET2", "ASXL1", "RUNX1", "CEBPA", "WT1", "SRSF2", "U2AF1", "SF3B1", "STAG2", "BCOR")
DAYS_PER_MONTH <- 30.4375

# ================================================================== helpers ==
xlsx_sheets <- function(path) {
  ns <- c(d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships")
  sh   <- xml_find_all(read_xml(unz(path, "xl/workbook.xml")), "//d1:sheets/d1:sheet", ns)
  rels <- xml_find_all(read_xml(unz(path, "xl/_rels/workbook.xml.rels")), "//*[local-name()='Relationship']")
  tgt  <- setNames(xml_attr(rels, "Target"), xml_attr(rels, "Id"))[xml_attr(sh, "r:id", ns = ns)]
  setNames(sub("^/?(xl/)?", "xl/", tgt), xml_attr(sh, "name"))
}
# one sheet as a data.frame of character columns (header = first row after `skip`)
xlsx_read <- function(path, sheet = 1, header = TRUE, skip = 0) {
  ns <- c(d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
  f  <- xlsx_sheets(path)[[sheet]]
  ss <- if ("xl/sharedStrings.xml" %in% unzip(path, list = TRUE)$Name) {
    vapply(xml_find_all(read_xml(unz(path, "xl/sharedStrings.xml")), "//d1:si", ns), function(si)
      paste(xml_text(xml_find_all(si, ".//d1:t[not(ancestor::d1:rPh)]", ns)), collapse = ""), "")
  } else character(0)
  cells <- xml_find_all(read_xml(unz(path, f)), "//d1:sheetData/d1:row/d1:c", ns)
  ref <- xml_attr(cells, "r"); typ <- xml_attr(cells, "t")
  val <- xml_text(xml_find_first(cells, "d1:v", ns))
  s <- !is.na(typ) & typ == "s"; val[s] <- ss[as.integer(val[s]) + 1]
  inl <- !is.na(typ) & typ == "inlineStr"
  if (any(inl)) val[inl] <- xml_text(xml_find_first(cells[inl], "d1:is", ns))
  col <- gsub("[0-9]", "", ref); row <- as.integer(gsub("[A-Z]", "", ref))
  colnum <- vapply(strsplit(col, ""), function(l) sum(match(l, LETTERS) * 26^(rev(seq_along(l)) - 1)), 0)
  M <- matrix(NA_character_, max(row), max(colnum)); M[cbind(row, colnum)] <- val
  if (skip) M <- M[-seq_len(skip), , drop = FALSE]
  if (!header) return(M)
  h <- trimws(gsub("\\s+", " ", M[1, ])); h[is.na(h) | h == ""] <- paste0("V", which(is.na(h) | h == ""))
  d <- as.data.frame(M[-1, , drop = FALSE], stringsAsFactors = FALSE); names(d) <- make.unique(h)
  d[rowSums(!is.na(d)) > 0, , drop = FALSE]
}
# very large single-sheet xlsx (too big for xml2): stream sheet1 to a TSV with perl, one row per line;
# shared-string cells come out as "s:<index>" and are resolved against sharedStrings.xml
xlsx_big_to_df <- function(path, n_cols) {
  tsv <- tempfile(fileext = ".tsv")
  perl <- paste0(
    'BEGIN{$/="</c>"; $cur=0; @row=()} ',
    'if (/<c r="([A-Z]+)(\\d+)"([^>]*)>(?:.*?<v>([^<]*)<\\/v>)?/s) { ',
    '  ($cl,$rn,$at,$v)=($1,$2,$3,$4); $v="" unless defined $v; ',
    '  if ($rn != $cur) { print join("\\t", map { defined $_ ? $_ : "" } @row[0..', n_cols - 1, ']), "\\n" if $cur; @row=(); $cur=$rn } ',
    '  $c=0; $c=$c*26+ord($_)-64 for split //, $cl; ',
    '  $row[$c-1] = ($at =~ /t="s"/) ? "s:$v" : $v } ',
    'END{ print join("\\t", map { defined $_ ? $_ : "" } @row[0..', n_cols - 1, ']), "\\n" if $cur }')
  status <- system(paste("unzip -p", shQuote(path), "xl/worksheets/sheet1.xml | perl -ne", shQuote(perl), ">", shQuote(tsv)))
  stopifnot(status == 0)
  ns <- c(d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
  ss <- vapply(xml_find_all(read_xml(unz(path, "xl/sharedStrings.xml")), "//d1:si", ns), function(si)
    paste(xml_text(xml_find_all(si, ".//d1:t[not(ancestor::d1:rPh)]", ns)), collapse = ""), "")
  d <- read.delim(tsv, header = FALSE, colClasses = "character", quote = "", comment.char = "", na.strings = "")
  unlink(tsv)
  res <- function(x) ifelse(!is.na(x) & startsWith(x, "s:"), ss[as.integer(sub("^s:", "", x)) + 1], x)
  h <- res(unlist(d[1, ])); d <- d[-1, , drop = FALSE]
  for (j in seq_len(ncol(d))) if (any(startsWith(d[[j]], "s:"), na.rm = TRUE)) d[[j]] <- res(d[[j]])
  names(d) <- make.unique(ifelse(is.na(h) | h == "", paste0("V", seq_along(h)), h))
  d
}
num <- function(x) suppressWarnings(as.numeric(x))
yn  <- function(x) { x <- toupper(trimws(x))
  ifelse(x %in% c("Y", "YES", "1", "TRUE", "+", "POSITIVE"), TRUE,
         ifelse(x %in% c("N", "NO", "0", "FALSE", "-", "NEGATIVE"), FALSE, NA)) }
# rows = features; average duplicated feature IDs (e.g. several array probes per gene)
collapse_rows <- function(X, ids) {
  keep <- !is.na(ids) & ids != ""
  X <- X[keep, , drop = FALSE]; ids <- ids[keep]
  if (!anyDuplicated(ids)) { rownames(X) <- ids; return(X) }
  S <- rowsum(ifelse(is.na(X), 0, X), ids); N <- rowsum(1 * !is.na(X), ids)
  S / ifelse(N == 0, NA, N)
}
num_matrix <- function(d, cols) { X <- as.matrix(d[, cols, drop = FALSE]); storage.mode(X) <- "double"; X }

# ---- gene annotation: current HGNC symbols (data/reference/hgnc_complete_set.txt)
# Cohorts use gene names from different annotation releases (e.g. HDGFRP3, now HDGFL3), so every gene name
# is mapped to its current approved symbol: approved symbol as is, else a previous symbol, else an alias,
# when that name points to a single approved gene. Protein-coding = HGNC locus group "protein-coding gene".
hgnc <- local({
  h <- read.delim("data/reference/hgnc_complete_set.txt", colClasses = "character", quote = "", na.strings = "")
  h <- h[h$status %in% "Approved", ]
  unq <- function(col) {                       # other name (upper case) -> approved symbol, if unique
    l <- strsplit(gsub('"', "", ifelse(is.na(h[[col]]), "", h[[col]])), "\\|")
    d <- data.frame(from = toupper(unlist(l)), to = rep(h$symbol, lengths(l)))
    d <- unique(d[d$from != "" & !d$from %in% toupper(h$symbol), ])
    d <- d[!d$from %in% d$from[duplicated(d$from)], ]
    setNames(d$to, d$from)
  }
  list(symbol = h$symbol, prev = unq("prev_symbol"), alias = unq("alias_symbol"),
       ensembl = setNames(h$symbol, h$ensembl_gene_id)[!is.na(h$ensembl_gene_id)],
       protein_coding = h$symbol[h$locus_group %in% "protein-coding gene"])
})
to_hgnc <- function(x) {
  out <- ifelse(x %in% hgnc$symbol, x, NA_character_)
  i <- is.na(out); out[i] <- hgnc$prev[toupper(x[i])]          # names are matched ignoring case
  i <- is.na(out); out[i] <- hgnc$alias[toupper(x[i])]
  unname(out)
}
# rows renamed to current symbols; rows of genes that now share a symbol are averaged; unknown names kept
# (keep_unmapped = TRUE) or dropped
hgnc_rows <- function(X, keep_unmapped = FALSE) {
  s <- to_hgnc(rownames(X))
  if (keep_unmapped) s <- ifelse(is.na(s), rownames(X), s)
  collapse_rows(X, s)
}
protein_coding <- hgnc$protein_coding
# Ensembl gene ID -> symbol for Beat AML and Malani 2022 (HGNC first, then the Beat AML annotation)
gene_ann <- local({
  b <- read.delim(file.path(dd, "BeatAML/beataml_waves1to4_norm_exp_dbgap.txt"), nrows = 1e5,
                  colClasses = c(rep("character", 2), rep("NULL", 709)), check.names = FALSE)
  setNames(b$display_label, b$stable_id)
})
ensembl_symbol <- function(id) { s <- unname(hgnc$ensembl[id]); ifelse(is.na(s), to_hgnc(unname(gene_ann[id])), s) }

# ---- drug annotation: DepMap compounds (names + synonyms -> HGNC targets)
compound_key <- function(x) {
  x <- toupper(x)
  x <- gsub("\\(.*?\\)", " ", x)                                    # bracketed codes / aliases
  x <- gsub("\\b(HCL|HYDROCHLORIDE|DIHYDROCHLORIDE|MESYLATE|DIMESYLATE|MALEATE|TOSYLATE|SODIUM|CITRATE|PHOSPHATE|SULFATE|ACETATE|FUMARATE|TARTRATE|BESYLATE|LACTATE|SUCCINATE|MALATE)\\b", " ", x)
  gsub("[^A-Z0-9]", "", x)
}
depmap_cmp <- read.csv("data/cell_lines/DepMap_26Q3/PortalCompounds.csv", stringsAsFactors = FALSE)
cmp_names <- local({
  syn <- strsplit(ifelse(is.na(depmap_cmp$Synonyms), "", depmap_cmp$Synonyms), ";")
  d <- data.frame(key = compound_key(c(depmap_cmp$CompoundName, unlist(syn))),
                  i = c(seq_len(nrow(depmap_cmp)), rep(seq_len(nrow(depmap_cmp)), lengths(syn))))
  d <- d[d$key != "" & !is.na(d$key), ]
  d <- d[nzchar(depmap_cmp$GeneSymbolOfTargets[d$i]) | !duplicated(d$key), ]   # prefer entries with targets
  d[!duplicated(d$key), ]
})
# names: one or more candidate names per drug (main name first); returns targets + mechanism
drug_targets <- function(names_list) {
  hit <- vapply(names_list, function(nm) {
    k <- compound_key(nm); k <- k[k != ""]
    i <- cmp_names$i[match(k, cmp_names$key)]; i <- i[!is.na(i)]
    if (length(i)) i[1] else NA_integer_ }, 0L)
  data.frame(depmap_compound = ifelse(is.na(hit), "", depmap_cmp$CompoundID[hit]),
             targets = ifelse(is.na(hit), "", gsub(" ", "", depmap_cmp$GeneSymbolOfTargets[hit])),
             stringsAsFactors = FALSE)
}
# split "Axitinib (AG-013736)" into its name and the bracketed alias
drug_aliases <- function(x) lapply(x, function(s) {
  al <- regmatches(s, gregexpr("\\(([^()]*)\\)", s))[[1]]
  c(trimws(sub("\\s*\\(.*$", "", s)), gsub("^\\(|\\)$", "", al), s) })

# ---- mutation columns from a sample x gene logical matrix (rows = sample_id; NA = not assayed)
add_mutations <- function(samples, M) {
  M <- M[match(samples$sample_id, rownames(M)), , drop = FALSE]
  for (g in aml_genes) samples[[paste0("mut_", g)]] <- if (g %in% colnames(M)) M[, g] else NA
  samples$mutations <- apply(M, 1, function(r) if (all(is.na(r))) NA else paste(colnames(M)[!is.na(r) & r], collapse = ";"))
  samples
}
mut_matrix <- function(sample_ids, genes, long_sample, long_gene) {
  M <- matrix(FALSE, length(sample_ids), length(genes), dimnames = list(sample_ids, genes))
  ok <- long_sample %in% sample_ids & long_gene %in% genes
  M[cbind(long_sample[ok], long_gene[ok])] <- TRUE
  M
}

FIELDS <- c("sample_id", "study_id", "participant_id", "sample_type", "population", "age_years", "sex",
            "disease_oncotree", "diagnosis", "tissue", "timepoint", "risk_group", "risk_classification",
            "KMT2A_rearranged", "KMT2A_partner", "treatment_response",
            "overall_survival_days", "os_event", "event_free_survival_days", "efs_event",
            paste0("mut_", aml_genes), "mutations")
GROUPINGS <- read.delim(text = paste(
  "column\tlabel\ttrue_label\tfalse_label",
  "KMT2A_rearranged\tKMT2A rearrangement\tKMT2A-r\tKMT2A wild-type",
  paste0("mut_", aml_genes, "\t", aml_genes, " mutation\t", aml_genes, "-mutant\t", aml_genes, "-wild-type", collapse = "\n"),
  "risk_group\tRisk group\t\t",
  "timepoint\tTimepoint\t\t",
  "treatment_response\tInduction response\t\t",
  "disease_oncotree\tOncoTree code\t\t",
  "tissue\tTissue\t\t",
  "sex\tSex\t\t",
  sep = "\n"), colClasses = "character", na.strings = NULL)
GROUPINGS$kind <- ifelse(GROUPINGS$true_label != "", "binary", "categorical")

# --------------------------------------------------------------- writer ------
# layers: list of list(key, X = features x samples, features = data.frame(label, gene_symbol, ...),
#                      label, modality, unit, source, lower_is_stronger)
write_cohort <- function(id, samples, layers, extra_groupings = NULL, notes = "") {
  out <- file.path(out_root, id)
  for (d in c("features", "values")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
  manifest <- list()
  add <- function(path, kind, desc) manifest[[length(manifest) + 1]] <<- data.frame(path = path, kind = kind, description = desc)
  wt <- function(x, rel, gz = FALSE) {
    con <- if (gz) gzfile(file.path(out, rel), "w") else file(file.path(out, rel), "w")
    write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); close(con)
  }
  stopifnot(!anyDuplicated(samples$sample_id))
  for (j in names(samples)) if (is.character(samples[[j]]))     # line breaks / tabs in text (e.g. Excel karyotype
    samples[[j]] <- trimws(gsub("\\s*(_x000D_|[\r\n\t])+\\s*", " ", samples[[j]]))   # cells) would split TSV rows
  samples$study_id <- id; samples$sample_type <- "patient"
  for (f in setdiff(FIELDS, names(samples))) samples[[f]] <- NA
  samples <- samples[, c(FIELDS, setdiff(names(samples), FIELDS))]

  idx <- list()
  for (L in layers) {
    stopifnot(!anyDuplicated(rownames(L$X)))                    # feature IDs must be unique in values/
    X <- L$X[, colnames(L$X) %in% samples$sample_id, drop = FALSE]
    ok <- rowSums(!is.na(X)) >= 3 & apply(X, 1, function(v) isTRUE(sd(v, na.rm = TRUE) > 0))
    X <- X[ok, , drop = FALSE]; feats <- L$features[ok, , drop = FALSE]
    X <- X[, order(colnames(X)), drop = FALSE]
    wt(data.frame(feature = rownames(X), feats, check.names = FALSE), file.path("features", paste0(L$key, ".tsv")))
    wt(data.frame(feature = rownames(X), signif(X, 5), check.names = FALSE),
       file.path("values", paste0(L$key, ".tsv.gz")), gz = TRUE)
    add(file.path("features", paste0(L$key, ".tsv")), "features", paste(L$label, "- feature annotation"))
    add(file.path("values", paste0(L$key, ".tsv.gz")), "values", paste(L$label, "- values per sample"))
    idx[[L$key]] <- data.frame(key = L$key, label = L$label, modality = L$modality, unit = L$unit,
                               source = L$source, lower_is_stronger = L$lower_is_stronger,
                               n_features = nrow(X), n_samples = sum(colSums(!is.na(X)) > 0),
                               duplicate_of = if (is.null(L$duplicate_of)) "" else L$duplicate_of)
    cat(sprintf("  %-14s %4d samples x %6d features\n", L$key, idx[[L$key]]$n_samples, nrow(X)))
  }
  samples <- samples[order(samples$sample_id), ]
  wt(samples, "samples.tsv"); add("samples.tsv", "samples", "Patient samples with harmonized annotation")
  wt(do.call(rbind, idx), "datasets.tsv"); add("datasets.tsv", "datasets", "Index of layers")

  g <- rbind(GROUPINGS, extra_groupings)
  g <- g[g$column %in% names(samples), ]
  g$n_ge3 <- vapply(g$column, function(c) sum(table(samples[[c]]) >= 3), 0L)
  g <- g[g$n_ge3 >= 2, setdiff(names(g), "n_ge3")]
  wt(g, "groupings.tsv"); add("groupings.tsv", "groupings", "Sample annotation columns offered as group comparisons")

  commit <- suppressWarnings(tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
                                      error = function(e) character(0)))
  wt(data.frame(key = c("prepared_at", "processing_commit", "r_version", "results_version", "notes"),
                value = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), if (length(commit)) commit else "",
                          R.version.string, "3", notes)), "provenance.tsv")
  add("provenance.tsv", "provenance", "How and when these files were prepared")
  wt(do.call(rbind, manifest), "manifest.tsv")
  cat(sprintf("  -> %s: %d samples, %d participants, %d groupings\n", out, nrow(samples),
              length(unique(samples$participant_id)), nrow(g)))
}

expr_layer <- function(E, label, unit, source) {
  E <- hgnc_rows(E)
  E <- E[rownames(E) %in% protein_coding, , drop = FALSE]
  list(key = "expression", X = E, features = data.frame(label = rownames(E), gene_symbol = rownames(E)),
       label = label, modality = "transcriptomics", unit = unit, source = source, lower_is_stronger = FALSE)
}
drug_layer <- function(key, X, names_list, label, unit, source, lower, mechanism = "", drug_class = "") {
  tg <- drug_targets(names_list)
  list(key = key, X = X,
       features = data.frame(label = vapply(names_list, `[`, "", 1), gene_symbol = "",
                             compound_name = rownames(X), compound_key = compound_key(vapply(names_list, `[`, "", 1)),
                             targets = tg$targets, depmap_compound = tg$depmap_compound,
                             mechanism = mechanism, drug_class = drug_class, stringsAsFactors = FALSE),
       label = label, modality = "drug_screen", unit = unit, source = source, lower_is_stronger = lower)
}
# proteomics: gene-level protein (rows = genes) and single-site phospho (rows = sites, "GENE_S123" labels)
protein_layer <- function(P, unit, source) {
  P <- P[!is.na(rownames(P)) & rownames(P) != "", , drop = FALSE]
  P <- hgnc_rows(P, keep_unmapped = TRUE)
  list(key = "protein", X = P, features = data.frame(label = rownames(P), gene_symbol = rownames(P)),
       label = "Protein abundance", modality = "proteomics", unit = unit, source = source, lower_is_stronger = FALSE)
}
phospho_layer <- function(X, site_label, window, unit, source) {
  gene <- sub("_[^_]*$", "", site_label)
  cur <- to_hgnc(gene); cur[is.na(cur)] <- gene[is.na(cur)]          # current HGNC symbol in the label
  site_label <- paste0(cur, substring(site_label, nchar(gene) + 1)); gene <- cur   # IDs (rownames) unchanged
  list(key = "phospho", X = X,
       features = data.frame(label = site_label, gene_symbol = gene, site = sub("^.*_", "", site_label),
                             sequence_window = window, stringsAsFactors = FALSE),
       label = "Phosphosite abundance", modality = "phosphoproteomics", unit = unit, source = source,
       lower_is_stronger = FALSE)
}

oncotree_of <- function(x) {
  x <- toupper(x)
  ifelse(grepl("PROMYELOCYTIC|^APL$|PML-RARA", x), "APLPMLRARA",
  ifelse(grepl("AMBIGUOUS|MPAL|MIXED PHENOTYPE", x), "MPAL",
  ifelse(grepl("MYELODYSPLASTIC/MYELOPROLIFERATIVE|MDS/MPN", x), "MDS/MPN",
  ifelse(grepl("MYELOPROLIFERATIVE|^MPN", x), "MPN",
  ifelse(grepl("MYELODYSPLASTIC|^MDS", x), "MDS",
  ifelse(grepl("AML|MYELOID LEUK", x), "AML", NA))))))
}
risk_of <- function(x) { x <- tolower(x)
  ifelse(grepl("^fav|^1$", x), "favorable", ifelse(grepl("^int|^2$", x), "intermediate",
         ifelse(grepl("^adv|^3$", x), "adverse", NA))) }

# ================================================================ BEAT AML ==
cat("BEAT_AML\n")
B <- file.path(dd, "BeatAML")
cl <- xlsx_read(file.path(B, "beataml_wv1to4_clinical.xlsx"), "summary")
cl$sample_id <- ifelse(is.na(cl$dbgap_rnaseq_sample), cl$dbgap_dnaseq_sample, cl$dbgap_rnaseq_sample)
cl <- cl[!is.na(cl$sample_id) & !duplicated(cl$sample_id), ]
fus <- cl$consensusAMLFusions
s <- data.frame(
  sample_id = cl$sample_id, participant_id = cl$dbgap_subject_id,
  population = ifelse(num(cl$ageAtDiagnosis) < 18, "pediatric", "adult"),
  age_years = num(cl$ageAtDiagnosis), sex = tolower(cl$consensus_sex),
  disease_oncotree = ifelse(grepl("PML-RARA", cl$specificDxAtAcquisition), "APLPMLRARA", oncotree_of(cl$dxAtSpecimenAcquisition)),
  diagnosis = cl$specificDxAtAcquisition,
  tissue = c("Bone Marrow Aspirate" = "bone_marrow", "Peripheral Blood" = "peripheral_blood",
             "Leukapheresis" = "peripheral_blood")[cl$specimenType],
  timepoint = c("Initial Diagnosis" = "diagnosis", "Relapse" = "relapse", "Remission" = "remission",
                "Residual" = "residual")[cl$diseaseStageAtSpecimenCollection],
  risk_group = risk_of(cl$ELN2017), risk_classification = ifelse(is.na(risk_of(cl$ELN2017)), NA, "ELN 2017"),
  KMT2A_rearranged = ifelse(fus %in% "N/A", NA, grepl("KMT2A", fus)),
  KMT2A_partner = ifelse(grepl("KMT2A", fus), sub("-?KMT2A-?", "", fus), NA),
  treatment_response = ifelse(grepl("^Complete Response", cl$responseToInductionTx), "CR",
                       ifelse(grepl("Refractory", cl$responseToInductionTx), "refractory", NA)),
  overall_survival_days = num(cl$overallSurvival),
  os_event = c(Dead = TRUE, Alive = FALSE)[cl$vitalStatus],
  stringsAsFactors = FALSE)
s$sex[!s$sex %in% c("female", "male")] <- "unknown"
s$KMT2A_partner[s$KMT2A_partner %in% c("re", "_re", "")] <- NA                 # "KMT2A_re": partner not given
# mutations: WES calls (dna sample) + clinical FLT3-ITD / NPM1 tests
mu <- read.delim(file.path(B, "beataml_wes_wv1to4_mutations_dbgap.txt"), colClasses = "character")[, c("dbgap_sample_id", "symbol")]
dna2s <- setNames(cl$sample_id, cl$dbgap_dnaseq_sample)
M <- mut_matrix(cl$sample_id, aml_genes, dna2s[mu$dbgap_sample_id], mu$symbol)
M[cl$`FLT3-ITD` %in% "positive", "FLT3"] <- TRUE
M[cl$NPM1 %in% "positive", "NPM1"] <- TRUE
has_wes <- !is.na(cl$dbgap_dnaseq_sample) & cl$dbgap_dnaseq_sample %in% mu$dbgap_sample_id |
           (!is.na(cl$dbgap_dnaseq_sample) & cl$exomeSeq %in% c("y", "yes", "Y"))
M[!has_wes, setdiff(aml_genes, c("FLT3", "NPM1"))] <- NA
M[!has_wes & is.na(cl$`FLT3-ITD`), "FLT3"] <- NA
M[!has_wes & is.na(cl$NPM1), "NPM1"] <- NA
s <- add_mutations(s, M)
# LSC-like phenotype class (Primitive / Committed), keyed by labId
map <- xlsx_read(file.path(B, "beataml_waves1to4_sample_mapping.xlsx"), "sample_map")
pvc <- read.csv(file.path(B, "beataml.pvc.csv"), stringsAsFactors = FALSE)
s$cell_state <- pvc$Classification[match(map$labId[match(s$sample_id, map$dbgap_rnaseq_sample)], pvc$ID)]
# expression (log2 normalized; RNA-seq samples)
ex <- read.delim(file.path(B, "beataml_waves1to4_norm_exp_dbgap.txt"), check.names = FALSE)
E <- collapse_rows(num_matrix(ex, 5:ncol(ex)), ensembl_symbol(ex$stable_id))   # by Ensembl ID
# drug response: probit-fit AUC, one value per specimen x inhibitor
dr <- read.delim(file.path(B, "beataml_probit_curve_fits_v4_dbgap.txt"), stringsAsFactors = FALSE)
dr$sample_id <- ifelse(is.na(dr$dbgap_rnaseq_sample) | dr$dbgap_rnaseq_sample == "", dr$dbgap_dnaseq_sample, dr$dbgap_rnaseq_sample)
D <- tapply(dr$auc, list(dr$inhibitor, dr$sample_id), mean)
write_cohort("BEAT_AML", s, list(
  expr_layer(E, "RNA expression", "log2 normalized expression", "Beat AML waves 1-4"),
  drug_layer("drug_auc", D, drug_aliases(rownames(D)), "Drug response (AUC)", "AUC (probit fit)",
             "Beat AML ex vivo inhibitor assay", TRUE)),
  extra_groupings = data.frame(column = "cell_state", label = "Cell state (Primitive / Committed)",
                               true_label = "", false_label = "", kind = "categorical"),
  notes = "sample_id = dbGaP RNA-seq sample (DNA sample when no RNA-seq)")
beat_s <- s; beat_D <- D; beat_map <- map

# ======================================================== PTRC (Beat AML proteomics) ==
# Proteomics of Beat AML specimens (OHSU / PNNL PTRC). sample_id = Beat AML labId. Specimens in Beat AML
# waves 1-4 take their annotation, survival and drug AUC from the matching Beat AML sample; RNA is not
# repeated here (it is in BEAT_AML), and the drug layer is flagged as a copy of BEAT_AML's.
cat("PTRC_BEAT_AML\n")
PT <- file.path(dd, "PTRC")
pt <- read.delim(file.path(PT, "PTRC_metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
m  <- beat_map[match(pt$Barcode.ID, beat_map$labId), ]
bsid <- ifelse(is.na(m$dbgap_rnaseq_sample), m$dbgap_dnaseq_sample, m$dbgap_rnaseq_sample)
bi <- match(bsid, beat_s$sample_id)
s <- beat_s[bi, setdiff(names(beat_s), c("study_id", "sample_type"))]
s$sample_id <- pt$Barcode.ID
s$participant_id <- ifelse(is.na(bi), paste0("PTRC_", pt$Barcode.ID), beat_s$participant_id[bi])
s$beat_aml_sample_id <- ifelse(is.na(bi), NA, bsid)
un <- is.na(bi)                                                   # not in Beat AML waves 1-4
s$disease_oncotree[un] <- "AML"
s$tissue[un] <- c("Bone Marrow Aspirate" = "bone_marrow", "Peripheral Blood" = "peripheral_blood",
                  "Leukapheresis" = "peripheral_blood")[pt$specimen.type[un]]
s$timepoint[un] <- ifelse(pt$InitialAMLDiagnosis[un] %in% TRUE, "diagnosis", NA)
fus <- pt$consensusAMLFusions
s$KMT2A_rearranged[un] <- ifelse(is.na(fus[un]), NA, grepl("KMT2A", fus[un]))
for (g in c("NPM1", "DNMT3A", "CEBPA", "IDH2", "IDH1", "TET2", "NRAS")) {
  v <- num(pt[[paste0(g, ".mutation.status")]]) > 0
  s[[paste0("mut_", g)]][un] <- v[un]
}
s$mut_FLT3[un] <- pt$FLT3.ITD[un] %in% c(TRUE, "TRUE")
s$mutations[un] <- apply(s[un, paste0("mut_", aml_genes)], 1, function(r)
  if (all(is.na(r))) NA else paste(aml_genes[!is.na(r) & r %in% TRUE], collapse = ";"))
s$cell_state <- ifelse(is.na(pt$pvc_classification), s$cell_state, pt$pvc_classification)
full2id <- setNames(pt$Barcode.ID, pt$SampleID.full)
pr <- read.delim(file.path(PT, "PTRC_protein_gene_ratio_MD.tsv"), check.names = FALSE)
P  <- collapse_rows(num_matrix(pr, 2:ncol(pr)), pr$GeneSymbol); colnames(P) <- full2id[colnames(P)]
ph <- read.delim(file.path(PT, "PTRC_phospho_single_site_ratio_MD.tsv"), check.names = FALSE)
X  <- num_matrix(ph, grep("^PTRC_", names(ph))); rownames(X) <- ph$Index; colnames(X) <- full2id[colnames(X)]
keep <- !is.na(bi) & bsid %in% colnames(beat_D)
Dp <- beat_D[, bsid[keep], drop = FALSE]; colnames(Dp) <- pt$Barcode.ID[keep]
dl <- drug_layer("drug_auc", Dp, drug_aliases(rownames(Dp)), "Drug response (AUC)", "AUC (probit fit)",
                 "Beat AML ex vivo inhibitor assay (same values as BEAT_AML)", TRUE)
dl$duplicate_of <- "BEAT_AML/drug_auc"
write_cohort("PTRC_BEAT_AML", s, list(
  protein_layer(P, "log2 ratio (TMT, median-centered)", "PTRC global proteomics"),
  phospho_layer(X, ph$short_name, ph$SequenceWindow, "log2 intensity (TMT)", "PTRC phosphoproteomics"),
  dl),
  extra_groupings = data.frame(column = "cell_state", label = "Cell state (Primitive / Committed)",
                               true_label = "", false_label = "", kind = "categorical"),
  notes = paste0("Beat AML specimens with proteomics; sample_id = Beat AML labId. ", sum(!un), " of ", nrow(pt),
                 " specimens are in Beat AML waves 1-4 and take their clinical data and drug AUC from there."))

# ======================================================== CPTAC adult AML (data freeze v8) ==
cat("CPTAC_AML_ADULT\n")
CF <- file.path(dd, "adult_AML_data_freeze")
cm <- read.delim(file.path(CF, "proteogenomics_metadata_20250112.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
st <- function(col) if (col %in% names(cm)) num(cm[[col]]) > 0 else rep(NA, nrow(cm))
kmt <- st("KMT2A-Fusions.mutation.status") | st("KMT2A-MLLT3.mutation.status") | cm$WHO_classification %in% "KMT2A-rearrangement"
s <- data.frame(
  sample_id = cm$patientID, participant_id = sub("^([^-]+-[^-]+).*$", "\\1", cm$Case_ID),
  population = ifelse(num(cm$Consent.age) < 18, "pediatric", "adult"), age_years = num(cm$Consent.age),
  sex = tolower(cm$Consent.sex), disease_oncotree = "AML", diagnosis = cm$WHO_classification,
  tissue = c(BMA = "bone_marrow", BM = "bone_marrow", PB = "peripheral_blood")[cm$protein_source],
  risk_group = risk_of(cm$ELN_classification), risk_classification = "ELN",
  KMT2A_rearranged = kmt, KMT2A_partner = ifelse(st("KMT2A-MLLT3.mutation.status") %in% TRUE, "MLLT3", NA),
  cell_state = cm$Primitive_Committed_Classification, stringsAsFactors = FALSE)
M <- sapply(aml_genes, function(g) st(paste0(g, ".mutation.status"))); rownames(M) <- cm$patientID
M[, "FLT3"] <- M[, "FLT3"] | st("FLT3_ITD.mutation.status") | st("FLT3_TKD.mutation.status")
M <- M[, colSums(!is.na(M)) > 0, drop = FALSE]                  # genes not reported stay NA (e.g. U2AF1)
s <- add_mutations(s, M)
pg <- read.delim(file.path(CF, "v8/prot_ratio_gene_MD_tidy_mapFile_proteogenomics_RBC_corrected.tsv"), check.names = FALSE)
P  <- collapse_rows(num_matrix(pg, 2:ncol(pg)), pg$Gene_name)
ps <- read.delim(file.path(CF, "v8/phospho_ratio_signle-site_MD_tidy_mapFile_proteogenomics_RBC_corrected.tsv"), check.names = FALSE)
X  <- num_matrix(ps, 4:ncol(ps)); rownames(X) <- ps$Index
write_cohort("CPTAC_AML_ADULT", s, list(
  protein_layer(P, "log2 ratio (TMT, RBC-corrected)", "CPTAC adult AML data freeze v8"),
  phospho_layer(X, ps$Gene_site, ps$SequenceWindow, "log2 ratio (TMT, RBC-corrected)", "CPTAC adult AML data freeze v8")),
  extra_groupings = data.frame(column = "cell_state", label = "Cell state (Primitive / Committed)",
                               true_label = "", false_label = "", kind = "categorical"),
  notes = "Proteogenomic cohort; protein and phosphosite ratios corrected for red-blood-cell contamination. No survival or drug response in this release.")

# ================================================ Zenodo 14160480 (AML proteogenomics) ==
# Protein (374 patients, batch-adjusted), phosphopeptides (217), RNA-seq TPM (361), single-drug IC50
# (102 non-M3 patients; Supp. Data 10-11) and clinical information (Supp. Data 1).
cat("ZENODO_14160480\n")
Z <- file.path(dd, "zenodo_14160480")
first_gene <- function(g) sub(";.*$", "", g)
pz <- read.csv(file.path(Z, "Protein abundance of 374 AML patients after adjusting batch effects.csv"), check.names = FALSE)
P  <- num_matrix(pz, 3:ncol(pz))
rmin <- apply(P, 1, min, na.rm = TRUE)                           # per-protein fill value used for missing
fill <- abs(P - rmin) < 1e-9 & rowSums(abs(P - rmin) < 1e-9, na.rm = TRUE) >= 3
P[fill] <- NA
P  <- collapse_rows(P, first_gene(pz$PG.Genes))
phz <- read.csv(file.path(Z, "Phosphopeptides of 217 AML patients.csv"), check.names = FALSE)
scol <- grep("^T[0-9]+$", names(phz))
site <- paste0(first_gene(phz$PG.Genes), "_", phz$PTM.SiteAA, phz$PTM.SiteLocation)
fid  <- paste0(sub(";.*$", "", phz$PTM.ProteinId), "_", phz$PTM.SiteAA, phz$PTM.SiteLocation)
X  <- collapse_rows(log2(num_matrix(phz, scol)), fid)            # raw intensity -> log2; duplicated sites averaged
lab <- site[match(rownames(X), fid)]; win <- phz$PTM.FlankingRegion[match(rownames(X), fid)]
ok <- !is.na(lab) & !startsWith(lab, "_")
X <- X[ok, , drop = FALSE]; lab <- lab[ok]; win <- win[ok]
rz <- xlsx_big_to_df(file.path(Z, "Gene expression matrix of 361 AML patients with RNA-seq data in transcript per million (TPM) format.xlsx"), 363)
E  <- collapse_rows(log2(num_matrix(rz, grep("^T[0-9]+$", names(rz))) + 1), rz$symbol)
d11 <- xlsx_read(file.path(Z, "Supp_Data_11.xlsx"), "Single drug", skip = 1)
D  <- tapply(log10(num(d11$IC50_measured)), list(d11$Compound_name, d11$Sample_ID), mean)
d10 <- xlsx_read(file.path(Z, "Supp_Data_10.xlsx"), 1, skip = 1)
k10 <- match(compound_key(rownames(D)), compound_key(d10$Drug))
# clinical information (Dataset S1): 374 newly diagnosed patients; OS / EFS in days
cz <- xlsx_read(file.path(Z, "Supp_Data_1.xlsx"), 1, skip = 1)
cz <- cz[grepl("^T[0-9]+$", cz$Sample_ID), ]
ids <- sort(unique(c(cz$Sample_ID, colnames(P), colnames(X), colnames(E), colnames(D))))
ci <- match(ids, cz$Sample_ID); cc <- cz[ci, ]
fus <- ifelse(cc$Gene_fusions %in% c("Negative", "NA"), NA, cc$Gene_fusions)
apl <- cc$FAB %in% "M3" | grepl("PML::RARA", cc$Gene_fusions)
s <- data.frame(
  sample_id = ids, participant_id = ids,
  population = ifelse(num(cc$Age) < 18, "pediatric", "adult"), age_years = num(cc$Age),
  sex = c(F = "female", M = "male")[cc$Sex],
  disease_oncotree = ifelse(apl, "APLPMLRARA", "AML"), diagnosis = cc$WHO_classification_2022,
  timepoint = ifelse(is.na(ci), NA, "diagnosis"),
  risk_group = risk_of(cc$ELN_risk_2017), risk_classification = ifelse(is.na(risk_of(cc$ELN_risk_2017)), NA, "ELN 2017"),
  KMT2A_rearranged = ifelse(is.na(ci), NA, grepl("KMT2A", cc$Gene_fusions) | cc$WHO_classification_2022 %in% "KMT2A rearrangement"),
  KMT2A_partner = ifelse(grepl("KMT2A::", cc$Gene_fusions), sub("^.*KMT2A::", "", cc$Gene_fusions), NA),
  overall_survival_days = num(cc$OS), os_event = c("1" = TRUE, "0" = FALSE)[cc$OS_status],
  event_free_survival_days = num(cc$EFS), efs_event = c("1" = TRUE, "0" = FALSE)[cc$EFS_status],
  proteomic_subtype = ifelse(cc$Proteomic_subtype %in% "NA", NA, cc$Proteomic_subtype),
  FAB = ifelse(cc$FAB %in% c("AML", "NA"), NA, cc$FAB), who_2022 = cc$WHO_classification_2022,
  gene_fusion = ifelse(is.na(ci), NA, ifelse(is.na(fus), "none detected", fus)),
  induction_group = ifelse(cc$Induction_group %in% "NA", NA, cc$Induction_group),
  HSCT = c("1" = TRUE, "0" = FALSE)[cc$HSCT],
  stringsAsFactors = FALSE)
zg <- data.frame(column = c("proteomic_subtype", "FAB", "who_2022", "gene_fusion", "induction_group"),
                 label = c("Proteomic subtype (S1-S8)", "FAB subtype", "WHO 2022 class", "Gene fusion", "Induction group"),
                 true_label = "", false_label = "", kind = "categorical")
write_cohort("ZENODO_14160480", s, list(
  expr_layer(E, "RNA expression", "log2(TPM+1)", "Zenodo 14160480"),
  protein_layer(P, "log2 abundance (batch-adjusted)", "Zenodo 14160480"),
  phospho_layer(X, lab, win, "log2 intensity", "Zenodo 14160480"),
  drug_layer("drug_ic50", D, drug_aliases(rownames(D)), "Drug response (IC50)", "log10 IC50 (nM)",
             "Zenodo 14160480 (Supp. Data 11)", TRUE,
             mechanism = ifelse(is.na(k10), "", d10$Targets[k10]), drug_class = ifelse(is.na(k10), "", d10$Category[k10]))),
  extra_groupings = zg,
  notes = paste("Proteogenomic cohort of 374 newly diagnosed AML patients deposited at Zenodo (record 14160480);",
                "clinical data from Dataset S1 (no per-gene mutation calls). Drug response for 102 non-M3 patients.",
                "Protein fill values for missing measurements set to NA; phosphosite intensities log2-transformed.",
                "Drug-combination CI values are not included."))

# ================================================================ LEE 2018 ==
cat("LEE_2018\n")
L <- file.path(dd, "Nat_Comm_2018")
cl <- xlsx_read(file.path(L, "Supp_Data_3.xlsx"), 1, skip = 2)
cl <- cl[grepl("^AML", trimws(cl$`Sample ID`)), ]
cl$sample_id <- trimws(cl$`Sample ID`)
kary <- cl$Karyotype
s <- data.frame(
  sample_id = cl$sample_id, participant_id = cl$sample_id,
  population = ifelse(num(cl$`Dx Age`) < 18, "pediatric", "adult"), age_years = round(num(cl$`Dx Age`), 1),
  sex = c(F = "female", M = "male")[toupper(trimws(cl$Sex))],
  disease_oncotree = "AML", diagnosis = "AML",
  timepoint = ifelse(yn(cl$`Relapse/ Refractory at sample`), "relapse", "diagnosis"),
  risk_group = risk_of(cl$`European Leukemia Net (ELN) cyto risk (Döhner et al 2010)`),
  risk_classification = "ELN 2010 (cytogenetic)",
  KMT2A_rearranged = ifelse(is.na(kary), NA, grepl("11\\)\\([^)]*q23|t\\(11;[0-9]+\\)\\(q23|11q23", kary)),
  treatment_response = c("1" = "CR", "0" = "refractory")[cl$`CR with initial treatment given after sample yes=1, no=0`],
  stringsAsFactors = FALSE)
M <- cbind(NPM1 = yn(cl$NPM1), FLT3 = yn(cl$`FLT3 clinic`)); rownames(M) <- cl$sample_id
s <- add_mutations(s, M)
ex <- xlsx_read(file.path(L, "Supp_Data_1.xlsx"), "patientExpression", skip = 2)
E <- collapse_rows(num_matrix(ex, grep("^AML", names(ex))), ex$UID)
au <- xlsx_read(file.path(L, "Supp_Data_1.xlsx"), "patientAUC")
D <- num_matrix(au, grep("^AML", names(au))); rownames(D) <- au$UID
colnames(D) <- sub("^AML\\.?", "AML", colnames(D))
dc <- xlsx_read(file.path(L, "Supp_Data_2.xlsx"), 1, skip = 1)
pd <- dc$Pharmacodynamics[match(rownames(D), dc$Name)]
write_cohort("LEE_2018", s, list(
  expr_layer(E, "RNA expression", "log2 expression (array)", "Lee et al. 2018"),
  drug_layer("drug_auc", D, lapply(rownames(D), function(n) c(n, dc$`Alternate Name`[match(n, dc$Name)])),
             "Drug response (AUC)", "AUC", "Lee et al. 2018 ex vivo panel", TRUE,
             mechanism = ifelse(is.na(pd), "", pd))),
  notes = "Survival is given without an event indicator, so it is not included. KMT2A from karyotype (11q23).")

# ======================================================= WANG 2022 (pediatric) ==
cat("WANG_2022_PED\n")
W  <- file.path(dd, "Blood_Cancer_Discovery_2022/Supp_Table.xlsx")
s2 <- xlsx_read(W, "Table S2-Clinical Summary")
s3 <- xlsx_read(W, "Table S3-Sample Dashboard")
s3 <- s3[grepl("^LEU[0-9]+$", trimws(s3$`Sample code`)), ]                       # drop footnote rows
names(s3)[names(s3) == "Diagnostic/ relapsed"] <- "stage"
s3$participant <- ifelse(grepl("\\(LEU", s3$stage), sub("^.*\\((LEU[0-9]+)\\).*$", "\\1", s3$stage), s3$`Sample code`)
p  <- s2[match(s3$participant, s2$`Patient code`), ]
s <- data.frame(
  sample_id = s3$`Sample code`, participant_id = s3$participant, population = "pediatric",
  age_years = num(p$`Age at Dx (y)`), sex = c(F = "female", M = "male")[p$Gender],
  disease_oncotree = oncotree_of(s3$Diagnosis), diagnosis = s3$Diagnosis,
  tissue = c(BM = "bone_marrow", PB = "peripheral_blood")[s3$Source],
  timepoint = ifelse(grepl("^REL", s3$stage), "relapse", "diagnosis"),
  risk_group = risk_of(p$`Cytogenetic risk group`), risk_classification = "cytogenetic",
  KMT2A_rearranged = yn(p$`KMT2A-rearranged`),
  treatment_response = c(Y = "CR", N = "refractory")[p$CR],
  event_free_survival_days = round(num(p$`Duration to event/ censor (months)`) * DAYS_PER_MONTH),
  efs_event = yn(p$Event),
  stringsAsFactors = FALSE)
# mutations: targeted panel (pathogenic / likely pathogenic) + clinical FLT3 / NPM1 / CEBPA
# panel variants: the sample code is a merged cell (first row only), so fill it down. Every protein-
# altering variant at VAF >= 5% counts, as for Beat AML WES (the source's pathogenicity labels call
# e.g. NRAS G12D "likely benign", so they are not used as a filter)
v <- xlsx_read(W, "Table S10-Variant List")
v <- v[!is.na(v$`Gene Name`), ]
for (i in seq_len(nrow(v))[-1]) if (is.na(v$`Sample code`[i])) v$`Sample code`[i] <- v$`Sample code`[i - 1]
v <- v[num(v$VAF) >= 0.05 & !is.na(v$Impact), ]
panel <- xlsx_read(W, "Table S9-NGS Gene Panel", skip = 1)$Gene
M <- mut_matrix(s$sample_id, intersect(aml_genes, panel), v$`Sample code`, v$`Gene Name`)
for (g in intersect(c("FLT3", "NPM1", "CEBPA"), colnames(M))) M[yn(p[[if (g == "CEBPA") "CEBPA" else g]]) %in% TRUE, g] <- TRUE
M[!yn(s3$`Targeted resequencing`) %in% TRUE, ] <- NA
s <- add_mutations(s, M)
rn <- xlsx_read(W, "Table S13-RNA-seq CPM")
E <- collapse_rows(num_matrix(rn, grep("^LEU", names(rn))), rn$`Gene Symbol`)
d5 <- xlsx_read(W, "Table S5-AML Drug Profile")
D <- tapply(num(d5$AUC), list(d5$Drug, d5$`Sample code`), mean)
d4 <- xlsx_read(W, "Table S4-Drug Panel", skip = 1)
d4 <- d4[!is.na(d4$`Drug name`), ]
cls4 <- d4$`Drug class`; for (i in seq_along(cls4)[-1]) if (is.na(cls4[i])) cls4[i] <- cls4[i - 1]   # merged cells
dclass <- cls4[match(rownames(D), d4$`Drug name`)]
write_cohort("WANG_2022_PED", s, list(
  expr_layer(E, "RNA expression", "log2 CPM", "Wang et al. 2022 (GSE192638)"),
  drug_layer("drug_auc", D, as.list(rownames(D)), "Drug response (AUC)", "AUC (viability, 0.1-10000 nM)",
             "Wang et al. 2022 ex vivo panel", TRUE, drug_class = ifelse(is.na(dclass), "", dclass))),
  notes = "Survival is event-free survival (relapse or death). Includes AML, APL, MDS and MPAL samples.")

# ================================================= MALANI 2022 (FPM tumor board) ==
cat("MALANI_2022\n")
F <- file.path(dd, "Functional_Precision_Medicine_Tumor_Board_AML")
an <- xlsx_read(file.path(F, "File_0_Common_sample_annotation_252S.xlsx"), 1)
cl <- xlsx_read(file.path(F, "File_1.1_Clinical_summary_186_Patients.xlsx"), 1)
an$pid <- sub("^(AML_[0-9]+)_.*$", "\\1", an$Sample_ID)
p <- cl[match(an$pid, cl$Patient_ID), ]
mll <- c(MLL_AF9 = "MLLT3", MLL_MLLT1 = "MLLT1", MLL_ELL = "ELL")
mll_pos <- sapply(names(mll), function(k) yn(p[[k]]) %in% TRUE)
s <- data.frame(
  sample_id = an$Sample_ID, participant_id = an$pid,
  population = ifelse(num(p$`Age at diagnosis`) < 18, "pediatric", "adult"), age_years = num(p$`Age at diagnosis`),
  sex = tolower(p$Gender), disease_oncotree = "AML", diagnosis = p$hema_diagnoses,
  tissue = c(Blood = "peripheral_blood", "Bone marrow" = "bone_marrow")[an$Tissue],
  timepoint = c(Diagnosis = "diagnosis", Relapse = "relapse", Refractory = "refractory")[trimws(an$`Disease status`)],
  risk_group = risk_of(p$eln2017_risk_class_qp), risk_classification = "ELN 2017",
  KMT2A_rearranged = ifelse(is.na(p$Patient_ID), NA, rowSums(mll_pos) > 0 | grepl("KMT2A|MLL", p$Fusion_gene)),
  KMT2A_partner = apply(mll_pos, 1, function(r) if (any(r)) paste(mll[r], collapse = ";") else NA),
  overall_survival_days = num(p$os_time), os_event = yn(p$os_event),
  event_free_survival_days = num(p$efs_time_pt), efs_event = yn(p$efs_event_pt),
  stringsAsFactors = FALSE)
s$sex[!s$sex %in% c("female", "male")] <- NA
bm <- xlsx_read(file.path(F, "File_6_Binary_mutation_225S_57G.xlsx"), 1)
M <- t(num_matrix(bm, 2:ncol(bm)) > 0); colnames(M) <- bm[[1]]
M <- M[, intersect(aml_genes, colnames(M)), drop = FALSE]
for (k in c("FLT3_ITD", "FLT3_TKD")) M[rownames(M) %in% an$Sample_ID[yn(p[[k]]) %in% TRUE], "FLT3"] <- TRUE
s <- add_mutations(s, M)
cpm <- read.csv(file.path(F, "File_7_RNA_seq_CPM_163S_4Healthy.csv"), check.names = FALSE, row.names = 1)
cpm <- as.matrix(cpm[, grepl("^AML_", colnames(cpm))])          # already log2 CPM (has negative values)
E <- collapse_rows(cpm, ensembl_symbol(rownames(cpm)))
dss <- xlsx_read(file.path(F, "File_3.2_Drug_response_DSS_sDSS_164S_17Healthy.xlsx"), 1, skip = 2)
dss <- dss[grepl("^AML_", dss$Sample_ID), ]
lib <- xlsx_read(file.path(F, "File_2_Drug_library_515D.xlsx"), 1)
li  <- match(sort(unique(dss$Compound_ID)), lib$`FIMM batch ID`)
nm  <- ifelse(is.na(li), sort(unique(dss$Compound_ID)), lib$`Preferred name`[li])
aliases <- lapply(seq_along(nm), function(i) c(nm[i], lib$`Non-proprietary name`[li[i]],
                                               trimws(strsplit(ifelse(is.na(lib$Alias[li[i]]), "", lib$Alias[li[i]]), ",")[[1]])))
aliases <- lapply(aliases, function(a) a[!is.na(a) & a != ""])
mech <- ifelse(is.na(li), "", lib$`Mechanism/Targets`[li]); dcl <- ifelse(is.na(li), "", lib$Drug_Class[li])
layers <- list(expr_layer(E, "RNA expression", "log2 CPM", "Malani et al. 2022"))
for (m in c("DSS", "sDSS")) {
  D <- tapply(num(dss[[m]]), list(dss$Compound_ID, dss$Sample_ID), mean)[sort(unique(dss$Compound_ID)), , drop = FALSE]
  layers[[m]] <- drug_layer(paste0("drug_", m), D, aliases,
                            paste0("Drug response (", m, ")"), if (m == "DSS") "DSS" else "selective DSS (vs healthy)",
                            "Malani et al. 2022 DSRT", FALSE, mechanism = mech, drug_class = dcl)
}
write_cohort("MALANI_2022", s, layers,
             notes = "DSS / sDSS: higher = more sensitive. Drug features are FIMM IDs; label = preferred name.")

# ===================================================== St. Jude pediatric AML (Umeda 2024) ==
# cBioPortal datahub study aml_stjude_2024 (Umeda et al. 2024 Nat Genet): 887 pediatric AML samples,
# RNA-seq (log2 CPM), WGS / WES mutations in an 87-gene panel, fusions / SVs; no survival or drug data
cat("AML_STJUDE_2024\n")
J <- file.path(dd, "aml_stjude_2024")
cbio <- function(f) read.delim(file.path(J, f), comment.char = "#", colClasses = "character",
                               na.strings = c("NA", ""), check.names = FALSE, quote = "")
cs <- cbio("data_clinical_sample.txt"); cp <- cbio("data_clinical_patient.txt")
p  <- cp[match(cs$PATIENT_ID, cp$PATIENT_ID), ]
sv <- cbio("data_sv.txt")
kfus <- sv[sv$Event_Info %in% "Fusion" & xor(sv$Site1_Hugo_Symbol %in% "KMT2A", sv$Site2_Hugo_Symbol %in% "KMT2A"), ]
kpart <- tapply(ifelse(kfus$Site1_Hugo_Symbol == "KMT2A", kfus$Site2_Hugo_Symbol, kfus$Site1_Hugo_Symbol),
                kfus$Sample_Id, function(x) paste(sort(unique(x)), collapse = ";"))
s <- data.frame(
  sample_id = cs$SAMPLE_ID, participant_id = cs$PATIENT_ID, population = "pediatric",
  age_years = num(p$DIAGNOSIS_AGE), sex = tolower(p$SEX),
  disease_oncotree = "AML", diagnosis = cs$WHO5_CATEGORY,
  tissue = c("Bone Marrow" = "bone_marrow", "Peripheral Blood" = "peripheral_blood")[cs$SAMPLE_SOURCE],
  timepoint = c(Diagnosis = "diagnosis", Relapse = "relapse", "Second Relapse" = "relapse")[cs$TIMEPOINT],
  KMT2A_rearranged = cs$MOLECULAR_CATEGORY %in% "KMT2Ar" | cs$SAMPLE_ID %in% names(kpart),
  KMT2A_partner = unname(kpart[cs$SAMPLE_ID]),
  molecular_category = cs$MOLECULAR_CATEGORY,
  FAB = cs$FAB_MORPHOLOGY,
  stringsAsFactors = FALSE)
s$sex[!s$sex %in% c("female", "male")] <- NA
# mutations: somatic calls (germline variants left out) in the samples with WGS / WES; FLT3-ITD are
# reported as FLT3 internal duplications in the SV file and are counted as FLT3 mutations
gp <- cbio("data_gene_panel_matrix.txt")
mf <- read.delim(file.path(J, "data_mutations.txt"), comment.char = "#", colClasses = "character", quote = "")
mf <- mf[mf$Mutation_Status != "GERMLINE", ]
itd <- sv$Sample_Id[sv$Site1_Hugo_Symbol %in% "FLT3" & sv$Site2_Hugo_Symbol %in% "FLT3" & sv$Event_Info %in% "Duplication"]
M <- mut_matrix(s$sample_id, aml_genes, c(mf$Tumor_Sample_Barcode, itd), c(mf$Hugo_Symbol, rep("FLT3", length(itd))))
M[!s$sample_id %in% gp$SAMPLE_ID[!is.na(gp$mutations)], ] <- NA
s <- add_mutations(s, M)
rn <- read.delim(file.path(J, "data_mrna_seq_cpm.txt"), check.names = FALSE, quote = "")
E <- collapse_rows(num_matrix(rn, setdiff(names(rn), "Hugo_Symbol")), rn$Hugo_Symbol)
stopifnot(all(colnames(E) %in% s$sample_id))
write_cohort("AML_STJUDE_2024", s, list(expr_layer(E, "RNA expression", "log2 CPM", "Umeda et al. 2024 (cBioPortal)")),
  extra_groupings = data.frame(column = c("molecular_category", "FAB"),
                               label = c("Molecular category", "FAB subtype"),
                               true_label = "", false_label = "", kind = "categorical"),
  notes = paste("Pediatric AML (St. Jude, Umeda et al. 2024 Nat Genet), from the cBioPortal datahub study aml_stjude_2024.",
                "Mutations from WGS / WES for 665 of 887 samples (somatic calls; FLT3-ITD from the structural variant",
                "calls). KMT2A from the molecular category or a KMT2A fusion. No risk group, survival or drug response.",
                "Includes TARGET and other published cohorts, so some patients may also appear in Wang 2022."))

# ================================================ AMLmap: TCGA, TARGET, Leucegene, LUMC ==
# AMLmap compendium (data/patients/AMLmap): RNA-seq counts corrected for study / library batch across five
# studies, clinical annotation, mutations and fusions. The Beat AML part is not used (Beat AML is processed
# above from its own release). Each of the other four studies becomes its own cohort.
cat("AMLmap\n")
A <- file.path(dd, "AMLmap")
am <- read.csv(file.path(A, "meta.csv"), check.names = FALSE, na.strings = c("NA", ""), colClasses = "character")
amlmap_studies <- c(TCGA_AML = "TCGA_LAML", TARGET_AML = "TARGET_AML", LEUCEGENE_AML = "LEUCEGENE", "100LUMC_AML" = "LUMC_AML")
am <- am[am$study %in% names(amlmap_studies), ]
# RNA: batch-corrected counts -> log2(CPM + 1); Ensembl IDs (versioned) -> current HGNC symbols
cnt <- read.csv(file.path(A, "counts/corrected_counts.csv"), check.names = FALSE, row.names = 1)
cnt <- as.matrix(cnt[, colnames(cnt) %in% am$sample_id]); storage.mode(cnt) <- "double"
cpm <- log2(t(t(cnt) / colSums(cnt, na.rm = TRUE)) * 1e6 + 1)
E_all <- collapse_rows(cpm, ensembl_symbol(sub("\\..*$", "", rownames(cpm)))); rm(cnt, cpm)
# mutations: protein-altering variants at VAF >= 5% (VAF missing: kept), plus hotspot / ITD columns
ab <- read.csv(file.path(A, "aberrations/aberrations_oh.csv"), check.names = FALSE, row.names = 1)
mu <- read.csv(file.path(A, "aberrations/mutations.csv"), colClasses = "character")
not_altering <- c("Silent", "RNA", "Intron", "3'UTR", "5'UTR", "Splice_Region", "3'Flank", "5'Flank", "IGR", "Splicing&Intron")
mu <- mu[!sub("\\(.*$", "", mu$variant) %in% not_altering & (is.na(num(mu$VAF)) | num(mu$VAF) >= 0.05), ]
hot <- list(IDH1 = "IDH1_R132", IDH2 = c("IDH2_R140", "IDH2_R172"), DNMT3A = "DNMT3A_R882",
            FLT3 = c("FLT3_ITD", "FLT3_TKD"), NPM1 = "NPM1")
sc <- read.csv(file.path(A, "scores.csv"), check.names = FALSE)
names(sc)[-1] <- gsub("[^a-z0-9]+", "_", tolower(names(sc)[-1]))        # e.g. hsc_like_score, hla_ii_score
for (st0 in names(amlmap_studies)) {
  id <- amlmap_studies[[st0]]; cat(id, "\n")
  x <- am[am$study == st0, ]
  kdx <- x$KMT2A_diagnosis
  kany <- num(ab[x$sample_id, "any_KMT2A"])
  s <- data.frame(
    sample_id = x$sample_id,                         # Leucegene: "X02H066" -> participant "02H066", the ID used
    participant_id = if (id == "LEUCEGENE") sub("^X", "", x$sample_id) else x$sample_id,   # by the Leucegene surfaceome
    population = if (id == "TARGET_AML") "pediatric" else ifelse(num(x$age) < 18, "pediatric", "adult"),
    age_years = round(num(x$age), 1), sex = tolower(x$sex), disease_oncotree = ifelse(grepl("PML.?RARA", x$WHO_2022) | x$FAB %in% "M3", "APLPMLRARA", "AML"),
    diagnosis = x$WHO_2022,
    tissue = c("bone marrow" = "bone_marrow", "peripheral blood" = "peripheral_blood")[x$tissue],
    timepoint = ifelse(x$prim_rec %in% "Primary", "diagnosis", NA),
    KMT2A_rearranged = ifelse(kdx %in% "no_KMT2A", FALSE, ifelse(grepl("^KMT2A_|^Other KMT2A$|^multi$", kdx), TRUE,
                              ifelse(is.na(kany), NA, kany > 0))),
    KMT2A_partner = ifelse(grepl("^KMT2A_", kdx), to_hgnc(sub("^KMT2A_", "", kdx)), NA),
    overall_survival_days = round(num(x$survival_years) * 365.25),
    os_event = c(Dead = TRUE, Alive = FALSE)[x$event],
    amlmap_cluster = x$clusters, who_2022 = x$WHO_2022, icc_2022 = x$ICC_2022, FAB = x$FAB,
    stringsAsFactors = FALSE)
  s$sex[!s$sex %in% c("female", "male")] <- NA
  s <- cbind(s, sc[match(s$sample_id, sc$sample_id), -1, drop = FALSE])
  m1 <- mu[mu$sample_id %in% s$sample_id, ]
  M <- mut_matrix(s$sample_id, aml_genes, m1$sample_id, m1$gene_id)
  for (g in names(hot)) M[, g] <- M[, g] | rowSums(as.matrix(ab[s$sample_id, hot[[g]], drop = FALSE]) > 0, na.rm = TRUE) > 0
  M[, colSums(M) == 0] <- NA                      # no call at all in this study: gene treated as not assayed
  s <- add_mutations(s, M)
  E <- E_all[, colnames(E_all) %in% s$sample_id, drop = FALSE]
  write_cohort(id, s, list(expr_layer(E, "RNA expression", "log2(CPM+1), batch-corrected counts", "AMLmap")),
    extra_groupings = data.frame(column = c("amlmap_cluster", "who_2022", "icc_2022", "FAB"),
                                 label = c("AMLmap cluster", "WHO 2022 class", "ICC 2022 class", "FAB subtype"),
                                 true_label = "", false_label = "", kind = "categorical"),
    notes = paste0(c(TCGA_LAML = "TCGA LAML", TARGET_AML = "TARGET AML (pediatric)", LEUCEGENE = "Leucegene",
                     LUMC_AML = "LUMC 100 AML")[[id]], ", from the AMLmap compendium. RNA from counts corrected for ",
                   "study and library batch, as log2(CPM+1). Mutations: protein-altering calls at VAF >= 5% plus ",
                   "FLT3-ITD / TKD, NPM1, IDH1/2 and DNMT3A R882 hotspot calls; genes with no call in the study are ",
                   "treated as not assayed. KMT2A from the AMLmap KMT2A diagnosis.",
                   if (id %in% c("TCGA_LAML", "TARGET_AML")) " Overall survival in days." else " No survival data.",
                   if (id == "TARGET_AML") " Some patients may also be in the St. Jude cohort (IDs cannot be matched)." else ""))
}
rm(E_all)

# ============================================ Jayavelu et al. 2022 Cancer Cell (proteogenomic subtypes) ==
# data/patients/AML_cancer_cell_2022: Data S1 (discovery cohort, 177 patients, label-free DIA proteome) and
# Data S3 (validation cohort, 75 patients, SILAC proteome); clinical data, targeted myeloid panels, OS.
# Each cohort is its own dataset. Discovery: the unimputed proteome (0 = not quantified -> NA). Validation:
# only an imputed proteome is provided (log2 SILAC ratio sample / standard).
JC <- file.path(dd, "AML_cancer_cell_2022")
jay_kmt2a <- function(k) ifelse(is.na(k), NA,             # balanced 11q23 translocations / deletions only
  grepl("t\\([0-9XY]+;11\\)\\([^;)]+;q23|t\\(11;[0-9XY]+\\)\\(q23|del\\(11\\)\\(q23|inv\\(11\\)\\([^)]*q23", k))
jay_samples <- function(cl, eln_col, eln_label) {
  num_or_na <- function(x) num(x)
  hsct <- cl$`Allogeneic HSCT`
  s <- data.frame(
    sample_id = cl$ID, participant_id = cl$ID,
    population = ifelse(num(cl$`Age at Diagnosis [Years]`) < 18, "pediatric", "adult"),
    age_years = round(num(cl$`Age at Diagnosis [Years]`), 1),
    sex = c(f = "female", m = "male")[trimws(tolower(cl[[grep("^Sex", names(cl))[1]]]))],
    disease_oncotree = "AML", diagnosis = c("de novo" = "de novo AML", sAML = "secondary AML",
                                            tAML = "therapy-related AML")[cl$Precursor],
    timepoint = "diagnosis",
    risk_group = risk_of(cl[[eln_col]]), risk_classification = ifelse(is.na(risk_of(cl[[eln_col]])), NA, eln_label),
    KMT2A_rearranged = jay_kmt2a(cl$Karyotype),
    treatment_response = c(CR1 = "CR", "no CR1" = "refractory")[cl$CR1],
    overall_survival_days = round(num(cl$`OS [months]`) * DAYS_PER_MONTH),
    os_event = c("1" = TRUE, "0" = FALSE)[cl$`Death Event`],
    relapse_free_survival_days = round(num(cl$`RFS [months]`) * DAYS_PER_MONTH),
    rfs_event = c("1" = TRUE, "0" = FALSE)[cl$`Relapse Event`],
    aml_type = cl$Precursor, FAB = cl$FAB,
    blasts_bm_pct = num_or_na(cl$`Bone marrow Blasts [%]`), blasts_pb_pct = num_or_na(cl$`Peripheral Blasts [%]`),
    karyotype = cl$Karyotype,
    HSCT = ifelse(is.na(hsct), NA, !grepl("^(no|non)\\b", tolower(hsct))),   # "TraNo ...t" = a find/replace-garbled "Transplant"
    stringsAsFactors = FALSE)
  s
}
jay_mutations <- function(s, cl, calls) {            # calls: sample, gene (already filtered); + clinical NPM1 / FLT3
  M <- mut_matrix(s$sample_id, aml_genes, calls$sample, calls$gene)
  M[cl$NPM1 %in% "Mut", "NPM1"] <- TRUE
  M[cl$FLT3 %in% c("ITD", "TKD"), "FLT3"] <- TRUE
  M[, colSums(M) == 0] <- NA                         # no call in the cohort: gene treated as not assayed
  s$FLT3_ITD <- ifelse(is.na(cl$FLT3), NA, cl$FLT3 == "ITD")
  add_mutations(s, M)
}
jay_groupings <- data.frame(column = c("aml_type", "FAB", "FLT3_ITD"), label = c("AML type (de novo / secondary)",
                            "FAB subtype", "FLT3-ITD"), true_label = c("", "", "FLT3-ITD"),
                            false_label = c("", "", "no FLT3-ITD"), kind = c("categorical", "categorical", "binary"))

cat("JAYAVELU_2022_DISCOVERY\n")
f1 <- file.path(JC, "Data S1 Discovery Cohort.xlsx")
cl <- xlsx_read(f1, "Clinical Characteristics")
s  <- jay_samples(cl, "ELN2017", "ELN 2017")
p1 <- xlsx_read(f1, "PanelSeq I"); p2 <- xlsx_read(f1, "PanelSeq II")
calls <- rbind(data.frame(sample = p1$Pat_ID, gene = p1$Gene, vaf = num(p1$Frequency)),      # pathological variants
               data.frame(sample = p2$Pat_ID, gene = p2$Gene, vaf = num(p2$`mutation load [%]`))[p2$`result variant` %in% "mutated", ])
calls <- calls[is.na(calls$vaf) | calls$vaf >= 5, ]
s <- jay_mutations(s, cl, calls)
up <- xlsx_read(f1, "Unimputed Proteome", header = FALSE)
ids <- up[1, ]; sc <- which(ids %in% s$sample_id)
P  <- matrix(num(up[-(1:2), sc]), ncol = length(sc)); colnames(P) <- ids[sc]
P[!is.na(P) & P <= 0] <- NA                                    # 0 = not quantified
P  <- collapse_rows(P, sub(";.*$", "", up[-(1:2), match("PG.Genes", up[1, ])]))   # row 1: IDs and PG.* headers
write_cohort("JAYAVELU_2022_DISCOVERY", s, list(
  protein_layer(P, "log2 intensity (label-free DIA, not imputed)", "Jayavelu et al. 2022 (Data S1)")),
  extra_groupings = jay_groupings,
  notes = paste("Discovery cohort of Jayavelu et al. 2022 Cancer Cell (177 adult AML patients at diagnosis): unimputed",
                "label-free proteome. Mutations from two targeted myeloid panels (pathological variants; panel II calls",
                "marked mutated) at VAF >= 5%, plus clinical NPM1 and FLT3-ITD / TKD. KMT2A-r from 11q23",
                "translocations in the karyotype. OS in days; relapse-free survival kept as relapse_free_survival_days."))

cat("JAYAVELU_2022_VALIDATION\n")
f3 <- file.path(JC, "Data S3 Validation Cohort.xlsx")
cl <- xlsx_read(f3, "Clinical Characteristics")
cl <- cl[!is.na(cl$ID), !grepl("^V[0-9]+$", names(cl))]
s  <- jay_samples(cl, "ELN", "ELN")
pv <- xlsx_read(f3, "PanelSeq")
calls <- data.frame(sample = pv$Pat_ID, gene = pv$Gene, vaf = num(sub("%.*$", "", pv$Coverage)))
calls <- calls[is.na(calls$vaf) | calls$vaf >= 5, ]
s <- jay_mutations(s, cl, calls)
pr <- xlsx_read(f3, "Imputed Proteome", header = FALSE)      # row 1: annotation headers; row 2: sample IDs
sc <- which(pr[2, ] %in% s$sample_id)
P  <- matrix(num(pr[-(1:2), sc]), ncol = length(sc)); colnames(P) <- pr[2, sc]
P  <- collapse_rows(P, sub(";.*$", "", pr[-(1:2), match("Gene names", pr[1, ])]))
write_cohort("JAYAVELU_2022_VALIDATION", s, list(
  protein_layer(P, "log2 ratio sample / SILAC standard (imputed)", "Jayavelu et al. 2022 (Data S3)")),
  extra_groupings = jay_groupings,
  notes = paste("Validation cohort of Jayavelu et al. 2022 Cancer Cell (75 AML patients at diagnosis): SILAC proteome,",
                "imputed by the authors (no unimputed version provided). Mutations from a targeted myeloid panel",
                "(pathological variants, VAF >= 5%) plus clinical NPM1 / FLT3; genes with no call are treated as not",
                "assayed. Risk group as reported (ELN). OS in days."))

# ======================================== Leucegene AML surfaceome (Cell Reports 2024) ==
# data/patients/AML_surfaceome_Cell_Reports_2024: Table S1 (clinical / molecular, 100 primary AML specimens) and
# Table S2 (surface proteome per specimen: MaxQuant, vsn-normalised, knn-imputed log2 intensity; SPAT score =
# the authors' surface-likelihood score, >= 8 = likely surface). 85 specimens are also in the AMLmap Leucegene
# RNA cohort; they share participant_id, so a patient is counted once.
cat("LEUCEGENE_SURFACEOME_2024\n")
LS <- file.path(dd, "AML_surfaceome_Cell_Reports_2024")
cl <- xlsx_read(file.path(LS, "Table_S1.xlsx"), "AML surfaceome cohort", skip = 2)
cl <- cl[!is.na(cl$`Specimen ID`), !grepl("^V[0-9]+$", names(cl))]
fab <- sub("^AML-", "", cl$`FAB classification`); fab[grepl("^Not classifiable", fab)] <- NA
fab <- sub("^M5([AB])$", "M5\\L\\1", fab, perl = TRUE)                       # M5A -> M5a, as other cohorts
s <- data.frame(
  sample_id = cl$`Specimen ID`, participant_id = cl$`Specimen ID`,
  population = ifelse(num(cl$`Age at diagnosis`) < 18, "pediatric", "adult"), age_years = num(cl$`Age at diagnosis`),
  sex = c(F = "female", M = "male")[cl$Sex],
  disease_oncotree = ifelse(grepl("PML::RARA", cl$`WHO 2022 classification at AML diagnosis`), "APLPMLRARA", "AML"),
  diagnosis = cl$`WHO 2022 classification at AML diagnosis`,
  timepoint = c(diagnosis = "diagnosis", refractory = "refractory", relapse = "relapse")[tolower(trimws(cl$`Status at sampling`))],
  risk_group = risk_of(trimws(cl$`Cytogenetic risk`)), risk_classification = "cytogenetic",
  KMT2A_rearranged = cl$`AML subgroup classification` %in% "KMT2A-r" | cl$`Cytogenetic group` %in% "KMT2A rearrangement",
  FLT3_ITD = c("1" = TRUE, "0" = FALSE)[cl$`FLT3-ITD`],
  aml_subgroup = cl$`AML subgroup classification`, who_2022 = cl$`WHO 2022 classification at AML diagnosis`, FAB = fab,
  cytogenetic_group = cl$`Cytogenetic group`, karyotype = cl$Karyotype,
  blasts_pct = num(cl$`Blasts (%) (blood or bone marrow smears)`), wbc = num(cl$`WBC (x 10E9/L)`),
  stringsAsFactors = FALSE)
M <- sapply(c("NPM1", "DNMT3A", "RUNX1"), function(g) c("1" = TRUE, "0" = FALSE)[cl[[g]]])   # only these were tested;
rownames(M) <- s$sample_id                                       # FLT3 = ITD only (TKD not tested): FLT3_ITD, not mut_FLT3
# TP53 was not tested in this study; take the call from the AMLmap Leucegene cohort (written above) for the same
# patient (85 of 100 specimens), so TP53-mutant vs wild-type can be compared with the cell lines. NA otherwise.
lg <- read.delim(file.path(out_root, "LEUCEGENE", "samples.tsv"), colClasses = "character", na.strings = "")
lg <- lg[!duplicated(lg$participant_id), ]
tp53 <- yn(lg$mut_TP53[match(s$participant_id, lg$participant_id)])
M <- cbind(M, TP53 = tp53)
s$TP53_source <- ifelse(is.na(tp53), NA, "Leucegene (AMLmap) mutation calls, same patient")
s <- add_mutations(s, M)
ms <- xlsx_read(file.path(LS, "Table_S2.xlsx"), "Mass spect data per specimen", header = FALSE)
ids <- ms[5, ]; sc <- which(ids %in% s$sample_id)               # row 4: headers; row 5: specimen IDs
body <- as.data.frame(ms[-(1:5), , drop = FALSE], stringsAsFactors = FALSE); body <- body[!is.na(body[, 1]), ]
X <- matrix(num(as.matrix(body[, sc])), ncol = length(sc), dimnames = list(NULL, ids[sc]))
body$spat <- num(body[, 3])
X <- collapse_rows(X, body[, 1])                     # rows = UniProt ID; a few HLA-DRB1 IDs are listed twice: averaged
body <- body[order(-ifelse(is.na(body$spat), -1, body$spat)), ]   # ... keeping the higher SPAT score
body <- body[!duplicated(body[, 1]), ]; body <- body[match(rownames(X), body[, 1]), ]
sym <- trimws(sub("[;,].*$", "", body[, 2])); cur <- to_hgnc(sym); cur[is.na(cur)] <- sym[is.na(cur)]
surf_layer <- list(key = "surface", X = X,
  features = data.frame(label = ifelse(is.na(cur) | cur == "", body[, 1], cur), gene_symbol = ifelse(is.na(cur), "", cur),
                        uniprot = body[, 1], spat_score = body$spat, likely_surface = body$spat >= 8,
                        stringsAsFactors = FALSE),
  label = "Surface protein abundance", modality = "surface_proteomics",
  unit = "log2 intensity (vsn-normalised, imputed)", source = "Leucegene AML surfaceome, Cell Reports 2024 (Table S2)",
  lower_is_stronger = FALSE)
write_cohort("LEUCEGENE_SURFACEOME_2024", s, list(surf_layer),
  extra_groupings = data.frame(column = c("aml_subgroup", "who_2022", "FAB", "FLT3_ITD"),
                               label = c("AML subgroup", "WHO 2022 class", "FAB subtype", "FLT3-ITD"),
                               true_label = c("", "", "", "FLT3-ITD"), false_label = c("", "", "", "no FLT3-ITD"),
                               kind = c("categorical", "categorical", "categorical", "binary")),
  notes = paste("Leucegene AML surfaceome (Bordeleau et al. 2024 Cell Rep 43:114260): surface proteome of 100 primary AML specimens",
                "(97 at diagnosis). Values as provided: vsn-normalised, knn-imputed log2 intensity; features keep the",
                "SPAT surface score (>= 8 = likely surface). Mutations tested: NPM1, DNMT3A, RUNX1 and FLT3-ITD only;",
                "TP53 from the AMLmap Leucegene mutation calls of the same patient (85 specimens; 7 of the 8 mutant have",
                "a complex karyotype).",
                "Risk group = cytogenetic risk. 85 specimens are also in the Leucegene RNA cohort (same participant_id).",
                "No survival data."))

# shared metadata the portal reads alongside the data
for (d in c("catalog", "schema")) { unlink(file.path("results", d), recursive = TRUE)
  file.copy(d, "results", recursive = TRUE) }
file.copy("dataset_sources.tsv", "results", overwrite = TRUE)
file.copy("docs/results_README.md", "results/README.md", overwrite = TRUE)
cat("done ->", out_root, "\n")
