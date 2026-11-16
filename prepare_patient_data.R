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
#
# Harmonization
#   genes     HGNC symbols; Ensembl IDs mapped with the Beat AML / Wang annotation; protein-coding only
#   drugs     compound_key = upper-case name without salts / brackets / punctuation, so the same drug
#             matches across cohorts; HGNC targets from DepMap PortalCompounds (name or synonym match)
#   mutations mut_<GENE> TRUE/FALSE for 20 recurrent AML genes (NA = gene not assayed for that sample)
#   MYC / mTOR RNA signature scores: mean z-score of the signature genes (data/*_signature_RNA.tsv)
#             within each cohort; *_high = above the cohort median. On DepMap AML lines this tracks the
#             published classification (r = 0.90 MYC, 0.82 mTOR; 96% / 89% same median split).
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

# ---- gene annotation: Ensembl ID -> HGNC symbol, biotype (Beat AML + Wang 2022 annotation)
gene_ann <- local({
  b <- read.delim(file.path(dd, "BeatAML/beataml_waves1to4_norm_exp_dbgap.txt"), nrows = 1e5,
                  colClasses = c(rep("character", 4), rep("NULL", 707)), check.names = FALSE)
  data.frame(ensembl = b$stable_id, symbol = b$display_label, biotype = b$biotype, stringsAsFactors = FALSE)
})
protein_coding <- unique(gene_ann$symbol[gene_ann$biotype == "protein_coding"])

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

# ---- MYC / mTOR RNA signature scores within a cohort
signatures <- list(MYC = read.delim("data/MYC_signature_RNA.tsv")$x,
                   mTOR = read.delim("data/mTOR_signature_RNA.tsv")$x)
add_signatures <- function(samples, E) {                          # E: genes x samples, log scale
  for (s in names(signatures)) {
    g <- intersect(signatures[[s]], rownames(E))
    Z <- t(scale(t(E[g, , drop = FALSE]))); Z <- Z[is.finite(rowSums(Z)), , drop = FALSE]
    sc <- colMeans(Z)
    samples[[paste0(s, "_score")]] <- round(unname(sc[samples$sample_id]), 4)
    samples[[paste0(s, "_high")]]  <- samples[[paste0(s, "_score")]] > median(sc)
  }
  samples$MYC_mTOR_category <- ifelse(is.na(samples$MYC_high), NA,
    paste0("MYC-", ifelse(samples$MYC_high, "High", "Low"), " / mTOR-", ifelse(samples$mTOR_high, "High", "Low")))
  samples
}

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
            "MYC_score", "mTOR_score", "MYC_high", "mTOR_high", "MYC_mTOR_category",
            paste0("mut_", aml_genes), "mutations")
GROUPINGS <- read.delim(text = paste(
  "column\tlabel\ttrue_label\tfalse_label",
  "MYC_high\tMYC RNA signature\tMYC-High\tMYC-Low",
  "mTOR_high\tmTOR RNA signature\tmTOR-High\tmTOR-Low",
  "MYC_mTOR_category\tMYC x mTOR RNA category\t\t",
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
  samples$study_id <- id; samples$sample_type <- "patient"
  for (f in setdiff(FIELDS, names(samples))) samples[[f]] <- NA
  samples <- samples[, c(FIELDS, setdiff(names(samples), FIELDS))]

  idx <- list()
  for (L in layers) {
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
  list(key = "protein", X = P, features = data.frame(label = rownames(P), gene_symbol = rownames(P)),
       label = "Protein abundance", modality = "proteomics", unit = unit, source = source, lower_is_stronger = FALSE)
}
phospho_layer <- function(X, site_label, window, unit, source) {
  gene <- sub("_[^_]*$", "", site_label)
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
E <- collapse_rows(num_matrix(ex, 5:ncol(ex)), ex$display_label)
s <- add_signatures(s, E[, colnames(E) %in% s$sample_id])
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
s <- add_signatures(s, E)
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
s <- add_signatures(s, E)
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
s <- add_signatures(s, E)
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
E <- collapse_rows(cpm, gene_ann$symbol[match(rownames(cpm), gene_ann$ensembl)])
s <- add_signatures(s, E)
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

# shared metadata the portal reads alongside the data
for (d in c("catalog", "schema")) { unlink(file.path("results", d), recursive = TRUE)
  file.copy(d, "results", recursive = TRUE) }
file.copy("dataset_sources.tsv", "results", overwrite = TRUE)
cat("done ->", out_root, "\n")
