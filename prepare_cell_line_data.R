# Prepare harmonized AML cell-line data for the portal.
#
#   data/cell_lines/   raw downloads (DepMap, DEMETER2, CCLE, ...; see data/README.md)
#   ->  results/cell_lines/
#         samples.tsv              one row per AML cell line, harmonized annotation (schema/)
#         groupings.tsv            annotation columns to compare (mTOR_high, MYC_high, mut_TP53, ...)
#         datasets.tsv             one row per omics layer
#         features/<key>.tsv       one row per feature: label, HGNC symbol, Entrez ID, drug annotation
#         values/<key>.tsv.gz      matrix: rows = features, columns = cell lines
#         manifest.tsv, provenance.tsv
#   and copies catalog/ and schema/ into results/.
#
# Cell lines = every DepMap model with OncotreePrimaryDisease "Acute Myeloid Leukemia".
# No statistics here: group comparisons are computed by the portal.
#
# Usage: Rscript prepare_cell_line_data.R     (from the repo root; base R + xml2)

suppressMessages(library(xml2))

dd      <- "data/cell_lines"
out_dir <- "results/cell_lines"
unlink(out_dir, recursive = TRUE)
for (d in c("features", "values")) dir.create(file.path(out_dir, d), recursive = TRUE)
manifest <- list()
add <- function(path, kind, desc) manifest[[length(manifest) + 1]] <<-
  data.frame(path = path, kind = kind, description = desc)
write_tsv <- function(x, rel, gz = FALSE) {
  f <- file.path(out_dir, rel)
  con <- if (gz) gzfile(f, "w") else file(f, "w")
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  close(con)
}

# ================================================================ cell lines ==
model <- read.csv(file.path(dd, "DepMap_26Q3/Model.csv"), stringsAsFactors = FALSE)
aml   <- model[model$OncotreePrimaryDisease == "Acute Myeloid Leukemia", ]
aml$CCLEName[aml$CCLEName == ""] <- paste0(aml$StrippedCellLineName[aml$CCLEName == ""],
                                           "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")
stopifnot(!anyDuplicated(aml$StrippedCellLineName))
line_of_model <- setNames(aml$StrippedCellLineName, aml$ModelID)
line_of_ccle  <- setNames(aml$StrippedCellLineName, aml$CCLEName)
cat("AML cell lines:", nrow(aml), "\n")

# read only the header + rows whose ID field is one of our ModelIDs (big DepMap CSVs)
read_rows <- function(path, ids) {
  hdr <- readLines(path, n = 1)
  tmp <- tempfile()
  pat <- paste0("(^|,)(", paste(ids, collapse = "|"), "),")
  system2("grep", c("-E", shQuote(pat), shQuote(path)), stdout = tmp)
  read.csv(text = c(hdr, readLines(tmp)), check.names = FALSE)
}
# read a features-in-rows table keeping only the needed sample columns
read_cols <- function(path, keep_cols, sep, skip = 0, id_cols = 1) {
  hdr <- strsplit(readLines(path, n = skip + 1)[skip + 1], sep, fixed = TRUE)[[1]]
  hdr <- gsub('^"|"$', "", hdr)
  cc  <- ifelse(hdr %in% keep_cols, "numeric", "NULL")
  cc[id_cols] <- "character"
  read.table(path, sep = sep, header = TRUE, skip = skip, colClasses = cc, quote = "\"",
             check.names = FALSE, na.strings = c("NA", "NaN", ""), strip.white = TRUE,
             comment.char = "")
}
# rows keyed by ModelID -> matrix with rownames = cell line
by_model <- function(df, id_col) {
  df <- df[df[[id_col]] %in% aml$ModelID, , drop = FALSE]
  X  <- as.matrix(df[, setdiff(seq_along(df), id_col), drop = FALSE])
  storage.mode(X) <- "double"
  rownames(X) <- line_of_model[df[[id_col]]]
  X
}

# ---------------------------------------------------------- annotation ------
# curated MYC/mTOR RNA classification (56 lines) and KMT2A status
cls <- read.delim(file.path(dd, "fig2_MYC_mTOR_DepMap_all_cell_line_classification.tsv"))
cls <- cls[grepl("^RNA", cls$modality), ]
mito <- read.delim("data/mito_classification.tsv")

# KMT2A fusions (all lines) from DepMap fusion calls
fu <- read.csv(file.path(dd, "DepMap_26Q3/OmicsFusionFiltered.csv"), check.names = FALSE)
fu <- fu[fu$IsDefaultEntryForModel == "Yes" & fu$ModelID %in% aml$ModelID, ]
g1 <- sub(" .*$", "", fu$Gene1); g2 <- sub(" .*$", "", fu$Gene2)
km <- fu[g1 == "KMT2A" | g2 == "KMT2A", ]
km$partner <- ifelse(sub(" .*$", "", km$Gene1) == "KMT2A", sub(" .*$", "", km$Gene2), sub(" .*$", "", km$Gene1))
# not translocations: KMT2A--KMT2A (partial tandem duplication) and lncRNA / clone-ID partners
# such as RP11-770J1.8 (likely read-through artifacts)
km <- km[km$partner != "KMT2A" & !grepl("^(RP[0-9]+-|CT[ABCD]-|A[CLP][0-9]{6}\\.|LINC[0-9])", km$partner), ]
km_partner <- tapply(km$partner, km$ModelID, function(p) paste(sort(unique(p)), collapse = ";"))

# mutations: hotspot or damaging, for recurrent AML genes
aml_genes <- c("TP53", "FLT3", "NPM1", "NRAS", "KRAS", "PTPN11", "KIT", "IDH1", "IDH2", "DNMT3A",
               "TET2", "ASXL1", "RUNX1", "CEBPA", "WT1", "SRSF2", "U2AF1", "SF3B1", "STAG2", "BCOR")
mut_any <- function(file) {
  m <- read_rows(file.path(dd, "DepMap_26Q3", file), aml$ModelID)
  m <- m[m$IsDefaultEntryForModel == "Yes", ]
  sy <- sub(" \\(.*$", "", names(m))
  M <- sapply(aml_genes, function(g) { j <- which(sy == g)
    if (!length(j)) return(rep(0, nrow(m))); rowSums(as.matrix(m[, j, drop = FALSE]) > 0) })
  rownames(M) <- m$ModelID; M
}
mh <- mut_any("OmicsSomaticMutationsMatrixHotspot.csv")
md <- mut_any("OmicsSomaticMutationsMatrixDamaging.csv")
ids <- union(rownames(mh), rownames(md))
mut <- matrix(0, length(ids), length(aml_genes), dimnames = list(ids, aml_genes))
mut[rownames(mh), ] <- mut[rownames(mh), ] + mh
mut[rownames(md), ] <- mut[rownames(md), ] + md

# FLT3-ITD (juxtamembrane in-frame insertion / duplication, codons ~560-640) and TKD (D835 / I836 /
# N841 / Y842) are in neither mutation matrix, so mut_FLT3 misses them; call them from the variant table
va <- read_rows(file.path(dd, "DepMap_26Q3", "OmicsSomaticMutations.csv"), aml$ModelID)
if ("IsDefaultEntryForModel" %in% names(va)) va <- va[va$IsDefaultEntryForModel == "Yes", ]
aa_pos <- suppressWarnings(as.numeric(sub("^p\\.[A-Za-z]+([0-9]+).*$", "\\1", va$ProteinChange)))
is_flt3 <- va$HugoSymbol == "FLT3"
itd <- is_flt3 & grepl("inframe_insertion|inframe_duplication", va$VariantInfo) & !is.na(aa_pos) &
  aa_pos >= 560 & aa_pos <= 640
tkd <- is_flt3 & grepl("^p\\.(D835|I836|N841|Y842)", va$ProteinChange)
itd_models <- unique(va$ModelID[itd])
flt3_models <- unique(va$ModelID[itd | tkd])
mut[rownames(mut) %in% flt3_models, "FLT3"] <- mut[rownames(mut) %in% flt3_models, "FLT3"] + 1
cat("FLT3-ITD lines:", paste(line_of_model[itd_models], collapse = ", "), "\n")

ic <- match(aml$StrippedCellLineName, cls$cell_line)
samples <- data.frame(
  cell_line        = aml$StrippedCellLineName,
  study_id         = "DEPMAP_26Q3",
  sample_type      = "cell_line",
  participant_id   = aml$ModelID,
  cell_line_name   = aml$CellLineName,
  CCLEName         = aml$CCLEName,
  RRID             = aml$RRID,
  disease_oncotree = aml$OncotreeCode,
  OncotreeSubtype  = aml$OncotreeSubtype,
  population       = c(Adult = "adult", Pediatric = "pediatric")[aml$AgeCategory],
  age_years        = suppressWarnings(as.numeric(aml$Age)),
  sex              = tolower(aml$Sex),
  PrimaryOrMetastasis = aml$PrimaryOrMetastasis,
  MYC_score        = cls$MYC_score[ic],
  mTOR_score       = cls$mTOR_score[ic],
  stringsAsFactors = FALSE)
# RNA signature classes as TRUE/FALSE groupings (High = above the score median; NA = not classified)
samples$MYC_high  <- cls$MYC_group[ic] == "High"
samples$mTOR_high <- cls$mTOR_group[ic] == "High"
samples$MYC_mTOR_category <- ifelse(is.na(ic), NA, paste0("MYC-", cls$MYC_group[ic], " / mTOR-", cls$mTOR_group[ic]))
samples$Mito_high <- mito$Mito_group[match(aml$StrippedCellLineName, mito$cell_line_as_published)] == "High"
samples$sex[!samples$sex %in% c("female", "male")] <- "unknown"
# KMT2A: curated status where available, else DepMap fusion calls
fusion_kmt2a <- aml$ModelID %in% names(km_partner)
samples$KMT2A_rearranged <- ifelse(is.na(ic), fusion_kmt2a, cls$KMT2A_rearranged[ic])
samples$KMT2A_partner    <- ifelse(is.na(ic), km_partner[aml$ModelID], cls$KMT2A_partner[ic])
samples$KMT2A_source     <- ifelse(is.na(ic), "DepMap fusion calls", "curated classification")
disagree <- !is.na(ic) & as.logical(cls$KMT2A_rearranged[ic]) != fusion_kmt2a
if (any(disagree)) cat("KMT2A curated vs fusion calls disagree for:",
                       paste(samples$cell_line[disagree], collapse = ", "), "\n")
# mutation columns: TRUE / FALSE, NA when the line has no mutation calls
mi <- match(aml$ModelID, rownames(mut))
for (g in aml_genes) samples[[paste0("mut_", g)]] <- ifelse(is.na(mi), NA, mut[pmax(mi, 1), g] > 0)
samples$mutations <- apply(mut[mi, , drop = FALSE], 1, function(r) paste(aml_genes[!is.na(r) & r > 0], collapse = ";"))
samples$FLT3_ITD <- ifelse(is.na(mi), NA, aml$ModelID %in% itd_models)
samples <- samples[order(samples$cell_line), ]

write_tsv(samples, "samples.tsv")
add("samples.tsv", "samples", "AML cell lines with harmonized annotation, groupings and mutation status")

# groupings the portal offers for two-group comparisons (volcano plots). Binary columns are
# TRUE/FALSE in samples.tsv; true_label / false_label are how the portal names the two groups.
grp <- read.delim(text = paste(
  "column\tlabel\ttrue_label\tfalse_label",
  "mTOR_high\tmTOR RNA signature\tmTOR-High\tmTOR-Low",
  "MYC_high\tMYC RNA signature\tMYC-High\tMYC-Low",
  "MYC_mTOR_category\tMYC x mTOR RNA category\t\t",
  "KMT2A_rearranged\tKMT2A rearrangement\tKMT2A-r\tKMT2A wild-type",
  paste0("mut_", aml_genes, "\t", aml_genes, " mutation\t", aml_genes, "-mutant\t", aml_genes, "-wild-type", collapse = "\n"),
  "FLT3_ITD\tFLT3-ITD\tFLT3-ITD\tno FLT3-ITD",
  "Mito_high\tMitochondrial group\tMito-High\tMito-Low",
  "disease_oncotree\tOncoTree code\t\t",
  "population\tAdult / pediatric\t\t",
  "sex\tSex\t\t",
  sep = "\n"), colClasses = "character", na.strings = NULL)
grp$kind <- ifelse(grp$true_label != "", "binary", "categorical")
# keep groupings with at least two groups of >= 3 lines
grp$n_groups_ge3 <- vapply(grp$column, function(c) sum(table(samples[[c]]) >= 3), 0L)
grp <- grp[grp$n_groups_ge3 >= 2, setdiff(names(grp), "n_groups_ge3")]
write_tsv(grp, "groupings.tsv")
add("groupings.tsv", "groupings", "Sample annotation columns offered as two-group comparisons")

# ============================================================ omics layers ==
layers <- list()
sym_entrez <- function(f) data.frame(
  gene_symbol = sub(" \\(.*$", "", f),
  entrez_id   = ifelse(grepl("\\(([0-9]+)\\)$", f), sub("^.*\\(([0-9]+)\\)$", "\\1", f), ""))
put <- function(key, X, features, label, modality, unit, source, lower_is_stronger = FALSE) {
  X <- X[rownames(X) %in% aml$StrippedCellLineName, , drop = FALSE]
  ok <- colSums(!is.na(X)) >= 3 & apply(X, 2, function(v) isTRUE(sd(v, na.rm = TRUE) > 0))
  X <- X[, ok, drop = FALSE]; features <- features[ok, , drop = FALSE]
  X <- X[order(rownames(X)), , drop = FALSE]
  write_tsv(data.frame(feature = colnames(X), features, check.names = FALSE),
            file.path("features", paste0(key, ".tsv")))
  write_tsv(data.frame(feature = colnames(X), t(signif(X, 5)), check.names = FALSE),
            file.path("values", paste0(key, ".tsv.gz")), gz = TRUE)
  add(file.path("features", paste0(key, ".tsv")), "features", paste(label, "- feature annotation"))
  add(file.path("values", paste0(key, ".tsv.gz")), "values", paste(label, "- values per cell line"))
  layers[[key]] <<- data.frame(key = key, label = label, modality = modality, unit = unit, source = source,
                               lower_is_stronger = lower_is_stronger, n_features = ncol(X),
                               n_lines = sum(rowSums(!is.na(X)) > 0))
  cat(sprintf("%-12s %3d lines x %6d features\n", key, layers[[key]]$n_lines, ncol(X)))
}

# --- RNA expression, DepMap 26Q3 (protein-coding genes, log2(TPM+1))
ex <- read_rows(file.path(dd, "DepMap_26Q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"), aml$ModelID)
ex <- ex[ex$IsDefaultEntryForModel == "Yes", ]
X  <- by_model(ex[, c("ModelID", names(ex)[7:ncol(ex)])], 1)
put("expression", X, data.frame(label = sub(" \\(.*$", "", colnames(X)), sym_entrez(colnames(X))),
    "RNA expression", "transcriptomics", "log2(TPM+1)", "DepMap 26Q3")

# --- CRISPR dependency, DepMap 26Q3 Chronos
cr <- read_rows(file.path(dd, "DepMap_26Q3/CRISPRGeneEffect.csv"), aml$ModelID)
X  <- by_model(cr, 1)
put("crispr", X, data.frame(label = sub(" \\(.*$", "", colnames(X)), sym_entrez(colnames(X))),
    "CRISPR dependency", "crispr_dependency", "Chronos gene effect", "DepMap 26Q3", TRUE)

# --- RNAi dependency, DEMETER2 v6 (combined)
d2 <- read_cols(file.path(dd, "DEMETER2_v6/D2_combined_gene_dep_scores.csv"), aml$CCLEName, ",")
X  <- t(as.matrix(d2[, -1])); colnames(X) <- d2[[1]]; rownames(X) <- line_of_ccle[rownames(X)]
put("rnai", X, data.frame(label = sub(" \\(.*$", "", colnames(X)), sym_entrez(colnames(X))),
    "RNAi dependency", "rnai_dependency", "DEMETER2 gene effect", "DEMETER2 v6", TRUE)

# --- Copy number, DepMap gene-level relative CN: log2(relative CN + 1), so 1 = ploidy-neutral.
# 26Q3 when downloaded, else the 24Q4 copy (source says which). Protein-coding genes with an NCBI
# cytoband only; chromosome + cytoband let the portal aggregate genes into genomic regions.
cn_file <- c(file.path(dd, "DepMap_26Q3/OmicsCNGene.csv"), file.path(dd, "DepMap_24Q4/OmicsCNGene.csv"))
cn_file <- cn_file[file.exists(cn_file)][1]
gi_file <- file.path(dd, "annotation/Homo_sapiens.gene_info.gz")
if (!is.na(cn_file) && file.exists(gi_file)) {
  cn <- read_rows(cn_file, aml$ModelID)
  if ("IsDefaultEntryForModel" %in% names(cn)) cn <- cn[cn$IsDefaultEntryForModel == "Yes", ]
  id_col <- if ("ModelID" %in% names(cn)) match("ModelID", names(cn)) else 1L   # 24Q4: unnamed first column
  gene_cols <- grep(" \\([0-9]+\\)$", names(cn))
  X  <- by_model(cn[, c(id_col, gene_cols)], 1)
  gi <- read.delim(gzfile(gi_file), quote = "", comment.char = "", stringsAsFactors = FALSE,
                   check.names = FALSE, colClasses = "character")
  gi <- gi[gi$type_of_gene == "protein-coding", ]
  i  <- match(sub("^.*\\(([0-9]+)\\)$", "\\1", colnames(X)), gi$GeneID)
  keep <- !is.na(i) & gi$map_location[i] != "-" & gi$chromosome[i] %in% c(1:22, "X", "Y")
  X <- X[, keep, drop = FALSE]; i <- i[keep]
  put("cn", X, data.frame(label = sub(" \\(.*$", "", colnames(X)), sym_entrez(colnames(X)),
                          chromosome = gi$chromosome[i], cytoband = gi$map_location[i]),
      "Copy number", "copy_number", "log2(relative CN + 1)",
      paste("DepMap", sub("^.*DepMap_([^/]+)/.*$", "\\1", cn_file)))
}

# --- Measured cell-surface proteome (George et al. 2026, bioRxiv 10.64898/2026.09.16.752145): surface
# biotinylation, mean abundance per line; 0 = not detected -> NA, then log2. Non-AML lines (ALL, CML) drop
# out in put(). surfaceome = UniProt places the protein at the cell surface (cell membrane, or signal
# peptide + transmembrane), when the annotation files are present; spc = the paper's surface score (0-4).
surface_genes <- local({
  files <- file.path(dd, "annotation", c("uniprot_cellmembrane.tsv", "uniprot_sp_tm.tsv"))
  files <- files[file.exists(files)]
  unique(unlist(lapply(files, function(f) read.delim(f, check.names = FALSE, quote = "")[["Gene Names (primary)"]])))
})
f_surf <- file.path(dd, "surfaceome_George_2026", "media-4_surface_mean_abundance.xlsx")
if (file.exists(f_surf) && requireNamespace("readxl", quietly = TRUE)) {
  x <- as.data.frame(readxl::read_excel(f_surf))
  lines <- setdiff(names(x), c("Accession", "Gene_Symbol", "SPC"))
  M <- as.matrix(x[, lines]); storage.mode(M) <- "double"; M[M <= 0] <- NA
  norm_line <- function(v) gsub("[^A-Z0-9]", "", toupper(v))
  alias <- c(WSU = "WSUAML")
  k <- norm_line(lines); k <- ifelse(k %in% names(alias), alias[k], k)
  X <- t(log2(M)); colnames(X) <- x$Accession
  rownames(X) <- aml$StrippedCellLineName[match(k, norm_line(aml$StrippedCellLineName))]
  ok <- !duplicated(x$Accession)
  X <- X[!is.na(rownames(X)), ok, drop = FALSE]
  sym <- sub("[;,].*$", "", ifelse(is.na(x$Gene_Symbol[ok]), "", x$Gene_Symbol[ok]))
  put("surface_george", X, data.frame(label = sym, gene_symbol = sym, entrez_id = "", spc = x$SPC[ok],
                                      surfaceome = sym %in% surface_genes, stringsAsFactors = FALSE),
      "Surface protein abundance (George 2026)", "proteomics", "log2 mean abundance",
      "George et al. 2026 bioRxiv 2026.09.16.752145")
}

# --- Drug response (AUC), DepMap harmonized screens
cmpd <- read.csv(file.path(dd, "DepMap_26Q3/PortalCompounds.csv"), stringsAsFactors = FALSE)
drug_files <- c(PRISM = "PRISM/REPURPOSINGAUCMatrix.csv", CTRP = "CTD2/CTRPAUCMatrix.csv",
                GDSC2 = "GDSC/GDSC2AUCMatrix.csv")
drug_src   <- c(PRISM = "PRISM Repurposing", CTRP = "CTD2 / CTRP", GDSC2 = "GDSC2")
for (s in names(drug_files)) {
  dr <- read_rows(file.path(dd, "DepMap_Portal_Harmonized_Drug_Screens", drug_files[[s]]), aml$ModelID)
  X  <- by_model(dr, 1)
  a  <- cmpd[match(colnames(X), cmpd$CompoundID), ]
  nm <- sub(" \\(.*\\)$", "", a$CompoundName)
  put(paste0("drug_", s), X,
      data.frame(label = ifelse(is.na(nm) | nm == "", colnames(X), nm), gene_symbol = "", entrez_id = "",
                 compound_name = a$CompoundName, targets = a$GeneSymbolOfTargets,
                 mechanism = a$TargetOrMechanism),
      paste0("Drug response (", s, ")"), "drug_screen", "AUC", drug_src[[s]], TRUE)
}

# --- RRBS methylation (CCLE): promoter (TSS +/- 1 kb) and CpG clusters
me <- read_cols(file.path(dd, "methylation_RRBS/CCLE_RRBS_TSS_1kb_20180614.txt"), aml$CCLEName, "\t", id_cols = 1:7)
X  <- t(as.matrix(me[, -(1:7)])); colnames(X) <- paste(me$gene, me$TSS_id, sep = "|")
rownames(X) <- line_of_ccle[rownames(X)]
put("meth_tss", X, data.frame(label = colnames(X), gene_symbol = me$gene, entrez_id = "", region = me$TSS_id),
    "Promoter methylation (TSS 1kb)", "methylation", "RRBS beta", "CCLE 2019")
mc <- read_cols(file.path(dd, "methylation_RRBS/CCLE_RRBS_TSS_CpG_clusters_20180614.txt"), aml$CCLEName, "\t", id_cols = 1:5)
X  <- t(as.matrix(mc[, -(1:5)])); colnames(X) <- paste(mc$gene_name, mc$cluster_id, sep = "|")
rownames(X) <- line_of_ccle[rownames(X)]
put("meth_cpg", X, data.frame(label = colnames(X), gene_symbol = mc$gene_name, entrez_id = "", region = mc$cluster_id),
    "CpG-cluster methylation", "methylation", "RRBS beta", "CCLE 2019")

# --- miRNA (CCLE 2019)
mr <- read_cols(file.path(dd, "CCLE_2019/CCLE_miRNA_20181103.gct"), aml$CCLEName, "\t", skip = 2, id_cols = 1:2)
X  <- log2(t(as.matrix(mr[, -(1:2)])) + 1); colnames(X) <- mr$Description
rownames(X) <- line_of_ccle[rownames(X)]
put("mirna", X, data.frame(label = colnames(X), gene_symbol = "", entrez_id = ""),
    "miRNA expression", "mirna", "log2(expression+1)", "CCLE 2019")

# --- Global chromatin profiling (CCLE 2019)
ch <- read.csv(file.path(dd, "CCLE_2019/CCLE_GlobalChromatinProfiling_20181130.csv"), check.names = FALSE)
ch <- ch[ch$CellLineName %in% aml$CCLEName, ]
X  <- as.matrix(ch[, -(1:2)]); rownames(X) <- line_of_ccle[ch$CellLineName]
put("chromatin", X, data.frame(label = colnames(X), gene_symbol = "", entrez_id = ""),
    "Global chromatin marks", "chromatin", "relative abundance", "CCLE 2019")

# --- Targeted metabolomics (Li et al. 2019, sheet "1-clean data")
read_xlsx_sheet <- function(dir, sheet_file) {
  ns <- c(d1 = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
  ss <- xml_text(xml_find_all(read_xml(file.path(dir, "xl/sharedStrings.xml")), "//d1:si", ns = ns))
  cells <- xml_find_all(read_xml(file.path(dir, "xl/worksheets", sheet_file)), "//d1:sheetData/d1:row/d1:c", ns)
  ref <- xml_attr(cells, "r"); typ <- xml_attr(cells, "t")
  val <- xml_text(xml_find_first(cells, "d1:v", ns))
  val[!is.na(typ) & typ == "s"] <- ss[as.integer(val[!is.na(typ) & typ == "s"]) + 1]
  col <- gsub("[0-9]", "", ref); row <- as.integer(gsub("[A-Z]", "", ref))
  colnum <- vapply(strsplit(col, ""), function(l) sum(match(l, LETTERS) * 26^(rev(seq_along(l)) - 1)), 0)
  M <- matrix(NA_character_, max(row), max(colnum)); M[cbind(row, colnum)] <- val
  M
}
xdir <- tempfile()
unzip(file.path(dd, "targeted_metabolomics/Supp_Tables.xlsx"),
      files = c("xl/sharedStrings.xml", "xl/worksheets/sheet3.xml"), exdir = xdir)
M <- read_xlsx_sheet(xdir, "sheet3.xml"); hdr <- M[1, ]; M <- M[-1, ]
X <- suppressWarnings(matrix(as.numeric(M[, -1]), nrow(M), dimnames = list(M[, 1], hdr[-1])))
X <- X[M[, 1] %in% aml$CCLEName, !is.na(hdr[-1]), drop = FALSE]
rownames(X) <- line_of_ccle[rownames(X)]
put("metabolomics", X, data.frame(label = colnames(X), gene_symbol = "", entrez_id = ""),
    "Metabolomics", "metabolomics", "log10 abundance", "Li et al. 2019 (CCLE)")

write_tsv(do.call(rbind, layers), "datasets.tsv")
add("datasets.tsv", "datasets", "Index of omics layers")

# ============================================================== provenance ==
commit <- suppressWarnings(tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
                                    error = function(e) character(0)))
dirty  <- length(suppressWarnings(system("git status --porcelain", intern = TRUE, ignore.stderr = TRUE))) > 0
write_tsv(data.frame(
  key   = c("prepared_at", "processing_commit", "working_tree_dirty", "r_version", "results_version"),
  value = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), if (length(commit)) commit else "",
            dirty, R.version.string, "3")), "provenance.tsv")
add("provenance.tsv", "provenance", "How and when these files were prepared")
write_tsv(do.call(rbind, manifest), "manifest.tsv")

# shared metadata the portal reads alongside the data
for (d in c("catalog", "schema")) { unlink(file.path("results", d), recursive = TRUE)
  file.copy(d, "results", recursive = TRUE) }
file.copy("dataset_sources.tsv", "results", overwrite = TRUE)
file.copy("docs/results_README.md", "results/README.md", overwrite = TRUE)
cat("done ->", out_dir, "\n")
