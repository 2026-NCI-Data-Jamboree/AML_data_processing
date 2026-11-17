# Convert a matrix or data frame saved as .rds to a gzipped TSV (first column "gene_id" = row names).
# Base R only; used by prepare_deconvolution.py for inputs Python cannot read (Beat AML's beatAML.rpkm.RDS).
#
# Usage: Rscript convert_rds_to_tsv.R <input.rds> <output.tsv.gz>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript convert_rds_to_tsv.R <input.rds> <output.tsv.gz>")
x <- readRDS(args[1])
if (is.null(rownames(x))) stop("The object has no row names to use as gene IDs")
con <- gzfile(args[2], "w")
write.table(data.frame(gene_id = rownames(x), as.data.frame(as.matrix(x)), check.names = FALSE),
            con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
