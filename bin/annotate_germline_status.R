#!/usr/bin/env Rscript

library(argparse)
library(logging)

parser <- ArgumentParser()
parser$add_argument("--counts", help = "Counts file", required = TRUE)
parser$add_argument("--variants", help = "Variants file", required = TRUE)
parser$add_argument("--output", help = "Output file", required = TRUE)
args <- parser$parse_args()

for (filename in c("counts", "variants")) {
    if (!file.exists(args[[filename]])) {
        stop(paste("File", args[[filename]], "does not exist"))
    }
}

library(data.table)
counts <- fread(args$counts, na.strings = "")
variants <- fread(args$variants, na.strings = "")

setkey(counts, CHROM, POS, REF, ALT)
setkey(variants, CHROM, POS, REF, ALT)
counts[variants, is_germline := i.is_germline]
setcolorder(counts, c("CHROM", "POS", "REF", "ALT", "is_germline"))
fwrite(counts, args$output, sep = "\t", na = "")
