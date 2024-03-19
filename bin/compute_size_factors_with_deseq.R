#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--counts",
    help="Path to counts files",
    nargs='+',
    type="character")
parser$add_argument("--output",
    help="Path to output file",
    type="character")
args <- parser$parse_args()

for (filename in args$counts) {
    if (!file.exists(filename)) {
        stop(paste("File", filename, "does not exist"))
    }
}

suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(dplyr)
})

dt <- lapply(args$counts, function(filename) {
    fread(filename, skip = 1)
})
dt <- Reduce(function(x, y) {
    merge(x, y,
        by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
}, dt)

mtx <- dt[, as.matrix(.SD), .SDcols = seq(7, ncol(dt))]
rownames(mtx) <- dt$Geneid
colnames(mtx) <- basename(colnames(mtx))
col_data <- as.data.frame(colnames(mtx))

dds <- DESeqDataSetFromMatrix(countData = mtx,
    colData = col_data,
    design = ~1)
dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)
fwrite(
    data.table(samplename = names(size_factors),
        size_factor = size_factors)[order(samplename)],
    args$output,
    sep = "\t")
