#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--variants", help = "Path to the variants file",
    required = TRUE)
parser$add_argument("--regions", help = "Path to the regions file",
    required = TRUE)
parser$add_argument("--output", help = "Path to the output file",
    required = TRUE)
args <- parser$parse_args()

if (!file.exists(args$variants)) {
    stop("Variants file not found")
}

if (!file.exists(args$regions)) {
    stop("Regions file not found")
}

library(logging)
basicConfig()
loginfo("Arguments are valid")

library(data.table)
library(dplyr)

loginfo("Reading variants from ", args$variants)
variants <- fread(args$variants, na.strings = c("", "NA"))

loginfo("Reading regions from ", args$regions)
regions <- fread(args$regions, na.strings = c("", "NA"),
    colClasses = c("character", "integer", "integer", "character"))

loginfo("Restricting variants to regions")
restricted <- variants %>%
    inner_join(regions, by = join_by(CHROM, POS >= START, POS <= END))

loginfo("Writing restricted variants to ", args$output)
fwrite(restricted, args$output, sep = "\t")
