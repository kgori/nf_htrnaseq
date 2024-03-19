#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--dataset",
    help = "Path to the dataset", required = TRUE)
parser$add_argument("--gtf",
    help = "Path to the GTF file", required = TRUE)
args <- parser$parse_args()

if (!dir.exists(args$dataset)) {
    stop("The dataset path does not exist")
}

if (!file.exists(args$gtf)) {
    stop("The GTF file does not exist")
}

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
    library(dplyr)
})

extract_value <- function(input_strings, key, default = "") {
    re <- regexpr(paste0(key, ' "(.*?)";'), input_strings)
    ix <- re == -1
    lengths <- attr(re, "match.length")
    lengths[ix] <- 0
    re[ix] <- nchar(input_strings[ix])
    attr(re, "match.length") <- lengths

    values <- ifelse(re == -1, default, regmatches(input_strings, re))
    values <- gsub(paste0(key, ' "|";'), "", values)
    values[ix] <- default
    return(values)
}

extract_gene_id <- function(input_strings) {
    gene_ids <- extract_value(input_strings, "gene_id")
    return(gene_ids)
}

extract_gene_name <- function(input_strings, default = "") {
    gene_names <- extract_value(input_strings, "gene_name", default)
    return(gene_names)
}

# This function will join the dataset with the exons or CDS
# `join_with` is a data.table with the columns CHROM, START, and END
# used to find overlaps
collect_variants <- function(chrom, join_with) { # nolint start
    print(chrom)
    by <- join_by(CHROM, POS >= START, POS < END)
    dataset %>%
        filter(CHROM == chrom) %>%
        filter(is_somatic_filter_fail == FALSE & is_indel == FALSE) %>%
        arrange(CHROM, POS, REF, ALT) %>%
        select(CHROM, POS, REF, ALT, is_germline) %>%
        unique() %>%
        as.data.table() %>%
        inner_join(join_with, by = by) %>%
        select(CHROM, POS, REF, ALT, is_germline, gene_id, gene_name) %>%
        unique() %>%
        as.data.table()
} # nolint end

gtf <- fread(args$gtf)
setnames(gtf,
    c("CHROM", "source", "feature", "START",
        "END", "score", "strand", "frame", "attribute"))
exons <- gtf %>%
    mutate(CHROM = sub("^chr", "", CHROM)) %>%
    filter(feature == "exon") %>%
    select(CHROM, START, END, attribute) %>%
    arrange(CHROM, START, END) %>%
    unique(by = c("CHROM", "START", "END")) %>%
    mutate(gene_id = extract_gene_id(attribute)) %>%
    mutate(gene_name = extract_gene_name(attribute))
exons[gene_name == "", gene_name := gene_id]

cds <- gtf %>%
    mutate(CHROM = sub("^chr", "", CHROM)) %>%
    filter(feature == "CDS") %>%
    select(CHROM, START, END, attribute) %>%
    arrange(CHROM, START, END) %>%
    unique(by = c("CHROM", "START", "END")) %>%
    mutate(gene_id = extract_gene_id(attribute)) %>%
    mutate(gene_name = extract_gene_name(attribute))
cds[gene_name == "", gene_name := gene_id]

dataset <- open_dataset(args$dataset)
chromosomes <- dataset %>%
    select(CHROM) %>%
    unique() %>%
    pull(as_vector = TRUE) %>%
    sort()

# This will make a list of variant positions that intersect any exon
# These include both transcribed or UTR variants.
exome_variants <- lapply(chromosomes, collect_variants, join_with = exons) %>%
    rbindlist() %>%
    unique() %>%
    arrange(CHROM, POS, REF, ALT)

# This makes a list of variants that intersect any coding exon sequence
# These will not include any UTR variants, and should be a strict subset
# of the exome_variants list.
cds_variants <- lapply(chromosomes, collect_variants, join_with = cds) %>%
    rbindlist() %>%
    unique() %>%
    arrange(CHROM, POS, REF, ALT)

exome_variants[, CDS := FALSE]
exome_variants[
    cds_variants,
    CDS := TRUE,
    on = .(CHROM, POS, REF, ALT, gene_id, gene_name)]

exome_variants[, c("cds_gene_id", "cds_gene_name") := .("", "")]
exome_variants[, c("utr_gene_id", "utr_gene_name") := .("", "")]
exome_variants[CDS == TRUE, cds_gene_id := gene_id]
exome_variants[CDS == TRUE, cds_gene_name := gene_name]
exome_variants[CDS == FALSE, utr_gene_id := gene_id]
exome_variants[CDS == FALSE, utr_gene_name := gene_name]

# Collapse any positions that intersect more than one gene
exome_variants_collapsed_1 <- exome_variants[,
    .(gene_id = paste(gene_id, collapse = ","),
        gene_name = paste(gene_name, collapse = ",")),
    by = .(CHROM, POS, REF, ALT, is_germline)]
exome_variants_collapsed_2 <- exome_variants[CDS == TRUE,
    .(cds_gene_id = paste(cds_gene_id, collapse = ","),
        cds_gene_name = paste(cds_gene_name, collapse = ",")),
    by = .(CHROM, POS, REF, ALT, is_germline)]
exome_variants_collapsed_3 <- exome_variants[CDS == FALSE,
    .(utr_gene_id = paste(utr_gene_id, collapse = ","),
        utr_gene_name = paste(utr_gene_name, collapse = ",")),
    by = .(CHROM, POS, REF, ALT, is_germline)]

common_cols <- c("CHROM", "POS", "REF", "ALT", "is_germline")
exome_variants <- Reduce(
    function(x, y) merge(x, y, all = TRUE, by = common_cols),
    list(exome_variants_collapsed_2,
        exome_variants_collapsed_3,
        exome_variants_collapsed_1)
)

fwrite(exome_variants, "exome_variants.tsv.gz", sep = "\t", na = "")
