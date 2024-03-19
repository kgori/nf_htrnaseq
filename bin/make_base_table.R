#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--dataset",
    help = "Path to the dataset", required = TRUE)
parser$add_argument("--counts",
    help = "Path to the counts file", required = TRUE)
parser$add_argument("--variants",
    help = "Path to the variants file", required = TRUE)
parser$add_argument("--ht_regions",
    help = "Path to the HT regions file", required = TRUE)
parser$add_argument("--ancestral_states",
    help = "Path to the ancestral states file", required = TRUE)
parser$add_argument("--ctvt_a",
    help = "Path to the list of CTVT-A sample IDs", required = TRUE)
parser$add_argument("--output",
    help = "Path to the output file", required = TRUE)
args <- parser$parse_args()

# args <- list(dataset = "../../data/dataset",
#     counts = "../../runs/202402_RNA_Ribo/out/combined_counts.tsv.gz",
#     variants = "../../runs/202402_RNA_Ribo/out/exome_variants.tsv.gz",
#     ht_regions = "../../data/ht_regions.tsv",
#     ancestral_states = "../../runs/202402_RNA_Ribo/out/ancestral_genotypes.tsv",
#     ctvt_a = "../../data/ctvt_a.tsv",
#     output = "chk.txt")

if (!dir.exists(args$dataset)) {
    stop("The dataset path does not exist")
}

if (!file.exists(args$counts)) {
    stop("The counts file does not exist")
}

if (!file.exists(args$variants)) {
    stop("The variants file does not exist")
}

if (!file.exists(args$ht_regions)) {
    stop("The HT regions file does not exist")
}

if (!file.exists(args$ancestral_states)) {
    stop("The ancestral states file does not exist")
}

if (!file.exists(args$ctvt_a)) {
    stop("The CTVT-A sample IDs file does not exist")
}

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
    library(dplyr)
})

#' Reduce a samplename (e.g. 3145T1a) to the sample ID (e.g. 3145)
sample_id <- function(x) {
    .inner <- function(elem) {
        pattern <- "\\d+"
        match <- regexpr(pattern, elem)
        if (match == -1) {
            return(NA_character_)
        }
        regmatches(elem, match)
    }
    sapply(x, .inner)
}

ctvt_a <- fread(args$ctvt_a, colClasses = c("sample_id" = "character"))
regions <- fread(args$ht_regions, colClasses = c("CHROM" = "character"))
counts <- fread(args$counts) %>%
    inner_join(regions, by = join_by(CHROM, POS >= START, POS <= END)) %>%
    mutate(sample_id = sample_id(samplename))
variants <- fread(args$variants) %>%
    inner_join(regions, by = join_by(CHROM, POS >= START, POS <= END))
ancestral <- fread(args$ancestral_states,
    colClasses = c("CHROM" = "character", "HT" = "character")) %>%
    inner_join(regions, by = join_by(CHROM, POS >= START, POS <= END))
ancestral[CTVT == "0/0", CTVT := "AA"]
ancestral[CTVT == "0/1", CTVT := "AB"]
ancestral[CTVT == "1/0", CTVT := "AB"]
ancestral[CTVT == "1/1", CTVT := "BB"]
ancestral[CTVT == "0", CTVT := "A"]
ancestral[CTVT == "1", CTVT := "B"]
ancestral[CTVT == "./.", CTVT := NA_character_]
ancestral[HT == "0", HT := "A"]
ancestral[HT == "1", HT := "B"]

ds <- open_dataset(args$dataset)
dt <- ds %>%
    filter(is_germline) %>%
    inner_join(variants, by = c("CHROM", "POS", "REF", "ALT")) %>%
    relocate(samplename, hostname, CHROM) %>%
    mutate(nAH = round(expected_host_cn * (1 - H_vaf)),
        nBH = round(expected_host_cn * H_vaf)) %>%
    arrange(samplename, CHROM, POS, REF, ALT) %>%
    as.data.table() %>%
    mutate(sample_id = sample_id(samplename)) %>%
    mutate(is_ctvt_a = sample_id %in% ctvt_a$sample_id) %>%
    inner_join(counts, by = c("sample_id", "CHROM", "POS", "REF", "ALT")) %>%
    inner_join(ancestral, by = c("CHROM", "POS", "REF", "ALT")) %>%
    select(Tumour = samplename.y, Chrom = CHROM, Position = POS,
        Ref = REF, Alt = ALT, is_ctvt_a, nAH, nBH, GeneID = gene_id,
        Symbol = gene_name, copynumber = pipeline_cn,
        eA = NREF, eB = NALT, CTVT, HT)

fwrite(dt, file = args$output, sep = "\t", na = "NA")
