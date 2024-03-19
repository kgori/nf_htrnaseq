#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--input", help = "Input file")
parser$add_argument("--output", help = "Output file")
args <- parser$parse_args()

if (!file.exists(args$input)) {
    stop("Input file does not exist")
}

# Taken directly from DESeq2, commit f6c515e R/methods.R:379
geomean <- function(x) {
    if (all(x == 0)) {
        NA_real_
    } else {
        exp(sum(log(x[x > 0])) / length(x))
    }
}

library(data.table)
dt <- fread(args$input)
dt[, empty := sum(TOTAL) == 0, by = .(CHROM, POS, REF, ALT)]

# To calculate size factors:
# (Optional) restrict to variants with presence in all tumours
# compute the geom.mean of totals for each variant (optional: add pseudocount)
# divide each variant's read total by the geometric mean for each variant,
# then take the median of these ratios for each sample

# Initial run to find any samples with a size factor of zero:
variants <- dt[empty == FALSE][, empty := NULL]
gmeans <- variants[, .(gmean = geomean(TOTAL + 1)),
    by = .(CHROM, POS, REF, ALT)]
setkey(gmeans, CHROM, POS, REF, ALT)
setkey(variants, CHROM, POS, REF, ALT)
variants[gmeans, gmean := i.gmean]
variants[, scaled_total := TOTAL / gmean]
size_factors <- variants[,
    .(size_factor = median(scaled_total, na.rm = TRUE)),
    by = samplename]
zero_factors <- size_factors[size_factor < 1e-6, samplename]
prev_size_factors <- copy(size_factors)
setorder(prev_size_factors, samplename)

# Recalculate size factors with these samples removed
variants <- dt[empty == FALSE & !(samplename %in% zero_factors)]
variants[, empty := NULL]
gmeans <- variants[, .(gmean = geomean(TOTAL + 1)),
    by = .(CHROM, POS, REF, ALT)]
setkey(gmeans, CHROM, POS, REF, ALT)
setkey(variants, CHROM, POS, REF, ALT)
variants[gmeans, gmean := i.gmean]
variants[, scaled_total := TOTAL / gmean]
size_factors <- variants[,
    .(size_factor = median(scaled_total, na.rm = TRUE)),
    by = samplename]
setorder(size_factors, samplename)

# Alternatively, restrict to variants with presence in all tumours
n_samples <- dt[, uniqueN(samplename)]
chosen <- dt[TOTAL > 0][, .N, by = .(CHROM, POS, REF, ALT)][N == n_samples]
variants <- dt[chosen, on = .(CHROM, POS, REF, ALT), nomatch = NULL]
gmeans <- variants[, .(gmean = geomean(TOTAL + 1)),
    by = .(CHROM, POS, REF, ALT)]
setkey(gmeans, CHROM, POS, REF, ALT)
setkey(variants, CHROM, POS, REF, ALT)
variants[gmeans, gmean := i.gmean]
variants[, scaled_total := TOTAL / gmean]
alt_size_factors <- variants[,
    .(size_factor = median(scaled_total, na.rm = TRUE)),
    by = samplename]
setorder(alt_size_factors, samplename)

output_stub <- tools::file_path_sans_ext(args$output)
fwrite(size_factors, paste0(output_stub, ".tsv"), sep = "\t")
fwrite(alt_size_factors, paste0(output_stub, ".alt.tsv"), sep = "\t")
