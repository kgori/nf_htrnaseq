#!/usr/bin/env Rscript
# Plotter of gene expressions of multiple categories
# ================================================== #
# Outputs: a couple of PDFs, for beeswarm and scatter plots
#          a few tables:
#          - normalised expressions (more or less a rehash of the input table)
#          - informative_sites (a table of the informative alleles)
#          - gene_stats (a table of statistics for each gene)

library(argparse)

parser <- ArgumentParser()
parser$add_argument("--table", help = "Input file",
    required = TRUE)
parser$add_argument("--size_factors", help = "Size factors file",
    required = TRUE)
parser$add_argument("--beeswarm_out",
    help = "Output file for beeswarm plots")
parser$add_argument("--scatterplot_out",
    help = "Output file for scatter plots")
parser$add_argument("--debug_table_out",
    help = "Output folder to write debug tables to")
args <- parser$parse_args()

if (!file.exists(args$table)) {
    stop("Input file does not exist")
}

if (!file.exists(args$size_factors)) {
    stop("Size factors file does not exist")
}

if (!endsWith(args$beeswarm_out, ".pdf")) {
    stop("Output file must be a PDF")
}

if (!endsWith(args$scatterplot_out, ".pdf")) {
    stop("Verbose output must be a PDF")
}

if (!dir.exists(args$debug_table_out)) {
    dir.create(args$debug_table_out)
}

library(dplyr)
library(data.table)

############################################
# Load the base table and the size factors #
############################################
deconv_table <- fread(args$table, na.strings = c("NA", "na", ""))
size_factors <- fread(args$size_factors, na.strings = c("NA", "na", ""))
size_factors[, Tumour := gsub(".bam$", "", samplename)]

##############
# eA and eB are observed counts of RNA reads
# eA_nrm and eB_nrm are normalised by size factor (subread + DESeq2)
deconv_table[size_factors,
    c("eA_nrm", "eB_nrm") := .(eA / size_factor, eB / size_factor),
    on = "Tumour"]

#########################################
# Plot expression of informative allele #
# for each relevant genotype            #
#########################################
find_informative_sites <- function(table) { # nolint start
    plot_data <- copy(table[,
        .(Tumour, Chrom, Position, Ref, Alt, is_ctvt_a,
            GeneID, Symbol, nAH, nBH,
            eA, eA_nrm, eB, eB_nrm, CTVT, HT)])

    # Generate host genotype
    plot_data[nAH == 0 & nBH == 0, Host := "-"]
    plot_data[nAH == 1 & nBH == 1, Host := "AB"]
    plot_data[nAH == 2 & nBH == 0, Host := "AA"]
    plot_data[nAH == 0 & nBH == 2, Host := "BB"]

    ###
    plot_data[is_ctvt_a == FALSE, HT := "-"]
    plot_data <- plot_data[!is.na(HT) & !is.na(CTVT) & HT != "" & CTVT != "?"]

    plot_data[, c("ht_only", "ht_and_ctvt", "ctvt", "host") := FALSE]
    plot_data[, paste("allele", c("ht_only", "ht_and_ctvt", "ctvt", "host"), sep = "_") := NA_character_]

    # Informative allele only in HT
    plot_data[HT == "A" & CTVT ==  "B" & Host == "BB", c("ht_only", "allele_ht_only") := .(TRUE, "A")]
    plot_data[HT == "A" & CTVT == "BB" & Host == "BB", c("ht_only", "allele_ht_only") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT ==  "A" & Host == "AA", c("ht_only", "allele_ht_only") := .(TRUE, "B")]
    plot_data[HT == "B" & CTVT == "AA" & Host == "AA", c("ht_only", "allele_ht_only") := .(TRUE, "B")]

    # Informative allele in HT and CTVT
    plot_data[HT == "A" & CTVT ==  "A" & Host == "BB", c("ht_and_ctvt", "allele_ht_and_ctvt") := .(TRUE, "A")]
    plot_data[HT == "A" & CTVT == "AA" & Host == "BB", c("ht_and_ctvt", "allele_ht_and_ctvt") := .(TRUE, "A")]
    plot_data[HT == "A" & CTVT == "AB" & Host == "BB", c("ht_and_ctvt", "allele_ht_and_ctvt") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT ==  "B" & Host == "AA", c("ht_and_ctvt", "allele_ht_and_ctvt") := .(TRUE, "B")]
    plot_data[HT == "B" & CTVT == "BB" & Host == "AA", c("ht_and_ctvt", "allele_ht_and_ctvt") := .(TRUE, "B")]
    plot_data[HT == "B" & CTVT == "AB" & Host == "AA", c("ht_and_ctvt", "allele_ht_and_ctvt") := .(TRUE, "B")]

    # Informative allele only in CTVT
    plot_data[HT == "-" & CTVT ==  "A" & Host == "BB", c("ctvt", "allele_ctvt") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT == "AA" & Host == "BB", c("ctvt", "allele_ctvt") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT == "AB" & Host == "BB", c("ctvt", "allele_ctvt") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT ==  "A" & Host == "BB", c("ctvt", "allele_ctvt") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT == "AA" & Host == "BB", c("ctvt", "allele_ctvt") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT == "AB" & Host == "BB", c("ctvt", "allele_ctvt") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT ==  "B" & Host == "AA", c("ctvt", "allele_ctvt") := .(TRUE, "B")]
    plot_data[HT == "-" & CTVT == "BB" & Host == "AA", c("ctvt", "allele_ctvt") := .(TRUE, "B")]
    plot_data[HT == "-" & CTVT == "AB" & Host == "AA", c("ctvt", "allele_ctvt") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT ==  "B" & Host == "AA", c("ctvt", "allele_ctvt") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT == "BB" & Host == "AA", c("ctvt", "allele_ctvt") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT == "AB" & Host == "AA", c("ctvt", "allele_ctvt") := .(TRUE, "B")]

    # Informative allele only in host
    plot_data[HT == "-" & CTVT ==  "B" & Host == "AA", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT ==  "B" & Host == "AB", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT == "BB" & Host == "AA", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT == "BB" & Host == "AB", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT ==  "B" & Host == "AA", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT ==  "B" & Host == "AB", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT == "BB" & Host == "AA", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "B" & CTVT == "BB" & Host == "AB", c("host", "allele_host") := .(TRUE, "A")]
    plot_data[HT == "-" & CTVT ==  "A" & Host == "AB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "-" & CTVT ==  "A" & Host == "BB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "-" & CTVT == "AA" & Host == "AB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "-" & CTVT == "AA" & Host == "BB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT ==  "A" & Host == "AB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT ==  "A" & Host == "BB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT == "AA" & Host == "AB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data[HT == "A" & CTVT == "AA" & Host == "BB", c("host", "allele_host") := .(TRUE, "B")]
    plot_data
} # nolint end

summarise_informative_alleles <- function(plot_data) { # nolint start
    meltonce <- melt(plot_data,
        id.vars = c("Tumour", "Chrom", "Position", "Ref", "Alt",
            "is_ctvt_a", "GeneID", "Symbol",
            "HT", "CTVT", "Host", "eA_nrm", "eB_nrm",
            "allele_ht_only", "allele_ht_and_ctvt", "allele_ctvt", "allele_host"),
        measure.vars = c("ht_only", "ht_and_ctvt", "ctvt", "host"),
        variable.name = "plot_category",
        value.name = "plot_is_active")[plot_is_active == TRUE][, plot_is_active := NULL]
    meltonce[plot_category == "ht_only", c("allele_ht_and_ctvt", "allele_ctvt", "allele_host") := NA]
    meltonce[plot_category == "ht_and_ctvt", c("allele_ht_only", "allele_ctvt", "allele_host") := NA]
    meltonce[plot_category == "ctvt", c("allele_ht_only", "allele_ht_and_ctvt", "allele_host") := NA]
    meltonce[plot_category == "host", c("allele_ht_only", "allele_ht_and_ctvt", "allele_ctvt") := NA]
    melttwice <- melt(meltonce, id.vars = c("Tumour", "Chrom", "Position", "Ref", "Alt",
        "is_ctvt_a", "GeneID", "Symbol",
        "HT", "CTVT", "Host",
        "eA_nrm", "eB_nrm", "plot_category"),
        measure.vars = c("allele_ht_only", "allele_ht_and_ctvt", "allele_ctvt", "allele_host"),
        variable.name = "allele_category",
        value.name = "informative_allele")[, allele_category := NULL][!is.na(informative_allele)]

    melttwice[, informative_allele_count := 0]
    melttwice[informative_allele == "A",
    informative_allele_count := ifelse(HT == "A", 1, 0) +
    ifelse(CTVT == "AB" | CTVT == "A", 1, ifelse(CTVT == "AA", 2, 0)) +
    ifelse(Host == "AB", 1, ifelse(Host == "AA", 2, 0))]

    melttwice[informative_allele == "B",
    informative_allele_count := ifelse(HT == "B", 1, 0) +
    ifelse(CTVT == "AB" | CTVT == "B", 1, ifelse(CTVT == "BB", 2, 0)) +
    ifelse(Host == "AB", 1, ifelse(Host == "BB", 2, 0))]

    melttwice[, expression := ifelse(informative_allele == "A",
        eA_nrm / informative_allele_count,
        eB_nrm / informative_allele_count)]
    melttwice
} # nolint end

infosites <- deconv_table %>%
    find_informative_sites()

selected_plot_data <- infosites %>%
    summarise_informative_alleles()

sample_mean_plot_data <- selected_plot_data[, .(expression = mean(expression)),
    by = .(plot_category, Tumour, GeneID, Symbol, informative_allele)]


###########################
# Homemade jitter functions
###########################
# Jitter within the kernel density (points go inside the lines of a violin plot)
# NB: no literal lines of a violin plot are drawn
violin_jitter <- function(data, dist = "uniform") {
    # KDE
    dens <- density(data)

    # Match data values to density bin
    bin <- findInterval(data, dens$x)

    # 'y' value of dens is the maximum amount of jitter to apply
    max_jitter <- dens$y[bin]

    # Sample from uniform [-1, 1] and multiply by max_jitter
    if (dist == "uniform") {
        jitter <- max_jitter * runif(length(data), -1, 1)
    } else {
        jitter <- max_jitter * rnorm(length(data))
    }

    data.frame(x = data, y = jitter)
}

beeswarm_jitter <- function(data, ptsize) { # nolint start
    # Dot plots, Leland Wilkinson, The American Statistician (1999)
    dt <- data.table(data = sort(data))
    dt[, grp := partition(data, ptsize)]
    dt[, y := 1:.N * ptsize, by = grp]
    dt[, x := sum(range(data)) / 2, by = grp]
    dt[, y := y - mean(y), by = grp]
    dt[, y := y[order(abs(y))], by = grp]
} # nolint end

# Beeswarm internal:
# inner recursive fn to split datapoints into groups
# point_j is grouped with point_i if:
#   point_i <= point_j < point_i + ptsize
#' @param v is vector of data points (assumed to be sorted
#'  in ascending order, but this is not checked. Since this
#'  function is only called from an outer partition function,
#'  which sorts its input, this should be OK)
#' @param ptsize numeric value
#' @param l list used to accumulate the groups. Gets passed into
#'  each recursive call.
.partition <- function(v, ptsize, l) {
    head <- v[v < v[1] + ptsize]
    tail <- v[v >= v[1] + ptsize]
    l <- append(l, list(head))
    if (length(tail) > 0) {
        return(.partition(tail, ptsize, l))
    } else {
        return(l)
    }
}

# Beeswarm internal:
# Transforms groups assigned by inner .partition fn into group index
#' @param v vector of data points.
partition <- function(v, ptsize) {
    sorted_order <- order(v)
    # Pass data in sorted order
    r <- .partition(v[sorted_order], ptsize, list())
    groups <- rep(seq_along(r), times = sapply(r, length))
    # Undo ordering, so group indices match the original order of the data
    groups[order(sorted_order)]
}

#' ptsize is passed to the beeswarm layout algorithm
plot_gene_data <- function(gene_data, axis_labels = FALSE, ptsize = 1) { # nolint start
    ymax <- max(10, gene_data[, max(expression)] * 1.05)
    ymax <- ceiling(ymax / 10) * 10
    gene <- gene_data[, Symbol[1]]
    plot(NA, xlim = c(0, 4), ylim = c(-ymax * 0.08, ymax),
        yaxt = "n", ylab = "Expression",
        xaxt = "n", xlab = NA,
        bty = "n", main = NA)

    r <- par()$usr
    w <- par()$fin[1] - sum(par()$mai[c(2, 4)])
    h <- par()$fin[2] - sum(par()$mai[c(1, 3)])

    title(main = gene, font.main = 3, adj = h * (1 / 10) / w,
        line = 0,
        cex.main = 0.9)

    if (axis_labels) {
        plot_labels <- c("HT", "CTVT", "Host")
    } else {
        plot_labels <- NA
    }

    axis(1, at = c(1:3), labels = plot_labels, mgp = c(3, 0.1, 0),
        lwd = 0, lwd.ticks = 1, cex.axis = .6, tck = -0.025)
    yticks <- c(0, ymax)
    axis(2, at = yticks, labels = yticks,
        lwd = 0, lwd.ticks = 1, las = 2)

    lines(c(r[1], r[2], r[1], r[1]), c(r[3], r[3], r[3], r[4]), lwd = 4)
    for (cat in sample_mean_plot_data[, sort(levels(plot_category))]) { #nolint start
        data <- gene_data[plot_category == cat, .(mean_expression = mean(expression)), by = .(Symbol, Tumour)][, sort(mean_expression)]
        if (length(data) > 0) {
            xpos <- if (cat == "ht_only") { 1 } else if (cat == "ctvt") { 2 } else { 3 }
            colour <- if (cat == "ht_only") {
                "gold"
            } else if (cat == "ctvt") {
                "grey70"
            } else {
                "#8ed3f3ff"
            }
            jitter <- beeswarm_jitter(sort(data), ptsize)$y
            jitter <- (jitter - mean(jitter))
            jitter <- jitter / (4 * diff(range(jitter)))
            jitter[is.na(jitter)] <- 0
            points(jitter + xpos, data, pch = 21, col = colour, bg = colour, cex = 0.9)
        }
    }
} # nolint end

genes <- intersect(sample_mean_plot_data[plot_category == "host", Symbol],
    sample_mean_plot_data[plot_category == "ht_only", Symbol])
genes <- sort(genes)

pdf(args$beeswarm_out, width = 8.25, height = 11.75 / 3)
par(mfrow = c(4, 5), mar = c(2.0, 1.75, 1.0, 1.25), oma = c(1,1,1,1))
for (gene in genes) {
    gene_data <- sample_mean_plot_data[Symbol == gene]
    plot_gene_data(gene_data, ptsize=1, axis_labels = TRUE)
}
dev.off()

#' `data` is sample_mean_plot_data object, from above
build_gene_stats_table <- function(gene, data, category = c("tumour", "host", "ht")) {
    category <- match.arg(category)
    gene_data <- data[Symbol == gene]
    expression_data <- if (category == "tumour") {
        gene_data[plot_category == "ctvt" | plot_category == "ht_and_ctvt",
            .(mean_expression = mean(expression)), by = .(Tumour, Symbol)][, mean_expression]
    } else if (category == "host") {
        gene_data[plot_category == "host", .(mean_expression = mean(expression)), by = .(Tumour, Symbol)][, mean_expression]
    } else {
        gene_data[plot_category == "ht_only", .(mean_expression = mean(expression)), by = .(Tumour, Symbol)][, mean_expression]
    }
    n_expr <- length(expression_data)
    mean_expr <- mean(expression_data)
    sd_expr <- sd(expression_data)
    sem_expr <- sd_expr / sqrt(n_expr)
    ci95_lower <- mean_expr - qnorm(0.975) * c(sem_expr)
    ci95_upper <- mean_expr + qnorm(0.975) * c(sem_expr)
    min_expr <- min(expression_data)
    max_expr <- max(expression_data)
    q25_expr <- quantile(expression_data, 0.25)
    q75_expr <- quantile(expression_data, 0.75)
    iqr_expr <- q75_expr - q25_expr
    data.table(gene = gene,
               category = category,
               n = n_expr,
               mean = mean_expr,
               sd = sd_expr,
               sem = sem_expr,
               ci95_l = ci95_lower,
               ci95_u = ci95_upper,
               min = min_expr,
               max = max_expr,
               q25 = q25_expr,
               q75 = q75_expr,
               iqr = iqr_expr)
}

gene_stats <- rbindlist(lapply(genes, function(gene) {
    rbindlist(lapply(c("tumour", "host", "ht"), function(cat) {
        build_gene_stats_table(gene, sample_mean_plot_data, cat)
    }))
}))

scatterplot <- function(data, xsample = "tumour", ysample = "ht") {
    ymax <- data[, max(ci95_u, na.rm = TRUE)]
    xy_data <- data[category %in% c(xsample, ysample)]
    xy_mean <- dcast(xy_data, gene ~ category, value.var = "mean")
    cc <- complete.cases(xy_mean)

    # Just plot the frame and axes
    xy_mean[, plot(get(ysample) ~ get(xsample), type = "n",
                xlim = c(0, ymax), ylim = c(0, ymax),
                xlab = xsample, ylab = ysample,
                bty = "n")]
    box(bty = "l", lwd = 3)
    
    # Fit a linear model and draw the regression line
    fit <- xy_mean[, lm(get(ysample) ~ get(xsample))]
    rsq <- summary(fit)$r.squared
    rho <- xy_mean[cc, cor(get(ysample), get(xsample))]
    print(summary(fit))
    do.call(abline, append(list(fit$coefficients),
                           list(col = "grey80")))
    
    # Add the points and error bars
    xy_mean[, points(get(ysample) ~ get(xsample),
        pch = 20, col = "grey20")]
    xy_data[, segments(x0 = mean[category == xsample],
        y0 = ci95_l[category == ysample],
        x1 = mean[category == xsample],
        y1 = ci95_u[category == ysample])]
    xy_data[, segments(x0 = ci95_l[category == xsample],
        y0 = mean[category == ysample],
        x1 = ci95_u[category == xsample],
        y1 = mean[category == ysample])]
    legend("topright",
        legend = c(sprintf("R-sq = %.3f", rsq),
            sprintf("rho = %.3f", rho)),
        bty = "n")
}

pdf(args$scatterplot_out, height = 8.25/2, width = 8.25/2)
scatterplot(gene_stats, "tumour", "ht")
scatterplot(gene_stats, "host", "ht")
scatterplot(gene_stats, "host", "tumour")
dev.off()

fwrite(deconv_table, file.path(args$debug_table_out, "debug_table_1_normalised.csv"))
fwrite(infosites, file.path(args$debug_table_out, "debug_table_2_infosites.csv"))
fwrite(selected_plot_data, file.path(args$debug_table_out, "debug_table_3_selected_plot_data.csv"))
fwrite(sample_mean_plot_data, file.path(args$debug_table_out, "debug_table_4_sample_mean_plot_data.csv"))
fwrite(gene_stats, file.path(args$debug_table_out, "debug_table_5_gene_stats.csv"))
