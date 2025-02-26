#!/usr/bin/env Rscript

# Author Chedly Kastally <ckastall@gmail.com>
# Version 1.1
# Copyright (C) 2024 Chedly Kastally <ckastall@gmail.com>
# Created on 2024-02-23 18:27
# Last updated on 2024-02-27 14:23

# Script to compute the Distribution of Spatial Variation in SNP Frequencies
# (DSVSF), the expected of heterozygosity and Fst (global and pairwise) from an
# output of PhyloGeoSim (v2.0.beta_210319)

# citation ----

# Kastally, C., Dellicour, S., Hardy, O. J., Gilbert, M., & Mardulyn, P. (2021).
# Estimating Migration of Gonioctena quinquepunctata (Coleoptera: Chrysomelidae) Inside a Mountain Range in a Spatially Explicit Context.
# Insect Systematics and Diversity, 5(5).
# https://doi.org/10.1093/isd/ixab019

# function ----

maf2dsvsf <- function(maf, thresh = 0.1) {
    dsvsf <- rep(NA, length(maf))
    dsvsf[maf <= thresh] <- 0
    dsvsf[maf > thresh & maf < (1 - thresh)] <- 1
    dsvsf[maf >= (1 - thresh)] <- 2
    return(dsvsf)
}

raw_pgs2mac <- function(data_df, retain_only_biallelic = T) {

    parse_pgs_singlepop_field <- function(x) {
        matrix(as.numeric(do.call(rbind, strsplit(as.character(x), "/"))), ncol = 4)
    }

    stopifnot("data.frame" %in% class(data_df))
    stopifnot(colnames(data_df)[1] == "ChromosomeID_locusID_fragmentID_snpID")

    result <- apply(data_df, 1, function(x) {
        locus_id <- x[1]
        mat_gt <- parse_pgs_singlepop_field(x[-1])
        totals <- colSums(mat_gt)
        res_mac <- NA
        if (sum(totals > 0) > 2) {
            cat("Non-bi-allelic locus detected: ", locus_id, " (Skipped)\n")
            res_mac <- NULL
        } else {
            index_alleles_secondIsMinor <- order(totals, decreasing = T)[2]
            res_mac <- mat_gt[, index_alleles_secondIsMinor]
        }
        return(list(locus_id, res_mac))
    })

    result_mac <- do.call(rbind, lapply(result, function(x) x[[2]]))
    list_ids <- sapply(result, function(x) x[[1]])
    which_retained_loci <- which(sapply(result, function(x) !(is.null(x[[2]]))))
    rownames(result_mac) <- list_ids[which_retained_loci]

    # We use the first row to have the sample size, we assume they are identical for all loci.
    # Note we do not care whether the locus considered is bi-allelic or not for computing the sample sizes.

    first_row <- parse_pgs_singlepop_field(data_df[1, -1])
    sample_size <- rowSums(first_row)

    all_results <- list()
    all_results$sample_size <- sample_size
    all_results$matrix_mac <- result_mac

    return(all_results)
}

fst_wc84 <- function(n_i = NA, p_i = NA, h_i = NA) {

    # Compute all pairwise Fst, following Weir & Cockerham 1984 (Evolution) page 1359
    # if *p_i* = the frequency of allele A in the sample of size *n_i* from population i (i = 1, 2, ..., r) and *h_i* is the observed proportion of individuals heterozygous for allele A, then:

    # NOTE: a better way would allow to send matrices
    stopifnot(is.numeric(n_i))
    stopifnot(is.numeric(p_i))
    stopifnot(length(n_i) > 1)
    stopifnot(length(n_i) == length(p_i))

    if (is.na(h_i)) {
        # if h_i is not provided, assumes HWE (Ho = He)
        h_i <- 2 * p_i * (1 - p_i)
    }

    r <- length(n_i)

    # n_bar = average sample size
    n_bar <- sum(c(n_i / r))

    n_c   <- ( (r * n_bar) - sum( c( n_i^2 / (r * n_bar))) ) / (r - 1)

    # p_bar = the average sample frequency of allele A
    p_bar <- sum( c( (n_i * p_i) / (r * n_bar) ) )

    # s_2 = the sample variance of allele A frequencies over populations
    s_2   <- sum( c( (n_i * (p_i - p_bar)^2) / ( (r - 1) *  n_bar)) )

    # h_bar = the average heterozygote frequency of allele A
    h_bar <- sum( c( (n_i * h_i) / (r * n_bar)) )

    a <- n_bar / n_c * (s_2 - ( (1 / (n_bar - 1)) * (p_bar * (1 - p_bar) - (r - 1) / r * s_2 - (1 / 4) * h_bar) ))
    b <- (n_bar / (n_bar - 1)) * ( p_bar * (1 - p_bar) - ( ((r - 1) / r) * s_2 - ( (2 * n_bar - 1) / (4 * n_bar) * h_bar ) ) )
    c <- h_bar / 2

    fst <- a / (a + b + c)
    result <- list(a, b, c, fst)
    return(result)
}

wrapper_fst <- function(matrix_maf, vector_sample_size) {

    if (length(vector_sample_size) == ncol(matrix_maf)) {
        list_res <- apply(matrix_maf, 1, function(x) {
                              fst_wc84(n_i = vector_sample_size, p_i = unname(x), h_i = NA)
                     })
    } else if (ncol(matrix_maf) == ncol(vector_sample_size)) {
        list_res <- lapply(seq_len(nrow(matrix_maf)), function(x) {

            c_af <- matrix_maf[x, ]
            c_ss <- vector_sample_size[x, ]

            fst_wc84(n_i = c_ss, p_i = unname(c_af), h_i = NA)
        })
    } else {
        stop("Problem during FST computation")
    }

    all_a <- (unlist(sapply(list_res, `[`, 1)))
    all_b <- (unlist(sapply(list_res, `[`, 2)))
    all_c <- (unlist(sapply(list_res, `[`, 3)))

    multi_locus_fst <- sum(all_a) / sum( (all_a + all_b + all_c) )

    return(multi_locus_fst)
}

process_vcf_lines <- function(x) {
    x <- strsplit(x, "\t")[[1]][-1:-9]
    x <- unlist(lapply(x, function(y) strsplit(y, ":")[[1]][1]))
    x <- gsub("\\|", "\\/", x)

    possible_gt <- c("./.", "0/0", "0/1", "1/1")

    if (!all(x %in% possible_gt)) {
        cat("locus with more than two alleles detected; dropping the allele\n")
        x <- NULL
    } else {
        x[x == "./."] <- NA
        x[x == "0/0"] <- "0"
        x[x == "0/1"] <- "1"
        x[x == "1/1"] <- "2"
        x <- as.numeric(x)
    }

    return(x)
}

get_vcf_locus_id <- function(x) {
    x <- strsplit(x, "\t")[[1]][1:2]
    locus_id <- paste0(x, collapse = "_")
    return(locus_id)
}

my_outer <- function(x) {
    Reduce("%o%", x)
}

hypergeom_projection_with_thresholds <- function(c_ac, c_an, c_post_n, threshold) {

    stopifnot(is.numeric(c(c_ac, c_an, c_post_n, threshold)))

    res <-
        dhyper(c(0:c_post_n),
               m = c_ac,
               n = (c_an - c_ac),
               k = c_post_n)

    which_0 <- which(0:c_post_n < threshold * c_post_n)
    which_1 <- which(0:c_post_n >= threshold * c_post_n & 0:c_post_n <= (1 - threshold) * c_post_n)
    which_2 <- which(0:c_post_n > (1 - threshold) * c_post_n)

    c("0" = sum(res[which_0]), "1" = sum(res[which_1]), "2" = sum(res[which_2]))
}

# Optparse ----

verbose_mode <- F

# Default parameters
default_threshold              <- 0.10
default_vcf_missing_data       <- "randomSample"
default_summary_statistic_mode <- TRUE

library("optparse")

option_list <-
    list(
         make_option(c("-i", "--input_file"),
         type    = "character",
         default = NA,
         help    = "Input file [Default: %default]"),
         make_option(c("-m", "--input_mode"),
         type    = "character",
         default = NA,
         help    = "Input mode: pgs or vcf [Default: %default]"),
         make_option(c("-n", "--posterior_n_file"),
         type    = "character",
         default = NA,
         help    = "File with posterior sample sizes, required if using a vcf [Default: %default]"),
         make_option(c("-p", "--population_file"),
         type    = "character",
         default = NA,
         help    = "population file, required if using a vcf [Default: %default]"),
         make_option(c("-s", "--vcf_missing_data"),
         type    = "character",
         default = "randomSample",
         help    = "Method to handle missing data: randomSample or proj [Default: %default]"),
         make_option(c("-t", "--threshold"),
         type    = "numeric",
         default = NA,
         help    = "Threshold used to compute the DSVSF [Default: %default]"),
         make_option(c("-o", "--output_prefix"),
         type    = "character",
         default = "Output_",
         help    = "Output prefix [Default: %default]")
         )

parser <-
    OptionParser(usage = "%prog -i input_file [-t 0.1] -m pgs [-o output_prefix]",
                 option_list = option_list)

opt <- parse_args(parser)

input_file       <- opt$input_file
output_prefix    <- opt$output_prefix
threshold        <- opt$threshold
input_mode       <- opt$input_mode
posterior_n_file <- opt$posterior_n_file
population_file  <- opt$population_file
vcf_missing_data <- opt$vcf_missing_data

## summary statistic mode
#
# The script estimate He, Fst and pwFst for a vcf file as well if the summary
# statistic mode is activated, the statistics are computed on all SNPs and all
# available data (ie, we treat correctly missing data, but we do not do any
# resampling; if we deactivate this mode, then we use the same random subsample
# as used for the computation of the DSVSF

summary_statistic_mode <- default_summary_statistic_mode

# setup ----

stopifnot(is.numeric(threshold))

if (!file.exists(input_file)) {
    cat("Input file is not valid, aborting.\n")
    stop("No or incorrect input file\n")
}

if (input_mode == "pgs") {
    pgs_mode <- T
    vcf_mode <- F
} else if (input_mode == "vcf") {
    pgs_mode <- F
    vcf_mode <- T
} else {
    cat("Input mode entered not recognized: ", input_mode, "\nAborting.\n")
    stop("Unrecognized input mode")
}

if (vcf_missing_data == "randomSample") {
    proj_mode <- F
    randomSample_mode <- T
} else if (vcf_missing_data == "proj") {
    proj_mode <- T
    randomSample_mode <- F
} else {
    cat("VCF missing data mode entered not recognized: ", vcf_missing_data, "\nAborting.\n")
    stop("Unrecognized VCF missing data mode")
}

# main ----

if (vcf_mode) {

    ## Only produces the DSVSF
    stopifnot(file.exists(population_file))
    stopifnot(file.exists(posterior_n_file))

    ## ploidy: assumed to be 2
    ploidy <- 2

    ### Handling the posterior n file
    post_n_dat <- read.table(posterior_n_file, sep = "\t", header = T, strip.white = T)

    expected_colnames_postn <- c("population_id", "posterior_n")
    stopifnot(colnames(post_n_dat) == expected_colnames_postn)
    stopifnot(is.numeric(post_n_dat$posterior_n))

    pop2post_n <- setNames(post_n_dat$posterior_n, post_n_dat$population_id)

    ### Handling the population data file
    pop_dat <- read.table(population_file, header = T, sep = "\t", strip.white = T)

    expected_pop_colnames <-
        c("sample_id", "population_id")

    stopifnot(colnames(pop_dat) == expected_pop_colnames)
    stopifnot(all(sort(unique(pop_dat$population_id)) == sort(post_n_dat$population_id)))

    ## Check the sample size in the population file: if it's below the one in the posterior sample size, there's an issue.
    sample_count_per_pop      <- as.data.frame(table(pop_dat$population_id) * ploidy)
    check_sample_size_in_post <- (pop2post_n[as.character(sample_count_per_pop$Var1)] <= sample_count_per_pop$Freq)

    if (!all(check_sample_size_in_post)) {
        cat(sprintf("The following populations have a posterior sample size above the number of samples included in the population file; that should not happen.\n%s\n",
                    paste0(names(which(!check_sample_size_in_post)), collapse = ", ")))
        stop("posterior sample sizes should be below the sample size inferred from the population data file\n")
    }

    sample2pop <-
        setNames(
            pop_dat$population_id,
            pop_dat$sample_id
        )

    ### Handling the vcf file

    expected_metadata_fields <-
        c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

    vlines <- readLines(input_file)

    header <- grep("^#CHR", vlines, value = T)
    header <- strsplit(header, "\t")[[1]]

    stopifnot(header[1:9] == expected_metadata_fields)

    sample_ids <- header[-1:-9]

    stopifnot(all(sample_ids %in% pop_dat$sample_id))
    unique_pops <- sort(unique(sample2pop[sample_ids]))

    list_index_pop <- lapply(unique_pops, function(x) which(sample2pop[sample_ids] == x))

    names(list_index_pop) <- unique_pops

    sample_size_pop_vcf <- sapply(list_index_pop, length) * ploidy

    check_sample_size_vcf <- sample_size_pop_vcf < pop2post_n[names(sample_size_pop_vcf)]

    if (any(check_sample_size_vcf)) {
        cat(sprintf("Number of samples of the following populations are smaller in the vcf than the requested posterior sample size; that should not be the case\n%s\n",
                paste0(names(check_sample_size_vcf), collapse = ", ")))
        stop("the number of samples in the vcf for each population should be higher or equal the posterior sample size requested for that population")
    }

    vlines <- grep("^#", vlines, invert = T, value = T)

    vlines_list <- lapply(vlines, process_vcf_lines)
    locus_id    <- unlist(lapply(vlines, get_vcf_locus_id))
    vcf_mat     <- do.call(rbind, vlines_list)

    stopifnot(length(sample_ids) == ncol(vcf_mat))

    colnames(vcf_mat) <- sample_ids

    if (randomSample_mode) {

        list_allSNPs_ac_pop <-
            lapply(seq_len(nrow(vcf_mat)), function(i) {

                if (verbose_mode) { cat("i: ", i, "\n") }

                x <- vcf_mat[i, ]

                x_pop <- split(x, sample2pop[names(x)])

                observed_sample_size <- sapply(x_pop, function(a) sum(!is.na(a)) * ploidy)
                check_sample_size_in_locus <- observed_sample_size < pop2post_n[names(observed_sample_size)]

                result <- NULL

                if (any(check_sample_size_in_locus)) {
                    cat(sprintf(
                        "post_n too high for locus %s in population(s): %s\nDropping locus %s\n", i, paste0(names(check_sample_size_in_locus), collapse = ", "), i
                        ))

                } else {

                    list_ac_pop <- lapply(unique(sample2pop), function(pop) {

                        # cat("pop: ", pop, "\n")

                        c_x_pop <- x_pop[[pop]]
                        post_n  <- pop2post_n[pop]

                        c_x_pop <- unlist(lapply(c_x_pop[!is.na(c_x_pop)], function(t) {

                            if (t == 0) {
                                res <- list(0, 0)
                            } else if (t == 1) {
                                res <- list(0, 1)
                            } else if (t == 2) {
                                res <- list(1, 1)
                            }

                            return(res)
                        }))


                        ac <- sum(sample(c_x_pop, post_n, replace = F))
                        names(ac) <- pop
                        return(ac)
                    })

                    result <- unlist(list_ac_pop)
                }

                return(result)
            })

        in_ac <- do.call(rbind, list_allSNPs_ac_pop)
        in_sample_size <- pop2post_n[colnames(in_ac)]
        total_sample_size <- sum(in_sample_size)
        in_af <- t(t(in_ac) / in_sample_size)

        stopifnot(all(in_af <= 1))

        # Use the minor allele at global scale
        change_polarity_i <- apply(in_ac, 1, function(x) sum(x) > total_sample_size / 2)

        in_maf <- in_af
        in_maf[which(change_polarity_i), ] <- 1 - in_maf[which(change_polarity_i), ]

        in_mac <- in_ac
        in_mac[which(change_polarity_i), ] <- t(in_sample_size - t(in_mac[which(change_polarity_i), ]))

        res_dsvsf                   <- t(apply(in_maf, 1, maf2dsvsf, thresh = threshold))
        res_dsvsf_2                 <- apply(res_dsvsf, 1, paste0, collapse = "")
        table_res_dsvsf_2           <- table(res_dsvsf_2)
        table_res_dsvsf_3           <- as.data.frame(sort(table_res_dsvsf_2, decreasing = T))
        colnames(table_res_dsvsf_3) <- c("DSVSF_pattern", "count")

        row.names(in_mac) <- locus_id
    }

    if (proj_mode) {

        # first "simplification", (but not used currently as it is fixed = 0), 
        # remove for a given group and given locus a projection that is below a threshold
        locus_level_simplification_threshold <- 0

        ## second simplification, after looking across groups for a single locus, remove patterns rarely seen
        second_level_thresh <- 1e-2

        list_pattern_counts <-
            lapply(seq_len(nrow(vcf_mat)), function(i) {

                x <- vcf_mat[i, ]

                list_proj <- lapply(unique(sample2pop), function(pop) {

                    i_pop <- which(sample2pop[names(x)] == pop)
                    x_pop <- x[i_pop]

                    an <- sum(!(is.na(x_pop))) * ploidy
                    ac <- sum(x_pop, na.rm = T)
                    post_n <- pop2post_n[pop]

                    res <- NULL
                    if (an >= post_n) {

                        res <- hypergeom_projection_with_thresholds(
                            c_ac = ac,
                            c_an = an,
                            c_post_n = post_n,
                            threshold = threshold
                        )
                        ## Add a limit to what we keep so that not all patterns are recorded

                        # Drop the empty indexes; it may be even wiser to drop anything "close" to 0; but right now I rely on a simplification later on
                        res <- res[res > locus_level_simplification_threshold]
                    } else {

                        cat("problem for i: ", i, " pop: ", pop, "\n")
                        cat("ac: ", ac, " an: ", an, "post_n: ", post_n, "\n")
                    }

                    return(res)
                })

                all_outer          <- my_outer(list_proj)

                patt_above_thresh_index <- which(all_outer >= second_level_thresh, arr.ind = T)

                pattern <- apply(patt_above_thresh_index, 1, paste0, collapse = "")
                count   <- all_outer[all_outer >= second_level_thresh]

                final_list <- list(pattern = unname(pattern), count = count)

                return(final_list)
            })

        all_patts <- sort(unique(unlist(sapply(list_pattern_counts, function(x) x[["pattern"]]))))

        list_all_patts_counts <- lapply(all_patts, function(c_pattern) {

                list_counts <- sapply(list_pattern_counts, function(x) {
                    c_count <- 0
                    index_pattern <- which(x[["pattern"]] == c_pattern)
                    if (length(index_pattern) > 0) {
                        c_count <- x[["count"]][index_pattern]
                    }
                    return(c_count)
                })

                sum(list_counts)
            })

        # adjust the patterns: they are here based on indexes, 1, 2, 3 whereas they should be 0, 1 and 2;

        all_patts <- gsub("2", "b", all_patts)
        all_patts <- gsub("3", "c", all_patts)

        all_patts <- gsub("1", "0", all_patts)
        all_patts <- gsub("b", "1", all_patts)
        all_patts <- gsub("c", "2", all_patts)

        final_table <- data.frame(pattern = all_patts, count = unlist(list_all_patts_counts))
        final_table <- final_table[order(final_table$count, decreasing = T), ]

        table_res_dsvsf_3           <- final_table
        colnames(table_res_dsvsf_3) <- c("DSVSF_pattern", "count")
    }

    ## Writing the dsvsf table
    dsvsf_output_table <- paste0(output_prefix, "dsvsf_output_table.tsv")
    write.table(table_res_dsvsf_3, dsvsf_output_table, sep = "\t", row.names = F, quote = F)


    if (summary_statistic_mode) {
    ## mode activated by default
    ## recompute the maf and mac on the entire sample instead of using the subsample (only done for random subsample actually)

        list_snps_ac_an <-
            lapply(seq_len(nrow(vcf_mat)), function(i) {

                # cat("i: ", i, "\n")

                x <- vcf_mat[i, ]

                x_pop <- split(x, sample2pop[names(x)])

                observed_sample_size <- sapply(x_pop, function(a) sum(!is.na(a)) * ploidy)
                observed_ac          <- sapply(x_pop, function(a) sum(a, na.rm = T))
                global_ac <- sum(observed_ac)
                global_an <- sum(observed_sample_size)
                observed_mac <- observed_ac
                if (global_ac / global_an > 0.5) {
                    observed_mac <- observed_sample_size - observed_ac
                }
                result <- list()
                result$mac <- observed_mac
                result$sample_size <- observed_sample_size

                return(result)
            })

        in_mac <- lapply(list_snps_ac_an, function(x) x$mac)
        in_mac <- do.call(rbind, in_mac)
        in_maf <- lapply(list_snps_ac_an, function(x) x$mac / x$sample_size)
        in_maf <- do.call(rbind, in_maf)
        in_sample_size <- lapply(list_snps_ac_an, function(x) x$sample_size)
        in_sample_size <- do.call(rbind, in_sample_size)

        stopifnot(all(in_maf <= 1))

    }
}

if (pgs_mode) {

    ## Handling the regular PGS file format

    in_dat <- read.table(input_file, header = T)

    in_parsed <- raw_pgs2mac(in_dat)
    in_mac <- in_parsed$matrix_mac
    in_sample_size <- in_parsed$sample_size
    in_maf <- t(t(in_mac) / in_sample_size)
    stopifnot(all(in_maf <= 1))

    # DSVSF ----

    res_dsvsf                   <- t(apply(in_maf, 1, maf2dsvsf, thresh = threshold))
    res_dsvsf_2                 <- apply(res_dsvsf, 1, paste0, collapse = "")
    table_res_dsvsf_2           <- table(res_dsvsf_2)
    table_res_dsvsf_3           <- as.data.frame(sort(table_res_dsvsf_2, decreasing = T))
    colnames(table_res_dsvsf_3) <- c("DSVSF_pattern", "count")

    dsvsf_output_table <- paste0(output_prefix, "dsvsf_output_table.tsv")
    write.table(table_res_dsvsf_3, dsvsf_output_table, sep = "\t", row.names = F, quote = F)
}

# Global He (across all populations) for each locus ----

matrix_of_He_per_pop <- 2 * (1 - in_maf) * in_maf
if (pgs_mode) {
    ## in pgs mode, we have the sample size fixed, no missing data
    loci_af              <- rowSums(in_mac) / sum(in_sample_size)

} else if (vcf_mode) {
    ## in vcf mode, we have missing data, unique sample size at each locus
    loci_af              <- rowSums(in_mac) / rowSums(in_sample_size)
    names(loci_af) <- locus_id

}

global_He            <- 2 * (1 - loci_af ) * loci_af

global_he_df <- data.frame(locus_id = names(global_He), He = unname(global_He))

if (pgs_mode) {
    colname_locus_id <- colnames(in_dat)[1]
} else {
    colname_locus_id <- "locus_id"
}

colnames(global_he_df)[1] <- colname_locus_id

global_he_output   <- paste0(output_prefix, "global_he_output.tsv")
write.table(global_he_df, global_he_output, sep = "\t", row.names = F, quote = F)

# Overall Fst ----

overall_fst <- wrapper_fst(in_maf, in_sample_size)

global_Fst_output  <- paste0(output_prefix, "global_Fst_output.txt")
cat("Overall Fst: ", overall_fst, "\n", file = global_Fst_output)

# Pairwise Fst ----

n_pops <- ncol(in_maf)
all_combn <- combn(n_pops, 2)

list_all_pwFst <- apply(all_combn, 2, function(x) {
          i <- x[1]
          j <- x[2]
          maf_ij <- in_maf[, c(i, j)]
          sample_size_ij <- in_sample_size[c(i, j)]
          overall_fst <- wrapper_fst(maf_ij, sample_size_ij)
          return(overall_fst)
                 })

result_pw_fst <- data.frame(cbind(t(all_combn), list_all_pwFst))
colnames(result_pw_fst) <- c("pop_i", "pop_j", "pwFst")

if (pgs_mode) {
    names_pop_i <- colnames(in_dat)[result_pw_fst$pop_i + 1]
    names_pop_j <- colnames(in_dat)[result_pw_fst$pop_j + 1]
} else {
    names_pop_i <- post_n_dat$population_id[result_pw_fst$pop_i]
    names_pop_j <- post_n_dat$population_id[result_pw_fst$pop_j]
}

result_pw_fst$pop_i <- names_pop_i
result_pw_fst$pop_j <- names_pop_j

all_pwFst_output   <- paste0(output_prefix, "all_pwFst_output.tsv")
write.table(result_pw_fst, all_pwFst_output, sep = "\t", row.names = F, quote = F)
