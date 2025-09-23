#!/usr/bin/env Rscript
# usage: Rscript hmp2gt.R <input.hmp.txt> <output.gt.tsv>
# description: Convert HapMap genotype file to numeric genotype format (Works only for SNPs)
# note: this script will only work for hapmap file in single nucl format (and only SNPs)
# date: 22-09-2025
# contact: c.s.sivasubramani@gmail.com

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

#---------------------------------------------
# Function to convert genotype alleles to numeric
#---------------------------------------------
convert_gt_to_numeric <- function(gt) {
  geno_cols <- 5:ncol(gt)        # columns after ID, CHR, START, alleles
  result <- gt[, 1:4]            # retain first 4 columns

  alleles <- gt$alleles
  gt_mat <- as.matrix(gt[, geno_cols])

  num_mat <- matrix(NA_integer_, nrow = nrow(gt), ncol = length(geno_cols))

  for (i in seq_len(nrow(gt))) {
    # Split ref and alt by '/'
    split_alleles <- strsplit(alleles[i], "/", fixed = TRUE)[[1]]
    ref <- split_alleles[1]
    alt <- split_alleles[2]

    row_vals <- gt_mat[i, ]

    # Map alleles to numeric values
    num_mat[i, ] <- vapply(row_vals, function(x) {
      if (x == ref) return(0L)
      else if (x == alt) return(2L)
      else if (x == "N") return(-1L)
      else return(1L)
    }, integer(1))
  }

  # Combine numeric matrix with original columns
  result <- cbind(result, as.data.frame(num_mat))
  colnames(result)[geno_cols] <- colnames(gt)[geno_cols]

  return(result)
}


args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript hmp2gt.R <input.hmp.txt> <output.gt.tsv>")
}

input_file <- args[1]
output_file <- args[2]

# Read input HapMap file
hmp <- read_tsv(input_file, show_col_types = FALSE)

# Select relevant columns and rename
gt <- hmp %>%
  # select(`rs#`, chrom, pos, alleles, everything()) %>%
  select(c(1:4), c(12:ncol(hmp))) %>%
  rename(ID = `rs#`, CHR = chrom, START = pos)

# Convert genotype to numeric
gt_numeric <- convert_gt_to_numeric(gt)

# Add END column, remove alleles, reorder columns
gt_numeric <- gt_numeric %>%
  mutate(END = START) %>%
  select(ID, CHR, START, END, everything(), -alleles)

# Write output
write_tsv(gt_numeric, output_file, na = "NA")
message("Conversion complete: ", output_file)

# EOF
