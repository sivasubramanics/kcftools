#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript gt2plink.R <input.gt.tsv> <output.prefix>")
}

input_file <- args[1]
output_prefix <- args[2]

# input_file <- "lsal_ucd.sal.50k.gt.tsv"

# Read data, skip comment line starting with "#"
dat <- read_tsv(input_file, na = c("", "NA"), comment = "#", show_col_types = FALSE) %>% 
  # filter CHR to keep only numeric chromosomes (1-9)
  filter(CHR <= 9)

# Extract marker and genotype columns
marker_info <- dat[,1:4]
geno <- dat[,5:ncol(dat)]

# Create SNP ID
marker_info <- marker_info %>%
  mutate(SNP = paste(CHR, START, END, sep="_")) %>%
  select(CHR, SNP, START, END)

# Write MAP
map <- marker_info
colnames(map) <- c("CHR", "SNP", "CM", "BP")
map$CM <- 0
write.table(map, paste0(output_prefix, ".map"), quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

# Genotype recode function (0=AA, 1=AB, 2=BB, -1 missing)
convert_to_alleles <- function(x){
  ifelse(x == 0, "A A",
         ifelse(x == 1, "A B",
                ifelse(x == 2, "B B", "0 0")))
}

geno_alleles <- apply(geno, 2, convert_to_alleles)

# Build PED (individuals are columns → transpose)
ped <- as.data.frame(t(geno_alleles))

ped_data <- data.frame(
  FID = rownames(ped),
  IID = rownames(ped),
  PID = 0,
  MID = 0,
  SEX = 0,
  PHENOTYPE = -9,
  ped,
  check.names = FALSE
)

write.table(ped_data,  paste0(output_prefix, ".ped"), quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)

message("Conversion completed: ", output_prefix, ".ped and ", output_prefix, ".map")
