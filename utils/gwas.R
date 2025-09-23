#!/usr/bin/env Rscript
# Script: gwas.R
# Description: Performs GWAS for the genotype table from kcftools and phenotype data.
# Author: c.s.sivsubramani@gmail.com
# Date: 2025-09-18

options(warn = -1)

# ----------- Logging Helpers -----------
log_info <- function(msg) {
  cat(sprintf("[%s] [INFO] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

log_warn <- function(msg) {
  cat(sprintf("[%s] [WARN] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = stderr())
}

log_error <- function(msg) {
  cat(sprintf("[%s] [ERROR] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = stderr())
  quit(status = 1)
}

# ----------- Utilities ----------

# Function to plot Manhattan plot
manhattan_plot <- function(data, trait_name, pval_threshold = 0.05) {
  data <- data %>%
    arrange(CHR, START) %>%
    mutate(
      SNP = paste(CHR, START, sep = ":"),
      BPcum = cumsum(c(0, diff(START))),
      logP = -log10(pval)
    )

  bonforroni_correction <- pval_threshold / nrow(data)

  # Calculate cumulative position for each chromosome
  chr_lengths <- data %>%
    group_by(CHR) %>%
    summarise(chr_len = max(START)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(CHR, tot)

  data <- data %>%
    inner_join(chr_lengths, by = "CHR") %>%
    mutate(BPcum = START + tot)

  # Create x-axis labels
  axis_set <- data %>%
    group_by(CHR) %>%
    summarise(center = (min(BPcum) + max(BPcum)) / 2)

  # Plot
  p <- ggplot(data, aes(x = BPcum, y = logP, color = as.factor(CHR))) +
    geom_point(alpha = 0.75) +
    scale_color_manual(values = rep(c("grey", "black"), length.out = length(unique(data$CHR)))) +
    # geom_hline(yintercept = -log10(bonforroni_correction), color = "red") +
    geom_segment(aes(x = min(BPcum),
                     xend = max(BPcum),
                     y = -log10(bonforroni_correction),
                     yend = -log10(bonforroni_correction)),
                 inherit.aes = FALSE,
                 color = "red") +
    scale_x_continuous(label = paste0("Chr", axis_set$CHR), breaks = axis_set$center, expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, max(data$logP) + 2)) +
    labs(x = "Chromosome", y = "-log10(P-value)", title = trait_name) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(.25, "cm"),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text = element_text(color = "black", size = 12, face = "bold"),
    ) +
    coord_capped_cart(
      left   = capped_vertical(capped = "both")
    )
  return(p)
}

# Function to run GWAS
run_gwas <- function(gt, pheno, trait, cutoff_maf = 0.05, cutoff_miss = 0.5, cutoff_het = 0.2, out_dir = "gwas_res"){

  log_info(strrep("-", nchar(paste0("Running GWAS for trait: ", trait))))
  log_info(paste0("Running GWAS for trait: ", trait))
  log_info(strrep("-", nchar(paste0("Running GWAS for trait: ", trait))))

  # filter phenotype data to keep only relevant trait and non-missing values
  trait_data <- pheno %>%
    select(1, all_of(trait)) %>%
    filter(!is.na(.data[[trait]]))

  # get list of common accessions between genotype and phenotype data
  accessions <- intersect(trait_data[[1]], colnames(gt)[-c(1:4)])
  log_info(paste0(trait, " - Number of common accessions: ", length(accessions)))

  if (length(accessions) < 50) {
    log_warn(paste0(trait, " - Not enough common accessions between genotype and phenotype data. Skipping GWAS."))
    return(NULL)
  }

  # filter genotype data to keep only common accessions
  use_gt <- gt %>%
    select(1:4, all_of(accessions))

  # add a column MAF to genotype data
  log_info(paste0(trait, " - Calculating MAF, MISS, HET ..."))
  count_all <- nrow(use_gt)
  use_gt <- use_gt %>%
    filter(apply(select(., -c(1:4)), 1, function(x) length(unique(x))) > 1)
  count_poly <- nrow(use_gt)
  # summarize the genotype data
  use_gt <- calc_summary(use_gt)

  # filter the genotype data
  use_gt <- use_gt %>%
    filter(MAF >= cutoff_maf, MISS <= cutoff_miss, HET <= cutoff_het)
  count_flt <- nrow(use_gt)

  log_info(strrep("-", nchar(sprintf("%-40s : %8d", paste0(trait, " - Markers - input"), count_all))))
  log_info(sprintf("%-40s : %8d", paste0(trait, " - Markers - input"), count_all))
  log_info(sprintf("%-40s : %8d", paste0(trait, " - Markers - polymorphic"), count_poly))
  log_info(sprintf("%-40s : %8d", paste0(trait, " - Markers - Post Filtering"), count_flt))
  log_info(strrep("-", nchar(sprintf("%-40s : %8d", paste0(trait, " - Markers - input"), count_all))))
  # write.table(use_gt, paste0(out_dir, "/", trait, ".gt.flt.tsv"), sep = "\t", quote = F, row.names = F)
  if (nrow(use_gt) < 100) {
    log_warn(paste0(trait, " - Not enough SNPs after filtering. Skipping GWAS."))
    return(NULL)
  }

  use_trait <- trait_data %>%
    filter(taxa %in% accessions) %>%
    arrange(match(taxa, accessions))

  # write.table(use_trait, paste0(out_dir, "/", trait, ".pheno.tsv"), sep = "\t", quote = F, row.names = F)

  use_trait <- use_trait %>%
    select(all_of(trait))

  # convert trait data to numeric list
  use_trait <- as.numeric(unlist(use_trait))

  # Extract genotype matrix (drop annotation columns)
  gt_matrix <- as.matrix(use_gt[, -c(1:7)])

  # Extract genotype info (keep annotation columns)
  gt_info <- use_gt[, 1:7]

  # Transpose genotype matrix
  gt_matrix_t <- t(gt_matrix)

  # calculate kinship matrix using zhang_kinship function
  log_info(paste0(trait, " - Calculating kinship matrix ..."))
  kinmat <- zhang_kinship(gt_matrix_t)
  log_info(paste0(trait, " - plotting kinship matrix ..."))
  pheatmap::pheatmap(kinmat,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     main = "Kinship matrix",
                     fontsize = 6,
                     fontsize_row = 4,
                     fontsize_col = 4,
                     filename = paste0(out_dir, "/", trait, ".kinship.pdf"),
                     width = 10,
                     height = 8)
  log_info(paste0(trait, " - Writing kinship matrix to file ..."))
  write.table(kinmat, paste0(out_dir, "/", trait, ".kinship.tsv"), sep = "\t", quote = F)

  log_info(paste0(trait, " - Fitting linear mixed model to estimate heritability ..."))
  ID <- rownames(kinmat) ; length(ID)
  cbind(ID, use_trait)

  # Fit linear mixed model with kinship matrix as random effect
  mod <- lme4qtl::relmatLmer(use_trait ~ (1|ID), relmat = list(ID = kinmat))

  # Calculate narrow sense heritability
  herit.mod <- lme4qtl::VarProp(mod)
  herit.mod <- herit.mod$prop[1]
  # log_info(paste0(trait, " - Narrow-sense heritability (h²)   : ", round(herit.mod, 3)))
  log_info(sprintf("%-40s : %8.3f", paste0(trait, " - Narrow-sense heritability (h²)"), herit.mod))
  # Plot heritability estimate
  herit_df <- data.frame(
    component = c("Additive genetic", "Residual"),
    proportion = c(herit.mod, 1 - herit.mod)
  )

  p <- ggplot(herit_df, aes(x = "", y = proportion, fill = component)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("#1f78b4", "#b2df8a")) +
    labs(title = paste("Narrow-sense heritability for", trait),
         subtitle = paste("h² =", round(herit.mod, 3))) +
    theme_void()
  ggsave(p, filename = paste0(out_dir, "/", trait, ".heritability.piechart.pdf"), width = 6, height = 6)

  V <- lme4qtl::varcov(mod, idvar = "ID")
  V_thr <- V
  V_thr[abs(V) < 1e-10] <- 0
  decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
  W <- decomp$transform

  # make data object for mapping without any extra factors
  nnn <- rep(1,ncol(kinmat))

  logfile <- paste0(out_dir, "/", trait, ".gwas.log")

  out <- capture.output(
    gassoc_gls <- matlm(
      as.numeric(use_trait) ~ nnn,
      nnn,
      pred =  t(gt_matrix),
      ids = rownames(W),
      transform = W,
      batch_size = 4000,
      verbose = 2,
      cores = 1,
      stats_full = TRUE
    ),
    type = "output"
  )

  writeLines(out, con = logfile)

  # create a dataframe with use_info
  gwas_res <- gt_info %>%
    mutate(nobs = length(accessions)) %>%
    mutate(pval = gassoc_gls$tab$pval,
           zscore = round(gassoc_gls$tab$zscore * 10000, 2),
           se = round(gassoc_gls$tab$se * 10000, 2),
           b = round(gassoc_gls$tab$b * 10000, 2))
  write.table(gwas_res, paste0(out_dir, "/", trait, ".gwas.res.tsv"), sep = "\t", quote = F, row.names = F)
  if (all(is.na(gwas_res$pval))) {
    log_warn(paste0(trait, " - All p-values are NA. Skipping plotting."))
    return(NULL)
  }

  # manhattan plot
  p <- manhattan_plot(gwas_res, trait)
  ggsave(p, filename = paste0(out_dir, "/", trait, ".manhattan.plot.png"), width = 18, height = 6)

  # qq plot
  qq(gwas_res$pval,
     main = paste("QQ Plot for", trait),
     col = "blue", pch = 20,
     filename = paste0(out_dir, "/", trait, ".qq.plot.pdf"))


  log_info(strrep("-", nchar(sprintf("%-40s : %s",paste0(trait, " - Heritability pie chart saved to"),paste0(out_dir, "/", trait, ".heritability.piechart.pdf")))))
  log_info(sprintf("%-40s : %s", paste0(trait, " - GWAS results saved to"), paste0(out_dir, "/", trait, ".gwas.res.tsv")))
  log_info(sprintf("%-40s : %s", paste0(trait, " - Manhattan plot saved to"), paste0(out_dir, "/", trait, ".manhattan.plot.pdf")))
  log_info(sprintf("%-40s : %s", paste0(trait, " - QQ plot saved to"), paste0(out_dir, "/", trait, ".qq.plot.pdf")))
  log_info(sprintf("%-40s : %s", paste0(trait, " - Kinship matrix saved to"), paste0(out_dir, "/", trait, ".kinship.tsv")))
  log_info(sprintf("%-40s : %s", paste0(trait, " - Kinship matrix heatmap saved to"), paste0(out_dir, "/", trait, ".kinship.pdf")))
  log_info(sprintf("%-40s : %s", paste0(trait, " - Heritability pie chart saved to"), paste0(out_dir, "/", trait, ".heritability.piechart.pdf")))
  log_info(strrep("-", nchar(sprintf("%-40s : %s",paste0(trait, " - Heritability pie chart saved to"),paste0(out_dir, "/", trait, ".heritability.piechart.pdf")))))
  log_info(paste0(trait, " - GWAS analysis completed."))
  return(p)
}

# site_summary
calc_summary <- function(geno){
  # genotype matrix (rows = variants, cols = samples)
  geno_mat <- as.matrix(geno[, -c(1:4)])

  # replace "N" with -1 to keep numeric coding
  geno_mat[geno_mat == "N"] <- -1
  storage.mode(geno_mat) <- "numeric"

  n <- ncol(geno_mat)  # number of individuals

  # mask missing genotypes (-1) for calculations
  is_obs <- geno_mat != -1 & !is.na(geno_mat)

  # Allele counts
  ref_count <- rowSums((geno_mat == 0) & is_obs, na.rm = TRUE) * 2 +
    rowSums((geno_mat == 1) & is_obs, na.rm = TRUE)
  alt_count <- rowSums((geno_mat == 2) & is_obs, na.rm = TRUE) * 2 +
    rowSums((geno_mat == 1) & is_obs, na.rm = TRUE)
  total_alleles <- 2 * rowSums(is_obs)

  freq_alt <- ifelse(total_alleles > 0, alt_count / total_alleles, NA)
  maf <- pmin(freq_alt, 1 - freq_alt)

  # Missing rate
  miss <- rowMeans(!is_obs)

  # Heterozygosity
  het <- rowSums((geno_mat == 1) & is_obs, na.rm = TRUE) / rowSums(is_obs)

  # Bind results back, but keep original geno_mat with -1 intact
  use_gt_summary <- geno %>%
    select(ID, CHR, START, END) %>%
    mutate(
      MAF  = round(maf, 3),
      MISS = round(miss, 3),
      HET  = round(het, 3)
    ) %>%
    bind_cols(as.data.frame(geno_mat))

  return(use_gt_summary)
}


# Function to calculate kinship matrix using Zhang method (GAPIT inspired)
zhang_kinship <- function(snps,hasInbred=TRUE) {
  # Object: To calculate ZHANG (Zones Harbored Adjustments of Negligent Genetic) relationship
  # Authors: Zhwiu Zhang
  # Last update: october 25, 2014
  ##############################################################################################
  #Remove invariants
  fa=colSums(snps)/(2*nrow(snps))
  index.non=fa>=1| fa<=0
  snps=snps[,!index.non]

  het=1-abs(snps-1)
  ind.sum=rowSums(het)
  fi=ind.sum/(2*ncol(snps))
  inbreeding=1-min(fi)

  nSNP=ncol(snps)
  nInd=nrow(snps)
  n=nInd
  snpMean= apply(snps,2,mean)   #get mean for each snp
  # print("substracting mean...")
  snps=t(snps)-snpMean    #operation on matrix and vector goes in direction of column
  # print("Getting X'X...")
  #K=tcrossprod((snps), (snps))
  K=crossprod((snps), (snps))
  if(is.na(K[1,1])) stop ("Missing data is not allowed for numerical genotype data")

  #Extract diagonals
  i =1:n
  j=(i-1)*n
  index=i+j
  d=K[index]
  DL=min(d)
  DU=max(d)
  floor=min(K)

  #Set range between 0 and 2
  top=1+inbreeding
  K=top*(K-floor)/(DU-floor)
  Dmin=top*(DL-floor)/(DU-floor)

  #Adjust based on expected minimum diagonal (1)
  if(Dmin<1) {
    # print("Adjustment by the minimum diagonal")
    K[index]=(K[index]-Dmin+1)/((top+1-Dmin)*.5)
    K[-index]=K[-index]*(1/Dmin)
  }

  #Limiting the maximum offdiagonal to the top
  Omax=max(K[-index])
  if(Omax>top){
    # print("Adjustment by the minimum off diagonal")
    K[-index]=K[-index]*(top/Omax)
  }

  # print("Calculating kinship with Zhang method: done")
  return(K)
}

# ----------- Packages -----------
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressMessages({
  pkgs <- c("optparse", "dplyr", "ggplot2", "ggthemes", "readr", "lme4qtl", "matlm", "wlm", "qqman", "lemon", "purrr")
  missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs)) {
    log_error(paste("Missing required packages:", paste(missing_pkgs, collapse = ", "),
               "\nPlease install them before running this script."))
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
})


option_list <- list(
  make_option(c("-g", "--genotype"), type = "character", default = NULL,
              help = "Genotype file in TSV format (output from GATK, PLINK, etc.) [required]",
              metavar = "character"),
  make_option(c("-p", "--phenotype"), type = "character", default = NULL,
              help = "Phenotype file in TSV format (samples in rows, traits in columns) [required]",
              metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = "gwas_results",
              help = "Output directory [default= %default]", metavar = "character"),
  make_option(c("-m", "--maf"), type = "numeric", default = 0.05,
              help = "Minor allele frequency threshold [default= %default]", metavar = "numeric"),
  make_option(c("-r", "--miss"), type = "numeric", default = 0.1,
              help = "Missing rate threshold [default= %default]", metavar = "numeric"),
  make_option(c("-t", "--het"), type = "numeric", default = 0.1,
              help = "Heterozygosity rate threshold [default= %default]", metavar = "numeric")
)

if (length(commandArgs(trailingOnly = TRUE)) == 0)
  print_help(OptionParser(option_list = option_list))

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$genotype) || is.null(opt$phenotype) || is.null(opt$outdir))
  log_error("Missing required inputs: genotype, phenotype, and outdir")


log_info(strrep("-", nchar("CMD Options")))
log_info("CMD Options")
log_info(strrep("-", nchar("CMD Options")))

gt_file <- opt$genotype
pheno_file <- opt$phenotype
out_dir <- opt$outdir
maf_thresh <- opt$maf
miss_thresh <- opt$miss
het_thresh <- opt$het

log_info(sprintf("%-40s : %s", "Genotype file", gt_file))
log_info(sprintf("%-40s : %s", "Phenotype file", pheno_file))
log_info(sprintf("%-40s : %s", "Output directory", out_dir))
log_info(sprintf("%-40s : %8.3f", "MAF threshold", maf_thresh))
log_info(sprintf("%-40s : %8.3f", "Missing rate threshold", miss_thresh))
log_info(sprintf("%-40s : %8.3f", "Heterozygosity threshold", het_thresh))
log_info(strrep("-", nchar(sprintf("%-40s : %8.3f", "Heterozygosity threshold", het_thresh))))

dir.create(out_dir, showWarnings = FALSE)

gt <- read_tsv(gt_file, show_col_types = FALSE, comment = "#") %>%
  filter(CHR < 10)

count_all <- nrow(gt)
# remove monomorphic markers (all values are identical across samples)
gt <- gt %>%
  filter(apply(select(., -c(1:4)), 1, function(x) length(unique(x))) > 1)
count_poly <- nrow(gt)
log_info(sprintf("%-40s : %8d", paste0("Markers - Input"), count_all))
log_info(sprintf("%-40s : %8d", paste0("Markers - Accessions"), ncol(gt) - 4))
log_info(sprintf("%-40s : %8d", paste0("Markers - Polymorphic"), count_poly))
log_info(strrep("-", nchar(sprintf("%-40s : %8d", paste0("Markers - Polymorphic"), count_poly))))

if (ncol(gt) - 4 < 50) {
  log_error("Not enough accessions in genotype data. At least 50 accessions are required.")
}

if (nrow(gt) < 100) {
  log_error("Not enough polymorphic markers in genotype data. At least 100 markers are required.")
}


# kinship
genotype_matrix <- t(as.matrix(gt[,-c(1:4)]))
log_info("Calculating kinship matrix...")
kinship_matrix <- zhang_kinship(genotype_matrix, hasInbred = TRUE)
write_tsv(as.data.frame(kinship_matrix), file = paste0(out_dir, "/input.kinship.tsv"), col_names = TRUE)
log_info(paste0("Kinship matrix written to ", out_dir, "/input.kinship.tsv"))
log_info("Plotting kinship matrix...")
# heatmap of kinship
pheatmap::pheatmap(kinship_matrix,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   main = "Kinship matrix",
                   fontsize = 6,
                   fontsize_row = 4,
                   fontsize_col = 4,
                   filename = paste0(out_dir, "/input.kinship.pdf"),
                   width = 10,
                   height = 8)
rm(genotype_matrix, kinship_matrix)
log_info(paste0("Kinship matrix heatmap saved to ", out_dir, "/input.kinship.pdf"))

pheno <- read_tsv(pheno_file, show_col_types = FALSE)
log_info(strrep("-", nchar(sprintf("%40s : %8d", paste0("Traits in phenotype"), ncol(pheno) - 1))))
log_info(sprintf("%-40s : %-s", paste0("Phenotype file"), pheno_file))
log_info(sprintf("%-40s : %8d", paste0("Samples in phenotype"), nrow(pheno)))
log_info(sprintf("%-40s : %8d", paste0("Traits in phenotype"), ncol(pheno) - 1))
log_info(strrep("-", nchar(sprintf("%40s : %8d", paste0("Traits in phenotype"), ncol(pheno) - 1))))

list_of_traits <- colnames(pheno)[-1]
log_info(paste0("Traits to analyze: ", paste(list_of_traits, collapse = ", ")))

chr_lengths <- gt %>%
  group_by(CHR) %>%
  summarise(chr_len = max(START))

# Dataframe for Manhattan of Manhattan plot
gwas_hits <- chr_lengths %>%
  mutate(
    grids = map2(CHR, chr_len, ~ tibble(
      CHR = .x,
      pos = seq(0, .y, by = 5e5) / 1e6
    ))
  ) %>%
  pull(grids) %>%
  bind_rows()

plist <- list()
for (trait in list_of_traits){
  # Run GWAS for the trait
  p <- run_gwas(gt, pheno, trait, out_dir = out_dir, cutoff_maf = maf_thresh, cutoff_miss = miss_thresh, cutoff_het = het_thresh)
  if (is.null(p)) next
  # Check if results file exists
  if (!file.exists(paste0(out_dir, "/", trait, ".gwas.res.tsv"))) next
  plist[[trait]] <- p
  # Read results and filter by Bonferroni
  gwas_res <- read_tsv(paste0(out_dir, "/", trait, ".gwas.res.tsv"), show_col_types = FALSE)
  if (nrow(gwas_res) == 0) next
  bf_cutoff <- -log10(0.05/nrow(gwas_res))
  gwas_res <- gwas_res %>%
    filter(-log10(pval) >= bf_cutoff) %>%
    mutate(pos = round(START/500000, 0)/2) %>%
    group_by(CHR, pos) %>%
    slice_min(order_by = pval, n = 1) %>%
    ungroup() %>%
    select(CHR, pos, pval) %>%
    distinct() %>%
    mutate(!!trait := -log10(pval)) %>%
    select(CHR, pos, !!sym(trait))
  if (nrow(gwas_res) == 0) next
  gwas_hits <- merge(gwas_hits, gwas_res, by = c("CHR", "pos"), all.x = TRUE)
}

# print plist for each 6 traits in one page in pdf
if (length(plist) > 0) {
  pdf(paste0(out_dir, "/all_traits.manhattan.plots.pdf"), width = 6, height = 18) # adjust height for 6 plots per page

  for (i in seq(1, length(plist), by = 6)) {
    # take a chunk of up to 6 plots
    plots_chunk <- plist[i:min(i+5, length(plist))]

    # arrange plots in 2 columns and 3 rows
    do.call(gridExtra::grid.arrange, c(plots_chunk, ncol = 1, nrow = 6))
  }

  dev.off()
  log_info(paste0("All Manhattan plots saved to ", out_dir, "/all_traits.manhattan.plots.pdf"))
} else {
  log_warn("No GWAS results to plot.")
}


gwas_hits <- gwas_hits %>%
  tidyr::pivot_longer(
    cols = -c(CHR, pos),
    names_to = "trait",
    values_to = "Pvalue"
  ) %>%
  # right_join(all_positions, by = c("CHR", "pos")) %>%
  arrange(CHR, pos, trait) %>%
  filter(!is.na(trait))

write_tsv(gwas_hits, file = paste0(out_dir, "/all_traits.gwas_hits.tsv"))
log_info(paste0("GWAS hits written to ", out_dir, "/all_traits.gwas_hits.tsv"))

pval_max <- max(gwas_hits$Pvalue, na.rm = TRUE)
pval_max_rounded <- ceiling(pval_max / 5) * 5
pval_breaks <- seq(5, pval_max_rounded, by = 5)

log_info("Plotting Manhattan of Manhattan plot...")
p <- ggplot(gwas_hits, aes(x = pos, y = trait, size = Pvalue)) +
  geom_point(aes(color = ifelse(!is.na(Pvalue) & Pvalue > 0, "darkred", NA))) +
  geom_point(
    data = subset(gwas_hits, !is.na(Pvalue) & Pvalue > 0),
    color = "darkred"
  ) +
  facet_wrap(
    ~ CHR,
    nrow = 1,
    scales = "free_x",
    labeller = labeller(CHR = function(x) paste0("Chr ", x))
  ) +
  scale_size_continuous(
    range = c(1, 10),
    breaks = pval_breaks,
    limits = c(5, pval_max_rounded)
  ) +
  labs(x = "Position (Mb)", y = "Trait", size = "-log10(P)") +
  guides(color = "none") +
  theme_base() +
  theme(
    legend.position = "bottom",
    text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.ticks = element_line(color = "black", size = 0.7),
    strip.background = element_rect(fill = "grey", color = "black", size = 1),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.4, linetype = "dotted", color = "black")
  )

num_chr <- length(unique(gwas_hits$CHR))
num_traits <- length(unique(gwas_hits$trait))

plot_width <- max(8, num_chr * 3)
plot_height <- max(5, num_traits * 0.8)

ggsave(paste0(out_dir, "/", "all_traits.gwas_hits.pdf"),
       plot = p, width = plot_width, height = plot_height, units = "in")
log_info(paste0("Manhattan plot of all traits saved to ", out_dir, "/all_traits.gwas_hits.pdf"))
log_info("GWAS analysis completed.")
rm(gwas_hits, gwas_res, pheno)
