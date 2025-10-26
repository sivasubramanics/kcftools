#!/usr/bin/env Rscript
# # Script: plotIBS.R
# Description: Plot IBS windows is chromosome wise rect plot
# Author: c.s.sivsubramani@gmail.com
# Date: 2025-09-23

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggh4x)
  library(ape)
  # library(ggtree)
  # library(patchwork)
  # library(phytools)
})

# ── Command line options ──
option_list <- list(
  make_option(c("-c", "--chrinfo"), type="character", help="Chromosome metadata file"),
  make_option(c("-i", "--ibs"), type="character", help="Space-separated list of IBS summary files"),
  make_option(c("-o", "--output"), type="character", help="Output PDF file"),
  make_option(c("-g", "--groups"), type="character", default=NULL,
              help="Optional sample-to-group TSV file"),
  make_option(c("-t", "--tree"), type="character", default=NULL,
              help="Optional NJ tree file in Newick format for sample ordering"),
  make_option(c("-m", "--minlen"), type="numeric", default=1e6,
              help="Minimum length of IBS regions to plot (default: 1e6)"),
  make_option(c("-C", "--chrom"), type="character", default=NULL,
               help="Optional specific chromosome to plot (default: all)"),
  make_option(c("-V", "--var"), action="store_true", default=FALSE,
              help="Plot variable regions instead of identical regions (default: FALSE)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$chrinfo) || is.null(opt$ibs) || is.null(opt$output)) {
  print_help(OptionParser(option_list=option_list))
  stop("Missing required arguments.", call.=FALSE)
}

# # Debug test inputs
# opt$chrinfo <- "/Users/selva001/projects/work/wp5/lser/lser.chrinfo.tsv"
# opt$ibs <- "/Users/selva001/projects/work/wp5/lser/lser.lser.50k.ibs.summary.tsv"
# opt$output <- "/Users/selva001/projects/work/wp5/lser/lser.lser.50k.ibs.summary.plot.pdf"
# # opt$groups <- "/Users/selva001/projects/work/wp5/lser/countries.tsv"
# opt$tree <- "/Users/selva001/projects/work/wp5/lser/nj_tree.nwk"
# opt$minlen <- 1e6


# Split multiple IBS files
ibs_files <- strsplit(opt$ibs, "\\s+")[[1]] %>%
  purrr::map(~ Sys.glob(.x)) %>%
  unlist()


# ── Functions ──
parse_chr_meta <- function(chr_meta_file) {
  read.table(chr_meta_file, header = FALSE, sep = "\t") %>%
    setNames(c("chrom_name", "len", "chrom_num", "cum_len")) %>%
    mutate(
      chrom_num = as.numeric(chrom_num),
      chr_start = cum_len - len,
      chr_end = cum_len,
      mid = cum_len - (len / 2)
    )
}

# ── Load data ──
chrinfo <- parse_chr_meta(opt$chrinfo)

if (!is.null(opt$chrom)) {
  if (!(opt$chrom %in% chrinfo$chrom_name || as.numeric(opt$chrom) %in% chrinfo$chrom_num)) {
    stop(paste("Chromosome", opt$chrom, "not found in chromosome metadata."))
  }
  chrinfo <- chrinfo %>%
    filter(chrom_name == opt$chrom | chrom_num == as.numeric(opt$chrom))
  cat("Plotting only chromosome:", paste(chrinfo$chrom_name, collapse=", "), "\n")
}

ibs <- map_dfr(ibs_files, read_tsv, show_col_types = FALSE) %>%
  filter(Chromosome %in% chrinfo$chrom_name)

accessions <- ibs %>%
  distinct(Sample)

ibs <- ibs %>%
  filter(Length >= opt$minlen)

if (nrow(ibs) == 0) {
  stop("No IBS segments found after filtering. Try lowering --minlen or check input files.")
}
write.table(ibs, file=sub(".pdf$", ".filtered.tsv", opt$output), sep="\t", quote=FALSE, row.names=FALSE)

# ── Sample positions ──
# samples <- sort(unique(ibs$Sample))
samples <- sort(accessions$Sample)
sample_positions <- tibble(Sample = samples, SamplePosition = seq_along(samples))

# If grouping is provided
if (!is.null(opt$groups)) {
  groups <- read_tsv(opt$groups, col_names = c("Sample","Group"), show_col_types = FALSE) %>%
    mutate(Group = factor(Group, levels = unique(Group)))
  sample_positions <- sample_positions %>%
    left_join(groups, by="Sample") %>%
    arrange(Group, Sample) %>%
    mutate(SamplePosition = row_number())
} else {
  groups <- NULL
}

ibs_plot <- ibs %>%
  left_join(sample_positions, by="Sample")

chrom_levels <- chrinfo$chrom_name

sample_backgrounds <- expand_grid(
  # Sample = unique(ibs$Sample),
  Sample = unique(accessions$Sample),
  Chromosome = chrom_levels
) %>%
  left_join(sample_positions, by="Sample") %>%
  left_join(chrinfo %>% select(chrom_name, len),
            by = c("Chromosome" = "chrom_name")) %>%
  mutate(Start = 0, End = len,
         Chromosome = factor(Chromosome, levels=chrom_levels))

ibs_plot <- ibs_plot %>%
  mutate(Chromosome = factor(Chromosome, levels=chrom_levels))


# ── Tree-based ordering ──
tree <- NULL
if (!is.null(opt$tree)) {
  if (!file.exists(opt$tree)) {
    stop(paste("Tree file not found:", opt$tree))
  }
  cat("Reading NJ tree from", opt$tree, "...\n")

  # Read tree
  tree <- ape::read.tree(opt$tree)

  # Ladderize to match typical iTOL display (larger clades on top)
  tree <- ladderize(tree, right = TRUE)

  # Extract tip order from ladderized tree
  tree_order <- tree$tip.label[tree$edge[tree$edge[,2] <= length(tree$tip.label),2]]

  # Keep only tips present in the IBS data
  tree_order <- tree_order[tree_order %in% sample_positions$Sample]

  # Append missing samples at the bottom
  missing_samples <- setdiff(sample_positions$Sample, tree_order)
  if (length(missing_samples) > 0) {
    warning("Samples missing from tree and will be appended at bottom: ",
            paste(missing_samples, collapse = ", "))
    tree_order <- c(tree_order, missing_samples)
  }

  # Reorder sample_positions based on tree_order
  sample_positions <- sample_positions %>%
    mutate(Sample = factor(Sample, levels = tree_order),
           SamplePosition = as.integer(factor(Sample, levels = tree_order)))

  # Update ibs_plot SamplePosition directly
  ibs_plot <- ibs_plot %>%
    mutate(Sample = factor(Sample, levels = tree_order),
           SamplePosition = match(Sample, tree_order))

  # Drop any tree tips not in the final order and enforce tree_order on labels
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tree_order))
  tree$tip.label <- tree_order
}


# ── Calculate group separators ──
group_separators <- NULL
if (!is.null(groups)) {
  group_boundaries <- sample_positions %>%
    group_by(Group) %>%
    summarise(
      group_start = min(SamplePosition),
      group_end = max(SamplePosition),
      .groups = "drop"
    ) %>%
    arrange(group_start)

  # Create separator lines between groups (exclude the last group)
  if (nrow(group_boundaries) > 1) {
    group_separators <- group_boundaries[-nrow(group_boundaries), ] %>%
      mutate(separator_y = group_end + 0.5) %>%
      select(separator_y)
  }
}

# ── Plot ──
bandwidth <- 0.4
n_samples <- length(unique(sample_positions$Sample))
n_chr <- length(chrom_levels)

if (opt$var) {
  low_col <- "darkred"
  high_col <- "red"
} else {
  low_col <- "green"
  high_col <- "darkgreen"
}

p <- ggplot() +
  geom_rect(
    data = sample_backgrounds,
    aes(xmin = Start/1e6, xmax = End/1e6,
        ymin = SamplePosition - bandwidth,
        ymax = SamplePosition + bandwidth),
    fill = "grey99", color = NA
  ) +
  geom_rect(
    data = ibs_plot,
    aes(xmin = Start/1e6, xmax = End/1e6,
        ymin = SamplePosition - bandwidth,
        ymax = SamplePosition + bandwidth,
        fill = MeanScore),
    color = NA
  ) +
  scale_fill_gradient(
    high = high_col,
    low = low_col,
    name = "Mean Score",
    limits = c(min(ibs_plot$IBSProportion, na.rm=TRUE), max(ibs_plot$IBSProportion, na.rm=TRUE)),
    oob = scales::squish
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = sample_positions$SamplePosition,
    labels = sample_positions$Sample,
    expand = c(0.01, 0.01),
    sec.axis = if (!is.null(groups)) {
      dup_axis(
        breaks = sample_positions %>%
          group_by(Group) %>%
          summarise(mid = mean(SamplePosition)) %>%
          pull(mid),
        labels = sample_positions %>%
          group_by(Group) %>%
          summarise(mid = mean(SamplePosition)) %>%
          pull(Group),
        name = NULL
      )
    } else waiver()
  ) +
  scale_x_continuous(
    name = "Position (Mb)",
    # if n_chr > 1, set breaks every 50 Mb, else every 10 Mb
    breaks = if (n_chr > 1) seq(0, max(chrinfo$len)/1e6, by=50) else seq(0, max(chrinfo$len)/1e6, by=10),
    # breaks = seq(0, max(chrinfo$len, na.rm=TRUE)/1e6, by=50),
    labels = function(x) round(x, 0),
    expand = c(0.01, 0.01)
  ) +
  # ggh4x::facet_nested(. ~ Chromosome, scales="free_x", space="free_x", switch="y") +
  ggh4x::facet_nested(
    . ~ Chromosome,
    scales = "free_x",
    space = "free_x",
    switch = "y",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 28, face="bold"),
    axis.text.y = element_text(size = 18, face="bold", color="black"),
    axis.text.y.right = element_text(size = 28, face = "bold", color = "black", hjust = 0),
    # axis.text.x = element_text(size = 8, face="bold", angle=90, hjust=1),
    axis.text.x = element_text(size = 24, face="bold"),
    axis.title.x = element_text(size = 28, face="bold", margin=margin(t=10)),
    axis.ticks = element_line(color="black"),
    axis.line = element_line(color="black"),
    panel.spacing.x = unit(0.6,"lines"),
    panel.spacing.y = unit(0.6,"lines"),
    panel.border = element_rect(color="black", fill=NA, linewidth=6),
    legend.position = "none",
    # legend.position = "bottom",
    # legend.title = element_text(size=12, face="bold"),
    # legend.text = element_text(size=10),
    # legend.key.size = unit(0.8,"cm"),
    plot.margin = margin(10,10,10,10)
  )

# Add group separator lines if groups are provided
if (!is.null(group_separators) && nrow(group_separators) > 0) {
  p <- p + geom_hline(
    data = group_separators,
    aes(yintercept = separator_y),
    color = "black",
    linewidth = 2  # Same thickness as panel border
  )
}

# ── Auto dimensions ──
if (n_samples < 10) {
  plot_height <- 2.3 + (n_samples * 0.3)
} else {
  plot_height <- 2.3 + (n_samples * 0.2)
}
# plot_height <- max(2.3, (n_samples * 0.2) + 5)
if (n_chr > 1){
  plot_width  <- max(6, (n_chr * 4) + (n_samples * 0.2))  # add width scaling with samples
} else {
  # fix plot with based on size of the crhomesome
  plot_width  <- max(12, (max(chrinfo$len, na.rm=TRUE)/1e6 * 0.1) + (n_samples * 0.2))  # add width scaling with samples
}
ggsave(opt$output, plot = p, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
