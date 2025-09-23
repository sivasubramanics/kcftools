#!/usr/bin/env Rscript
# Script: plotIBS.R
# Description: Plot IBS windows is chromosome wise rect plot
# Author: c.s.sivsubramani@gmail.com
# Date: 2025-09-23

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggh4x)
})

# ── Command line options ──
option_list <- list(
  make_option(c("-c", "--chrinfo"), type="character", help="Chromosome metadata file"),
  make_option(c("-i", "--ibs"), type="character", help="Space-separated list of IBS summary files"),
  make_option(c("-o", "--output"), type="character", help="Output PDF file"),
  make_option(c("-g", "--groups"), type="character", default=NULL,
              help="Optional sample-to-group TSV file"),
  make_option(c("-m", "--minlen"), type="numeric", default=1e6,
              help="Minimum length of IBS regions to plot (default: 1e6)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$chrinfo) || is.null(opt$ibs) || is.null(opt$output)) {
  print_help(OptionParser(option_list=option_list))
  stop("Missing required arguments.", call.=FALSE)
}

# Split multiple IBS files
ibs_files <- strsplit(opt$ibs, "\\s+")[[1]] %>%
  map(~ Sys.glob(.x)) %>%
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

ibs <- map_dfr(ibs_files, read_tsv, show_col_types = FALSE) %>%
  filter(Chromosome %in% chrinfo$chrom_name) %>%
  filter(Length >= opt$minlen)

write.table(ibs, file=sub(".pdf$", ".filtered.tsv", opt$output), sep="\t", quote=FALSE, row.names=FALSE)

# ── Sample positions ──
samples <- sort(unique(ibs$Sample))
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
  Sample = unique(ibs$Sample),
  Chromosome = chrom_levels
) %>%
  left_join(sample_positions, by="Sample") %>%
  left_join(chrinfo %>% select(chrom_name, len),
            by = c("Chromosome" = "chrom_name")) %>%
  mutate(Start = 0, End = len,
         Chromosome = factor(Chromosome, levels=chrom_levels))

ibs_plot <- ibs_plot %>%
  mutate(Chromosome = factor(Chromosome, levels=chrom_levels))

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
    low = "white", high = "#8B0000",
    name = "Mean Score", limits = c(0,100),
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
    breaks = seq(0, max(chrinfo$len, na.rm=TRUE)/1e6, by=50),
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
    strip.text = element_text(size = 48, face="bold"),
    axis.text.y = element_text(size = 14, face="bold", color="black"),
    axis.text.y.right = element_text(size = 48, face = "bold", color = "black", hjust = 0),
    # axis.text.x = element_text(size = 8, face="bold", angle=90, hjust=1),
    axis.text.x = element_text(size = 24, face="bold", hjust=1),
    axis.title.x = element_text(size = 48, face="bold", margin=margin(t=10)),
    axis.ticks = element_line(color="black"),
    axis.line = element_line(color="black"),
    panel.spacing.x = unit(0.6,"lines"),
    panel.spacing.y = unit(0.6,"lines"),
    panel.border = element_rect(color="black", fill=NA, linewidth=2.4),
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
    linewidth = 0.8  # Same thickness as panel border
  )
}

# ── Auto dimensions ──
n_samples <- length(unique(sample_positions$Sample))
n_chr <- length(chrom_levels)

plot_height <- max(2.3, (n_samples * 0.2) + 2)
plot_width  <- max(6, (n_chr * 4) + (n_samples * 0.2))  # add width scaling with samples

ggsave(opt$output, plot = p, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
