#!/usr/bin/env Rscript
# Script: scoreHistogram.R
# Description: Plot IBS scores histogram that could help in defining score thresholds
# Author: c.s.sivsubramani@gmail.com
# Date: 2025-09-23

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggh4x)
  library(cowplot)
  library(scales)
})

# function to parse chromosome metadata
parse_chr_meta <- function(chr_meta_file) {
  chr_meta <- read.table(chr_meta_file, header = FALSE, sep = "\t") %>%
    setNames(c("chrom_name", "len", "chrom_num", "cum_len")) %>%
    mutate(
      chrom_num = as.numeric(chrom_num),
      chr_start = cum_len - len,
      chr_end = cum_len,
      mid = cum_len - (len / 2)
    )
  return(chr_meta)
}


option_list <- list(
  make_option(c("-c", "--chrinfo"), type="character",
              help="Chromosome metadata file"),
  make_option(c("-k", "--kcf"), type="character",
              help="KCF file"),
  make_option(c("-s", "--sample"), type="character",
              help="Sample name to plot"),
  make_option(c("-o", "--output"), type="character", default="score_histograms.pdf",
              help="Output PDF file (default: score_histograms.pdf)")
)


opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$chrinfo) || is.null(opt$kcf) || is.null(opt$sample)) {
  print_help(OptionParser(option_list=option_list))
  stop("Missing required arguments.", call.=FALSE)
}

# load and parse chromosome metadata
chr_meta <- parse_chr_meta(opt$chrinfo)


kcf_file <- opt$kcf
sample_name <- opt$sample
lines <- readLines(kcf_file)
header_line_index <- which(grepl("^#CHROM", lines))

# read data from the header line onward
kcf <- read.delim(kcf_file,
                  header = TRUE,
                  sep = "\t",
                  skip = header_line_index - 1,
                  comment.char = "")

# fix column name
colnames(kcf)[1] <- "CHROM"

if (!(sample_name %in% colnames(kcf))) {
  stop(paste("Sample", sample_name, "not found in KCF file."), call.=FALSE)
}

kcf_df <- kcf %>%
  select(CHROM, START, END, TOTAL_KMERS, !!sym(sample_name)) %>%
  rename(
    chrom_name = CHROM,
    start = START,
    end = END,
    total_kmers = TOTAL_KMERS
  ) %>%
  left_join(chr_meta %>% select(chrom_name, chrom_num), by="chrom_name") %>%
  filter(!is.na(chrom_num)) %>%
  separate(!!sym(sample_name), into = c("gt", "va", "ob", "id", "ld", "rd", "sc"), sep = ":") %>%
  mutate(across(c(va, ob, id, ld, rd, sc), as.numeric)) %>%
  mutate(dist = id + ld + rd) %>%
  select(chrom_name, chrom_num, start, end, total_kmers, va, ob, dist, sc)

kcf_df$chrom_num <- factor(kcf_df$chrom_num, levels = sort(unique(kcf_df$chrom_num)))

# Top plot (overall score distribution)
p1_clean <- ggplot(kcf_df, aes(x = sc)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "grey85", color = "black") +
  geom_density(color = "red", size = 0.8) +
  labs(
    title = NULL,
    x = "Score",
    y = "Density"
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    text = element_text(size = 14, face = "bold"),
    axis.title.y.right = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
  )

# Convert chrom_num to factor and filter
kcf_df$chrom_num <- as.factor(kcf_df$chrom_num)
filtered_df <- subset(kcf_df, sc >= 40 & sc <= 97)

# Define score boundaries for interpretation zones
score_boundaries <- c(60, 80, 97, 99)

# Bottom plot (per-chromosome density) with region dividers
p2 <- ggplot(filtered_df, aes(x = sc, color = chrom_num)) +
  geom_density(size = 1) +
  geom_point(
    data = filtered_df[0, ],
    aes(x = sc, y = 0, color = chrom_num),
    shape = 16, size = 3, show.legend = TRUE
  ) +
  # Add vertical dashed lines for interpretation boundaries
  geom_vline(xintercept = score_boundaries, linetype = "dashed", color = "grey40", size = 0.6) +

  scale_color_manual(
    name = "Chromosome",
    values = hue_pal()(length(levels(filtered_df$chrom_num))),
    guide = guide_legend(
      override.aes = list(shape = 16, size = 4, linetype = 0)
    )
  ) +
  scale_x_continuous(breaks = seq(60, 95, 5)) +
  labs(
    title = NULL,
    x = "Score",
    y = "Density"
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    text = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
  )

# Extract legend
legend <- get_legend(p2)
p2_clean <- p2 + theme(legend.position = "none")

# Combine plots
bottom_row <- plot_grid(p2_clean, legend, rel_widths = c(4, 1))

final_plot <- plot_grid(
  p1_clean, bottom_row,
  ncol = 1,
  rel_heights = c(1, 1),
  labels = c("a", "b"),
  label_size = 16,
  label_fontface = "bold"
)

# Save the final plot
ggsave(
  filename = opt$output,
  plot = final_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
# EOF
