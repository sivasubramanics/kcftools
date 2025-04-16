#!/usr/bin/env Rscript

suppressMessages({
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  pkgs <- c("optparse", "ggplot2", "dplyr", "ggforce", "ggh4x")
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  lapply(pkgs, library, character.only = TRUE)
})

# ---------- Option Parsing ----------
option_list <- list(
  make_option(c("-r", "--reference"), type="character", help="Reference .fai file", metavar="file"),
  make_option(c("-i", "--ibs"), type="character", help="IBS summary file", metavar="file"),
  make_option(c("-c", "--chromosomes"), type="character", help="Chromosome name(s) comma-separated or file with list", metavar="chr1[,chr2,...|file.txt]"),
  make_option(c("-s", "--samples"), type="character", default=NULL, help="Optional sample list file", metavar="file"),
  make_option(c("-o", "--out"), type="character", help="Output image file (png/pdf/svg)", metavar="file"),
  make_option("--legend", action="store_true", default=FALSE, help="Enable legend in output")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$reference) || is.null(opt$ibs) || is.null(opt$chromosomes) || is.null(opt$out)) {
  stop("Missing required arguments. Use -h for help.")
}

# ---------- Input Parsing ----------
fai_file <- opt$reference
ibs_file <- opt$ibs
out_file <- opt$out
legend_enabled <- opt$legend

# Chromosome parsing: file or comma list
if (file.exists(opt$chromosomes)) {
  chromosomes <- readLines(opt$chromosomes)
} else {
  chromosomes <- unlist(strsplit(opt$chromosomes, ","))
}

# Sample list: optional
if (!is.null(opt$samples)) {
  sample_list <- readLines(opt$samples)
} else {
  sample_list <- NULL
}

# ---------- Read Data ----------
fai <- read.table(fai_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1:2]
colnames(fai) <- c("seqname", "seqlen")

# Filter chromosomes & preserve order
fai <- fai[fai$seqname %in% chromosomes, ]
fai$seqname <- factor(fai$seqname, levels = chromosomes)
fai$seqlen <- fai$seqlen / 1e6  # Convert to Mb

ibs <- read.table(ibs_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(Chromosome %in% chromosomes) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chromosomes),
    Start = Start / 1e6,
    End = End / 1e6,
    Length = Length / 1e6,
    MeanScore = as.numeric(gsub(",", ".", MeanScore))
  )

# Samples
if (!is.null(sample_list)) {
  samples <- sample_list
} else {
  samples <- sort(unique(ibs$Sample))
}

# Create SamplePosition and background
sample_backgrounds <- expand.grid(Sample = samples, Chromosome = chromosomes, stringsAsFactors = FALSE) %>%
  left_join(fai, by = c("Chromosome" = "seqname")) %>%
  mutate(
    SamplePosition = as.numeric(factor(Sample, levels = samples)),
    Start = 0,
    End = seqlen
  )

ibs <- ibs %>%
  mutate(SamplePosition = as.numeric(factor(Sample, levels = samples)))

bandwidth <- 0.25

# ---------- Plot ----------
chrom_widths <- fai %>%
  arrange(factor(seqname, levels = chromosomes)) %>%
  pull(seqlen)

p <- ggplot() +
  geom_rect(
    data = sample_backgrounds,
    aes(
      xmin = Start,
      xmax = End,
      ymin = SamplePosition - bandwidth,
      ymax = SamplePosition + bandwidth
    ),
    fill = "grey95", color = NA
  ) +
  geom_rect(
    data = ibs,
    aes(
      xmin = Start,
      xmax = End,
      ymin = SamplePosition - bandwidth,
      ymax = SamplePosition + bandwidth,
      fill = MeanScore
    ),
    color = NA
  ) +
  scale_fill_gradient(
    low = "yellow",
    high = "red",
    name = "Mean Score",
    limits = c(0, 100),
    oob = scales::squish
  ) +
  scale_y_continuous(
    breaks = 1:length(samples),
    labels = samples,
    expand = c(0.01, 0.01)
  ) +
  labs(
    x = "Position (Mb)",
    y = NULL
  ) +
  ggh4x::facet_grid2(
    . ~ Chromosome,
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 1),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Add plot border
    panel.spacing.x = unit(0.2, "lines"),
    panel.spacing.y = unit(0.1, "lines"),
    legend.position = ifelse(legend_enabled, "bottom", "none"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ---------- Save Image ----------
n_samples <- length(samples)
n_chroms <- length(chromosomes)
file_ext <- tolower(tools::file_ext(out_file))

# Dynamic size
per_chr_width <- 2.5
per_sample_height <- 0.25
plot_width <- max(6, per_chr_width * sum(chrom_widths) / mean(chrom_widths))  # Adjust width by total chromosome span
plot_height <- max(4, per_sample_height * n_samples)

message("Saving plot (", toupper(file_ext), ") ", plot_width, " x ", plot_height, " inches")

if (file_ext == "png") {
  ggsave(out_file, p, width = plot_width, height = plot_height, dpi = 600, limitsize = FALSE)
} else if (file_ext == "svg") {
  ggsave(out_file, p, width = plot_width, height = plot_height, device = "svg", limitsize = FALSE)
} else if (file_ext == "pdf") {
  ggsave(out_file, p, width = plot_width, height = plot_height, device = "pdf", limitsize = FALSE)
} else {
  stop("Unsupported output format: use .png, .svg or .pdf")
}
