#!/usr/bin/env Rscript
# script: plot_introgressions.R
# description: this script can be used to plot IBS (Identity By State) data from a summary file generated using kcftools findIBS
# contact: c.s.sivsubramani@gmail.com
# date: 16-04-2025

# ---------- install and load packages ----------
suppressMessages({
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  pkgs <- c("optparse", "ggplot2", "dplyr", "ggforce", "ggh4x", "svglite")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install)
  lapply(pkgs, library, character.only = TRUE)
})

# ---------- parse cmdline options ----------
option_list <- list(
  make_option(c("-r", "--reference"), type = "character", help = "Reference .fai file", metavar = "file"),
  make_option(c("-i", "--ibs"), type = "character", help = "IBS summary file (from findIBS plugin)", metavar = "file"),
  make_option(c("-c", "--chromosomes"), type = "character", help = "Chromosome name(s) comma-separated or file with list", metavar = "chr1[,chr2,...|file.txt]"),
  make_option(c("-s", "--samples"), type = "character", default = NULL, help = "Optional sample list file", metavar = "file"),
  make_option(c("-o", "--out"), type = "character", help = "Output image file (png/pdf/svg)", metavar = "file"),
  make_option("--legend", action = "store_true", default = FALSE, help = "Enable legend in output")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$reference) || is.null(opt$ibs) || is.null(opt$chromosomes) || is.null(opt$out)) {
  stop("Missing required arguments. Use -h for help.")
}

# ---------- mandatory inputs ----------
fai_file <- opt$reference
ibs_file <- opt$ibs
out_file <- opt$out
legend_enabled <- opt$legend

chromosomes <- if (file.exists(opt$chromosomes)) {
  readLines(opt$chromosomes)
} else {
  unlist(strsplit(opt$chromosomes, ","))
}

sample_list <- if (!is.null(opt$samples)) readLines(opt$samples) else NULL

# ---------- process chromosomes and IBS windows ----------
fai <- read.table(fai_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1:2]
colnames(fai) <- c("seqname", "seqlen")
fai <- fai[fai$seqname %in% chromosomes, ]
fai$seqname <- factor(fai$seqname, levels = chromosomes)
fai$seqlen <- fai$seqlen / 1e6

ibs <- read.table(ibs_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(Chromosome %in% chromosomes) %>%
  mutate(
    Chromosome = factor(Chromosome, levels = chromosomes),
    Start = Start / 1e6,
    End = End / 1e6,
    Length = Length / 1e6,
    MeanScore = as.numeric(gsub(",", ".", MeanScore))
  )

samples <- if (!is.null(sample_list)) sample_list else sort(unique(ibs$Sample))

# if the sample is not mentioned in the IBS file, add an empty band
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

# ---------- plotting ----------
p <- ggplot() +
  geom_rect(
    data = sample_backgrounds,
    aes(xmin = Start, xmax = End, ymin = SamplePosition - bandwidth, ymax = SamplePosition + bandwidth),
    fill = "grey95", color = NA
  ) +
  geom_rect(
    data = ibs,
    aes(xmin = Start, xmax = End, ymin = SamplePosition - bandwidth, ymax = SamplePosition + bandwidth, fill = MeanScore),
    color = NA
  ) +
  # scale_fill_gradient(low = "yellow", high = "#3b3b3b", name = "Mean Score", limits = c(0, 100), oob = scales::squish) +
  scale_fill_gradient(low = "yellow", high = "darkred", name = "Mean Score", limits = c(0, 100), oob = scales::squish) +
  scale_y_continuous(
    breaks = 1:length(samples),
    labels = samples,
    expand = c(0.01, 0.01)
  ) +
  labs(x = "Position (Mb)", y = NULL) +
  ggh4x::facet_nested(. ~ Chromosome, scales = "free_x", space = "free_x", switch = "y") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),  # ❌ Remove grid lines
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 8, face = "bold", angle = 90, hjust = 1),
    axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    axis.line = element_line(color = "black"),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.1, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # ✅ Plot border
    legend.position = ifelse(legend_enabled, "bottom", "none"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ---------- export the plot ----------
n_samples <- length(samples)
plot_height <- max(4, 0.25 * n_samples)
plot_width <- max(6, sum(fai$seqlen) * 0.25)
ext <- tolower(tools::file_ext(out_file))

ggsave(
  filename = out_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  dpi = if (ext == "png") 300 else NA,
  device = if (ext %in% c("pdf", "svg", "png")) ext else stop("Unsupported format"),
  limitsize = FALSE
)
#EOF