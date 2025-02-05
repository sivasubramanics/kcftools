#!/usr/bin/env Rscript

# check if packages ggplot2, dplyr, ggforce are installed, if not install them and load them silently
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggforce", quietly = TRUE)) {
  install.packages("ggforce")
}


# user input and usage message plotIBS.R <reference.fna.fai> <ibs.summary.tsv> [var]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("Usage: plotIBS.R <reference.fna.fai> <ibs.summary.tsv> <chr_name> <out_file>")
  quit(status = 1)
}

fai <- args[1]
ibs <- args[2]
chr_name <- args[3]
out <- args[4]

# fai <- "/Users/selva001/projects/work/wp4/kcftools/gene_test/genome.fna.fai"
# ibs <- "/Users/selva001/projects/work/wp4/kcftools/gene_test/LK085_lsatv11.window50kb.ibs.summary.tsv"
# out <- "/Users/selva001/projects/work/wp4/kcftools/gene_test/LK085_lsatv11.window50kb.ibs.chr1.png"
# chr_name <- "NC_056628.2"
# chr_name <- "OX465082.1"

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggforce))

# Read fai file
fai <- read.table(fai, header = FALSE, sep = "\t", stringsAsFactors = FALSE)[, 1:2]
fai <- fai %>% 
  setNames(c("seqname", "seqlen")) %>%
  filter(seqname == chr_name) %>% 
  mutate(seqlen = seqlen / 1e6)

# Read ibs file
ibs <- read.table(ibs, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ibs <- ibs %>% 
  filter(Chromosome == chr_name) %>%
  mutate(Start = Start / 1e6, End = End / 1e6, Length = Length / 1e6)

# Convert MeanScore to numeric (replace commas with dots)
ibs$MeanScore <- as.numeric(gsub(",", ".", ibs$MeanScore))

# Get unique samples
samples <- unique(ibs$Sample)

# Add sample positions for stacking
ibs <- ibs %>%
  mutate(SamplePosition = as.numeric(factor(Sample)))

# Create a data frame for the full chromosome length rectangles
sample_backgrounds <- data.frame(
  Sample = samples,
  SamplePosition = as.numeric(factor(samples)),
  Start = 0,
  End = fai$seqlen
)

bandwidth <- 0.3
# Create the plot
p <- ggplot() +
  # Add chromosome-wide rectangles for each sample
  geom_rect(
    data = sample_backgrounds,
    aes(
      xmin = Start,
      xmax = End,
      ymin = SamplePosition - bandwidth,
      ymax = SamplePosition + bandwidth
    ),
    fill = "grey90",
    color = NA
  ) +
  
  # Add IBS regions for each sample with fill based on MeanScore
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
  
  # Customize color scale for MeanScore
  # scale_fill_gradient(low = "yellow", high = "red", name = "Mean Score") +
  scale_fill_gradient(
    low = "yellow",
    high = "red",
    name = "Mean Score",
    limits = c(0, 100),
    oob = scales::squish
  ) +
  
  # Set axis labels and scales
  scale_y_continuous(
    breaks = sample_backgrounds$SamplePosition,
    labels = sample_backgrounds$Sample,
    expand = c(0.02, 0.02)
  ) +
  
  # Axis and title customization
  labs(
    x = "Position (Mb)",
    y = "Samples",
    title = paste("IBS Regions on Chromosome", chr_name)
  ) +
  
  # Custom theme and appearance
  theme(
    # plot title in the middle
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black"),
    panel.background = element_rect(fill = "transparent", colour = NA), # Transparent background for the plot
    plot.background = element_rect(fill = "transparent", colour = NA), # Transparent overall background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "bold", color = "black"),  # X axis labels
    axis.text.y = element_text(size = 24, face = "bold", color = "black"),  # Y axis labels
    axis.title.x = element_text(size = 30, face = "bold", color = "black"),  # X axis title
    axis.title.y = element_text(size = 30, face = "bold", color = "black"),  # Y axis title
    axis.line.x = element_line(color = "black", linewidth = 0.5),  # X axis line color
    axis.line.y = element_line(color = "black", linewidth = 0.5),  # Y axis line color
    axis.ticks = element_line(color = "black", linewidth = 0.5),  # Axis ticks color and size
    axis.ticks.length = unit(0.2, "cm"),  # Axis tick length
    legend.position = "bottom",  # Place legend at the bottom
    legend.text = element_text(size = 12, vjust = 0.5),  # Legend text size and alignment
    legend.title = element_text(size = 16, face = "bold"),  # Legend title styling
    legend.key.size = unit(1, "cm"),  # Size of legend keys
    legend.background = element_rect(fill = NA),  # Transparent legend background
    legend.box.background = element_rect(fill = NA, color = "black"),  # Border around legend box
    plot.margin = margin(10, 10, 10, 10)  # Add some margin around the plot
  )

plot_height <- (2 * nrow(sample_backgrounds))+1
ggsave(out, p, height = plot_height, width = 24, dpi = 300)
 

