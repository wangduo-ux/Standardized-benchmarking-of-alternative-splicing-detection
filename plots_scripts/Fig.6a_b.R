#!/usr/bin/env Rscript

# =========================
# Load libraries
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexUpset)
})

# =========================
# Parse arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: script.R <input_file>")
}

input_file <- args[1]

### input_file contains the intersection of six event tools, this file is avaiable upon request (18801232285@163.com)

# =========================
# Load data
# =========================
# Expected columns:
# set (e.g., "SE_SUPPA2"), event_id, event_type, ID
df <- read.table(input_file,
                 header = TRUE,
                 sep = "\t",
                 check.names = FALSE,
                 stringsAsFactors = FALSE)



# =========================
# Ensure required columns exist
# =========================
required_tools <- c("SUPPA2", "rMATS", "PSI-Sigma", "MAJIQ", "Spladder", "Whippet")

missing_tools <- setdiff(required_tools, colnames(df_binary))
if (length(missing_tools) > 0) {
  stop(paste("Missing tools in input:", paste(missing_tools, collapse = ", ")))
}

# =========================
# Event type colors
# =========================
event_colors <- c(
  "SE" = "#1f77b4",
  "A3SS" = "#ff7f0e",
  "A5SS" = "#2ca02c",
  "AF" = "#d62728",
  "AL" = "#9467bd",
  "MX" = "#8c564b",
  "RI" = "#e377c2"
)

# =========================
# Generate UpSet plot
# =========================
p <- upset(
  df_binary,
  intersect = required_tools,
  n_intersections = 40,
  min_size = 5,
  width_ratio = 0.1,
  annotations = list(
    "Event type" = (
      ggplot(mapping = aes(fill = event_type)) +
        geom_bar() +
        scale_fill_manual(values = event_colors)
    )
  )
)

p