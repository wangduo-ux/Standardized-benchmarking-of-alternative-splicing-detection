#!/usr/bin/env Rscript

# =========================
# Load libraries
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
})

# =========================
# Parse arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: script.R <tpm_dir> <annotation_file> ")
}


### tpm_dir contain the tpm data from all 21 isoform quantification pipelines for 19 high-quality laboratories, these data are available upon request (18801232285@163.com)
### annotation_file contain correspondence of each gene and their isoforms. This file is available upon request

tpm_dir <- args[1]
annotation_file <- args[2]


# =========================
# Load annotation (transcript → gene)
# =========================
anno <- read.table(annotation_file,
                   sep = "\t",
                   header = TRUE,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)

# =========================
# Define tools
# =========================
tools <- c(
  "assembly_star_stringtie","assembly_hisat2_stringtie","assembly_subjunc_stringtie",
  "star_cuffdiff_tpm","hisat2_cuffdiff_tpm",
  "star_stringtie","hisat2_stringtie","subjunc_stringtie",
  "star_featurecounts","hisat2_featurecounts","subjunc_featurecounts",
  "star_rsem","hisat2_rsem","bowtie2_rsem",
  "star_express","bowtie2_express",
  "star_salmon","bowtie2_salmon","salmon","kallisto","sailfish"
)

tool_labels <- c(
  "STA_Stri_A","HIS_Stri_A","Sub_Stri_A",
  "STA_Cuff","HIS_Cuff",
  "STA_Stri","HIS_Stri","Sub_Stri",
  "STA_feat","HIS_feat","Sub_feat",
  "STA_RS","HIS_RS","Bow_RS",
  "STA_eXpr","Bow_eXpr",
  "STA_Sal","Bow_Sal","Sal","Kall","Sail"
)

# =========================
# Parameters
# =========================
lab <- "lab1_"
time_points <- c("202207","202208","202209")

# =========================
# Function: compute gene-level CV
# =========================
compute_cv <- function(file_path) {
  
  if (!file.exists(file_path)) {
    warning(paste("Missing:", file_path))
    return(NULL)
  }
  
  df <- read.table(file_path, sep = "\t",
                   header = TRUE, check.names = FALSE)
  
  # Select relevant columns
  df <- df[, grep("transcript_id|202207|202208|202209",
                  colnames(df), value = TRUE)]
  
  # Normalize transcript IDs
  df$transcript_id <- sub("\\..*$", "", df$transcript_id)
  
  # Extract lab-specific columns
  sample_cols <- paste0(lab, time_points)
  
  if (!all(sample_cols %in% colnames(df))) {
    warning(paste("Missing columns in:", file_path))
    return(NULL)
  }
  
  df$mean_expr <- rowMeans(df[, sample_cols], na.rm = TRUE)
  
  # Merge with annotation
  df_merge <- anno %>%
    left_join(df[, c("transcript_id", "mean_expr")],
              by = "transcript_id")
  
  # Compute gene-level CV
  df_cv <- df_merge %>%
    group_by(gene_id) %>%
    summarise(
      mean_expr = mean(mean_expr, na.rm = TRUE),
      sd_expr   = sd(mean_expr, na.rm = TRUE),
      cv = ifelse(mean_expr > 0,
                  sd_expr / mean_expr,
                  NA)
    ) %>%
    select(gene_id, cv)
  
  return(df_cv)
}

# =========================
# Process all tools
# =========================
cv_list <- list()

for (i in seq_along(tools)) {
  
  tool <- tools[i]
  label <- tool_labels[i]
  
  file_path <- file.path(tpm_dir, paste0(tool, ".txt"))
  
  df_cv <- compute_cv(file_path)
  
  if (is.null(df_cv)) next
  
  colnames(df_cv)[2] <- label
  cv_list[[label]] <- df_cv
}

# =========================
# Merge all tools
# =========================
merged_df <- Reduce(function(x, y)
  merge(x, y, by = "gene_id", all = TRUE),
  cv_list
)

# =========================
# Transform to long format
# =========================
df_long <- merged_df %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Group",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# Log-transform CV (recommended)
df_long$Value <- log2(df_long$Value)

# =========================
# Color palette
# =========================
my_colors <- c(
  "#669999","#669999","#669999",
  "#EEB8C3","#EEB8C3",
  "#9ECAE1","#9ECAE1","#9ECAE1",
  "#CBB7FF","#CBB7FF","#CBB7FF",
  "#FFEBAC","#FFEBAC","#FFEBAC",
  "#FEB381","#FEB381",
  "#99DFB9","#99DFB9","#99DFB9",
  "#FB6501","#8A6800"
)

# =========================
# Plot
# =========================
p <- ggplot(df_long, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "log2(CV)", y = "Method") +
  scale_fill_manual(values = my_colors) +
  scale_y_discrete(limits = rev(unique(df_long$Group))) +
  coord_cartesian(xlim = c(-1, 5))

p