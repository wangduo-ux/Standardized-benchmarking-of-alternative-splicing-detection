#!/usr/bin/env Rscript

#The scripts used for generated the Figure 3j and 5g

# =========================
# Load required libraries
# =========================
suppressPackageStartupMessages({
  library(ggplot2)
})

# =========================
# Input file
# =========================
input_file <- "input_data"  #avaiable in the Source Data file

# =========================
# Load data
# =========================
# Expected columns: Pipeline, Factor, P_value, η (effect size)
df <- read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE)

# =========================
# Define pipeline order
# =========================
software_levels <- c(
  "STA_Stri (A)_SUP","HIS_Stri (A)_SUP","Sub_Stri (A)_SUP",
  "STA_Cuff_SUP","HIS_Cuff_SUP",
  "STA_Stri_SUP","HIS_Stri_SUP","Sub_Stri_SUP",
  "STA_feat_SUP","HIS_feat_SUP","Sub_feat_SUP",
  "STA_RS_SUP","HIS_RS_SUP","Bow_RS_SUP",
  "STA_eXpr_SUP","Bow_eXpr_SUP",
  "STA_Sal_SUP","Bow_Sal_SUP","Sal_SUP","Kall_SUP","Sail_SUP"
)

df$Pipeline <- factor(df$Pipeline, levels = software_levels)

# =========================
# Define factor order
# NOTE: only one definition should be used
# =========================
df$Factor <- factor(
  df$Factor,
  levels = c("AS event number","Isoform number","GC content",
             "Mappability","PSI level","Isoform expression")
)

# =========================
# Transform P-values
# =========================
# Use -log10(P) for better interpretability (larger = more significant)
df$logP <- -log10(df$P_value)

# =========================
# Generate bubble plot
# =========================
ggplot(df, aes(x = Factor, y = Pipeline,
                    size = logP, color = η)) +
  geom_point(alpha = 0.85) +
  
  # Point size represents statistical significance
  scale_size_continuous(
    range = c(1, 10),
    name = "-log10(P-value)"
  ) +
  
  # Color represents effect size (η²)
  scale_color_gradient2(
    low = "#ED3640",
    mid = "#FFFFCB",
    high = "#296CA9",
    midpoint = 0.47,
    name = expression(eta^2)
  ) +
  
  # Reverse y-axis order for visual hierarchy
  scale_y_discrete(limits = rev(levels(df$Pipeline))) +
  
  # Minimal theme
  theme_minimal(base_size = 14) +
  labs(x = "Factor", y = "Pipeline") +
  
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
