#!/usr/bin/env Rscript

# =========================
# Load libraries
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(RColorBrewer)
  library(ggpubr)
})

# =========================
# Parse arguments
# =========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: script.R <base_dir> <truth_file> <output_prefix>")
}

base_dir <- args[1]
truth_file <- args[2]
output_prefix <- args[3]

# =========================
# Load truth data
# =========================
truth <- read.table(truth_file, header = TRUE, sep = ",",
                    stringsAsFactors = FALSE)

truth$ASE_compare <- paste0(truth$compare, "_", truth$ASE)

# =========================
# Function: compute PSI per condition
# =========================
compute_psi <- function(lab, tool_path) {
  
  result_df <- data.frame()
  conditions <- c("M8", "F7", "D6")
  
  for (cond in conditions) {
    
    file_path <- file.path(base_dir, tool_path,
                           paste0(lab, cond, ".psi"))
    
    if (!file.exists(file_path)) {
      warning(paste("Missing:", file_path))
      next
    }
    
    psi <- read.table(file_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
    
    # Remove rows with all NA
    psi <- psi[rowSums(!is.na(psi)) > 0, ]
    
    # Compute mean PSI per event
    psi$psi_mean <- rowMeans(psi, na.rm = TRUE)
    
    psi$ASE <- rownames(psi)
    psi$ASE_compare <- paste0(cond, psi$ASE)
    
    psi <- psi[, c("ASE_compare", "psi_mean")]
    psi <- psi[!is.na(psi$psi_mean), ]
    
    result_df <- rbind(result_df, psi)
  }
  
  return(result_df)
}

# =========================
# Define tools
# =========================

tool_labels <- c(
  "STA_Stri(A)_SUP","HIS_Stri(A)_SUP","Sub_Stri(A)_SUP",
  "STA_Cuff_SUP","HIS_Cuff_SUP",
  "STA_Stri_SUP","HIS_Stri_SUP","Sub_Stri_SUP",
  "STA_feat_SUP","HIS_feat_SUP","Sub_feat_SUP",
  "STA_RS_SUP","HIS_RS_SUP","Bow_RS_SUP",
  "STA_eXpr_SUP","Bow_eXpr_SUP",
  "STA_Sal_SUP","Bow_Sal_SUP","Sal_SUP",
  "Kall_SUP","Sail_SUP"
)

# =========================
# Process all tools
# =========================
lab <- "lab1_"
all_data <- list()

for (i in seq_along(tool_labels)) {
  
  label <- tool_labels[i]
  
  df_tool <- compute_psi(lab, label)
  
  if (nrow(df_tool) == 0) next
  
  df_tool$tool <- label
  colnames(df_tool)[2] <- "PSI"
  
  all_data[[label]] <- df_tool
}

# Combine all tools
df_long <- bind_rows(all_data)

# =========================
# Add reference
# =========================
ref_df <- data.frame(
  ASE_compare = truth$ASE_compare,
  PSI = truth$mean_delta_psi_mean,
  tool = "Quartet reference"
)

df_long <- bind_rows(df_long, ref_df)

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
  "#FB6501","#8A6800","#CC3333"
)

# =========================
# Plot 1: Ridge plot
# =========================
p_ridge <- ggplot(df_long, aes(x = PSI, y = tool, fill = tool)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "ΔPSI", y = "Method") +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(limits = c(-0.1, 1.1)) +
  scale_y_discrete(limits = rev(unique(df_long$tool)))

# =========================
# Plot 2: Violin plot
# =========================
p_violin <- ggplot(df_long, aes(x = tool, y = PSI, fill = tool)) +
  geom_violin(alpha = 0.6, trim = TRUE, size = 0.3) +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# =========================
# Combine plots
# =========================
p_combined <- ggarrange(
  p_ridge, p_violin,
  ncol = 1,
  labels = c("a", "b")
)

# =========================
# Save outputs
# =========================
Cairo::Cairo(
  file = paste0(output_prefix, ".pdf"),
  type = "pdf",
  width = 12,
  height = 8
)
print(p_combined)
dev.off()

ggsave(paste0(output_prefix, ".png"),
       p_combined,
       width = 12, height = 8, dpi = 300)