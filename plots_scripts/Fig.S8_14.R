#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# =========================
# 参数
# =========================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: script.R <junction_dir> <annotation_file>")
}

junction_dir <- args[1]  # for example ./data/junction/
annotation_file <- args[2]  # for example ./data/coverage/

### junction_dir contain the TP, FN, FP of junction detection, these files are avaiable upon request
### annotations such as reads coverage, coverage uniformity, length of asseciated isoforms, exon count distribution of isoforms, GC content, mappability distribution. these files are avaiable upon request


labs <- c(1,4,5,8,11,13,16,17,21,22,23,24,29,30,31,34,36,38,45)

anno <- read.table(anno_file, header=TRUE, sep="\t", stringsAsFactors = TRUE)

plot_list <- list()
summary_list <- list()

for (i in labs) {
  cat("Processing lab", i, "\n")
  
  # read junction file
  junction_file <- file.path(junction_dir, paste0("lab", i, "_annotated.txt"))
  if (!file.exists(junction_file)) next
  df <- read.table(junction_file, header=TRUE, sep="\t", stringsAsFactors=TRUE)
  df$sample <- str_extract(df$ID_2, "^[^_]+")
  df$ID <- sub("^[^_]+_", "", as.character(df$ID_2))
  
  # read annotation file
  
  dff <- read.table(annotation_file, header=TRUE, sep="\t", stringsAsFactors=TRUE)

  
  dff <- merge(df, dff, by="transcript_id")
  
  # 绘图
  axis_y_title <- if (i %% 5 == 1) element_text() else element_blank()
  axis_y_text  <- if (i %% 5 == 1) element_text() else element_blank()
  axis_x_text  <- if (i > 14) element_text() else element_blank()
  
  p <- ggplot(dff, aes(x=classification, y=mean, fill=classification)) +
    geom_violin(trim=FALSE, scale="width", size=0.3) +
    geom_boxplot(width=0.1, outlier.shape=NA, color="black", size=0.2) +
    stat_summary(fun=mean, geom="crossbar", width=0.1, color="white", size=0.2) +
    scale_fill_manual(values=c("TP"="#1b9e77","FN"="#d95f02","FP"="#7570b3")) +
    coord_cartesian(ylim=c(0,1)) +
    theme_minimal(base_family="Arial") +
    theme(
      legend.position="none",
      axis.title.x=element_blank(),
      axis.title.y=axis_y_title,
      axis.text.x=axis_x_text,
      axis.text.y=axis_y_text,
      text=element_text(color="black")
    )
  
  plot_list[[length(plot_list)+1]] <- p
}

