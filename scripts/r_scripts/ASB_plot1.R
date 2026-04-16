#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexUpset)
  library(patchwork)
})

# 读取命令行参数
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

if (is.na(input_file) || is.na(output_file)) {
  stop("Usage: Rscript plot_upset_auto.R input.csv output.pdf")
}

# 读取数据
df <- fread(input_file)

# 获取事件类型和工具列表
event_types <- unique(df$event_type)
tools <- setdiff(colnames(df), c("event_id", "event_type"))

n_events <- length(event_types)
n_tools <- length(tools)

########################
# 自动布局
########################
if (n_events <= 4) {
  ncol <- 1
} else {
  ncol <- 2
}
nrow <- ceiling(n_events / ncol)

########################
# 子图尺寸
########################
subplot_width <- 11
subplot_height <- 6
final_width <- ncol * subplot_width
final_height <- nrow * subplot_height

########################
# 字体自动调整
########################
base_size <- 11
if (n_tools >= 5) {
  base_size <- 9
}

########################
# 生成 UpSet 图
########################
plots <- lapply(event_types, function(ev) {
  sub <- df[event_type == ev]
  
  n_intersections <- min(20, nrow(sub))
  
  ylim_val <- max(
    colSums(
      sub[, which(sapply(sub[, -1, with = FALSE], is.numeric)) + 1, with = FALSE],
      na.rm = TRUE
    )
  ) * 1.5
  
  upset(
    sub,
    intersect = tools,
    width_ratio = 0.1,
    set_sizes = (
      upset_set_size() + 
        geom_text(
          aes(label = after_stat(count)), 
          hjust = 1.1, 
          size = 3, 
          stat = 'count'
        ) + 
        expand_limits(y = ylim_val) + 
        theme(axis.text.x = element_text(angle = 90, size = 8))
    ),
    min_size = 1,
    n_intersections = n_intersections,
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts = TRUE,
        text = list(
          size = 3,
          angle = 90,
          vjust = 0.5,
          hjust = 0
        )
      )
    )
  ) +
    ggtitle(paste0(ev, " (n=", nrow(sub), ")")) +
    theme(
      text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2, face = "bold")
    )
})

########################
# 拼图
########################
final_plot <- wrap_plots(plots, ncol = ncol)

########################
# 输出
########################
ggsave(
  filename = output_file,
  plot = final_plot,
  width = final_width,
  height = final_height
)