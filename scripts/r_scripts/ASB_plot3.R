#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ===== 命令行参数 =====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript ASB_plot_pairwise_correct.R input.csv output.pdf")
}

input_path <- args[1]
output_path <- args[2]

# ===== 读取数据 =====
dt <- fread(input_path)

required_cols <- c("value1", "value2", "tool1", "tool2", "type")
if (!all(required_cols %in% colnames(dt))) {
  stop("Missing required columns")
}

# ===== 主题 =====
theme_clean <- function() {
  theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(size = 0.3),
      axis.ticks = element_line(size = 0.3)
    )
}

# ===== 核心函数（完全 pairwise）=====
make_matrix_plot <- function(sub_dt, title_name) {
  tools <- sort(unique(c(sub_dt$tool1, sub_dt$tool2)))
  k <- length(tools)
  if (k < 2) return(NULL)
  
  plist <- list()
  
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      t1 <- tools[j]
      t2 <- tools[i]
      
      # ===== pairwise 抽取（关键）=====
      df <- sub_dt[(tool1 == t1 & tool2 == t2) | (tool1 == t2 & tool2 == t1)]
      if (nrow(df) > 0) {
        df[, x := ifelse(tool1 == t1, value1, value2)]
        df[, y := ifelse(tool1 == t1, value2, value1)]
        df <- df[is.finite(x) & is.finite(y)]
      }
      
      show_x <- (i == k)
      show_y <- (j == 1)
      xlab <- if (show_x) t1 else NULL
      ylab <- if (show_y) t2 else NULL
      
      # ===== 对角线 =====
      if (i == j) {
        df_diag <- rbind(
          sub_dt[tool1 == t1, .(v = value1)],
          sub_dt[tool2 == t1, .(v = value2)]
        )
        df_diag <- df_diag[is.finite(v)]
        if (nrow(df_diag) == 0) {
          p <- ggplot() + theme_void()
        } else {
          p <- ggplot(df_diag, aes(v)) +
            geom_histogram(bins = 40, fill = "#3B6EA8") +
            scale_x_continuous(breaks = pretty_breaks(3))
        }
        
        # ===== 左下 scatter =====
      } else if (i > j) {
        if (nrow(df) == 0) {
          p <- ggplot() + theme_void()
        } else {
          p <- ggplot(df, aes(x, y)) +
            geom_point(size = 0.4, alpha = 0.5, color = "#3B6EA8") +
            scale_x_continuous(breaks = pretty_breaks(3)) +
            scale_y_continuous(breaks = pretty_breaks(3))
        }
        
        # ===== 右上 r/n =====
      } else {
        n <- nrow(df)
        r <- if (n > 2) cor(df$x, df$y) else NA
        text_size <- 2
        label <- sprintf("r = %.2f\nn = %d", r, n)
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = label, size = text_size, lineheight = 0.9) +
          xlim(0, 1) + ylim(0, 1) + theme_void()
        plist[[length(plist) + 1]] <- p
        next
      }
      
      p <- p + labs(x = xlab, y = ylab) + theme_clean() +
        theme(
          axis.title.x = if (show_x) element_text(size = 5) else element_blank(),
          axis.title.y = if (show_y) element_text(size = 5) else element_blank(),
          axis.text.x = if (show_x) element_text(size = 4) else element_blank(),
          axis.text.y = if (show_y) element_text(size = 4) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.ticks.y = if (show_y) element_line() else element_blank()
        )
      
      plist[[length(plist) + 1]] <- p
    }
  }
  
  matrix_plot <- wrap_plots(plist, ncol = k)
  title_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = title_name, size = 2.5, fontface = "bold") +
    theme_void()
  
  return(title_plot / matrix_plot + plot_layout(heights = c(0.05, 1)))
}

# ===== 主流程 =====
plot_list <- list()
for (tp in unique(dt$type)) {
  sub_dt <- dt[type == tp]
  p <- make_matrix_plot(sub_dt, tp)
  if (!is.null(p)) {
    plot_list[[tp]] <- p
  }
}

plot_list <- Filter(Negate(is.null), plot_list)
if (length(plot_list) == 0) stop("No valid plots")

# ===== 排版 =====
n_plot <- length(plot_list)
ncol <- if (n_plot < 2) 1 else 2
nrow <- ceiling(n_plot / ncol)

max_k <- max(sapply(plot_list, function(p) { length(p$patches$plots)^(1/2) }))
panel_size <- if (max_k <= 4) 2.2 else if (max_k <= 6) 2.6 else 3

final_plot <- wrap_plots(plot_list, ncol = ncol)

ggsave(
  output_path, final_plot,
  width = panel_size * max_k * ncol,
  height = panel_size * max_k * nrow * 1.05,
  device = cairo_pdf
)

message("Done: ", output_path)