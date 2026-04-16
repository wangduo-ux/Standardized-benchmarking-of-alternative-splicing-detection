#!/usr/bin/env Rscript
.libPaths("/home/wangduo/anaconda3/envs/r451/lib/R/library")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)

input_file  <- args[1]
output_file <- args[2]

if (is.na(input_file) || is.na(output_file)) {
  stop("Usage: Rscript plot_dse_panel_bar.R input.csv output.pdf")
}

df <- read.delim(input_file, check.names = FALSE)

event_types <- unique(df$Event_type)

n_events <- length(event_types)

########################
# 自动布局
########################

if (n_events <= 2) {
  ncol <- 1
} else if (n_events <= 4) {
  ncol <- 2
} else {
  ncol <- 3
}

nrow <- ceiling(n_events / ncol)

########################
# subplot尺寸
########################

subplot_width  <- 5
subplot_height <- 5

final_width  <- ncol * subplot_width
final_height <- nrow * subplot_height

########################
# 生成图
########################

plots <- lapply(event_types, function(ev){
  
  sub <- df[df$Event_type == ev, ]
  
  p1 <- ggplot(sub, aes(x=Software, y=DSE, fill=Software)) +
    geom_bar(stat="identity") +
    geom_text(
      aes(label=DSE),
      vjust=-0.3,
      size=3
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      legend.position="none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    labs(
      title = ev,
      y = "DSE"
    )
  
  p2 <- ggplot(sub, aes(x=Software, y=nonDSE, fill=Software)) +
    geom_bar(stat="identity") +
    geom_text(
      aes(label=nonDSE),
      vjust=-0.3,
      size=3
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      legend.position="none",
      axis.text.x = element_text(angle=45, hjust=1)
    ) +
    labs(
      y = "nonDSE",
      x = "Software"
    )
  
  p1 / p2 + plot_layout(heights=c(1,1))
  
})

########################
# 拼图
########################

final_plot <- wrap_plots(
  plots,
  ncol = ncol
)

########################
# 输出
########################

ggsave(
  filename = output_file,
  plot = final_plot,
  width = final_width,
  height = final_height
)