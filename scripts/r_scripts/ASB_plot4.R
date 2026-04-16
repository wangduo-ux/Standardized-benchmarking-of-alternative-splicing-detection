#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ===== 命令行参数 =====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_whippet_overlap.R input.csv output.pdf")
}

input_file <- args[1]
output_file <- args[2]

# ===== 读取数据 =====
dt <- fread(input_file)

# ===== 计算未交集 =====
no_overlap <- dt[, .(
  overlap_whippet_n = max(total_whippet_n) - sum(overlap_whippet_n)
), by = event_type]

no_overlap[, software := "No_overlap"]

# ===== 合并数据 =====
plot_dt <- rbind(
  dt[, .(event_type, software, overlap_whippet_n)],
  no_overlap
)

# ===== 排序软件顺序 =====
plot_dt[, software := factor(
  software,
  levels = c("No_overlap", setdiff(unique(software), "No_overlap"))
)]

# ===== 自适应图宽 =====
n_type <- uniqueN(plot_dt$event_type)
plot_width <- max(4, min(12, n_type * 1.2))

# ===== 配色 =====
cols <- scales::hue_pal()(length(unique(plot_dt$software)))
names(cols) <- unique(plot_dt$software)
cols["No_overlap"] <- "grey70"

# ===== 绘图 =====
p <- ggplot(plot_dt, aes(x = event_type, y = overlap_whippet_n, fill = software)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) +
  geom_text(aes(label = overlap_whippet_n), position = position_stack(vjust = 0.5), size = 3) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Event type",
    y = "Whippet count",
    fill = "Software"
  )

# ===== 保存图形 =====
ggsave(output_file, p, width = plot_width, height = 5, units = "in")
