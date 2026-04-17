library(precrec)
library(data.table)
library(ggplot2)


input_dir <- "your data"


## input_dir contains all differential splicing event analysis results from six event analysis tools. These data are available upon request


labs <- c(1,4,5,8,11,13,16,17,21,22,23,24,30,31,34,36,38,45)
scores_list <- list()
labels_list <- list()
names_list  <- character()
for (lab in labs) {
  
  path <- paste0("/input_dir/SUPPA2/", lab, ".txt")
  
  if (!file.exists(path)) next
  
  df <- fread(path)
  
  scores_list[[length(scores_list)+1]] <- -log10(df$tool_sig)
  labels_list[[length(labels_list)+1]] <- df$truth_sig
  names_list <- c(names_list, paste0("Lab", lab))
}

mm <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = names_list
)

res <- evalmod(mm)

p1 <- suppressWarnings(
  autoplot(res, "PRC") +
    theme_bw() +
    xlab("Recall") +
    ylab("Precision") +
    ggtitle("SUPPA2 PR curves") +
    theme(legend.position = "none") +
    scale_color_manual(values = rep("#CC855E", length(names_list)))
)




#################################################

scores_list <- list()
labels_list <- list()
names_list  <- character()
for (lab in labs) {
  
  path <- paste0("/input_dir/rMATS", lab, ".txt")
  
  if (!file.exists(path)) next
  
  df <- fread(path)
  
  scores_list[[length(scores_list)+1]] <- -log10(pmax(df$tool_sig, 1e-300))
  labels_list[[length(labels_list)+1]] <- df$truth_sig
  names_list <- c(names_list, paste0("Lab", lab))
}

mm <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = names_list,
  dsids = seq_along(scores_list)
)

res <- evalmod(mm)

p2 <- autoplot(res, "PRC") +
  theme_bw() +
  xlab("Recall") +
  ylab("Precision") +
  ggtitle("rMATS PR curves") +
  theme(legend.position="none") +
  scale_color_manual(values = rep("#994F45", length(names_list)))


#################################################

scores_list <- list()
labels_list <- list()
names_list  <- character()
for (lab in labs) {
  
  path <- paste0("/input_dir//PSI-Sigma", lab, ".txt")
  
  if (!file.exists(path)) next
  
  df <- fread(path)
  
  scores_list[[length(scores_list)+1]] <- -log10(pmax(df$tool_sig, 1e-300))
  labels_list[[length(labels_list)+1]] <- df$truth_sig
  names_list <- c(names_list, paste0("Lab", lab))
}

mm <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = names_list,
  dsids = seq_along(scores_list)
)

res <- evalmod(mm)

p3 <- autoplot(res, "PRC") +
  theme_bw() +
  xlab("Recall") +
  ylab("Precision") +
  ggtitle("PSI-Sigma PR curves") +
  theme(legend.position = "none") +
  scale_color_manual(values = rep("#456E8C", length(names_list)))


#################################################

scores_list <- list()
labels_list <- list()
names_list  <- character()
for (lab in labs) {
  
  path <- paste0("/input_dir/MAJIQ", lab, ".txt")
  
  if (!file.exists(path)) next
  
  df <- fread(path)
  
  scores_list[[length(scores_list)+1]] <- df$tool_sig
  labels_list[[length(labels_list)+1]] <- df$truth_sig
  names_list <- c(names_list, paste0("Lab", lab))
}

mm <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = names_list,
  dsids = seq_along(scores_list)
)

res <- evalmod(mm)

p4 <- autoplot(res, "PRC") +
  theme_bw() +
  xlab("Recall") +
  ylab("Precision") +
  ggtitle("MAJIQ PR curves") +
  theme(legend.position = "none") +
  scale_color_manual(values = rep("#828C61", length(names_list)))


#################################################

scores_list <- list()
labels_list <- list()
names_list  <- character()
for (lab in labs) {
  
  path <- paste0("/input_dir/Spladder", lab, ".txt")
  
  if (!file.exists(path)) next
  
  df <- fread(path)
  
  scores_list[[length(scores_list)+1]] <- -log10(pmax(df$tool_sig, 1e-300))
  labels_list[[length(labels_list)+1]] <- df$truth_sig
  names_list <- c(names_list, paste0("Lab", lab))
}

mm <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = names_list,
  dsids = seq_along(scores_list)
)

res <- evalmod(mm)

p5 <- autoplot(res, "PRC") +
  theme_bw() +
  xlab("Recall") +
  ylab("Precision") +
  ggtitle("Spladder PR curves") +
  theme(legend.position = "none") +
  scale_color_manual(values = rep("#25828D", length(names_list)))


#################################################

scores_list <- list()
labels_list <- list()
names_list  <- character()
for (lab in labs) {
  
  path <- paste0("/input_dir/Whippet", lab, ".txt")
  
  if (!file.exists(path)) next
  
  df <- fread(path)
  
  scores_list[[length(scores_list)+1]] <- df$tool_sig
  labels_list[[length(labels_list)+1]] <- df$truth_sig
  names_list <- c(names_list, paste0("Lab", lab))
}

mm <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = names_list,
  dsids = seq_along(scores_list)
)

res <- evalmod(mm)

p6 <- autoplot(res, "PRC") +
  theme_bw() +
  xlab("Recall") +
  ylab("Precision") +
  ggtitle("Whippet PR curves") +
  theme(legend.position = "none") +
  scale_color_manual(values = rep("#EDC780", length(names_list)))


p <- ggarrange(p1, p2,p3,p4,p5,p6,
               ncol = 6, nrow = 1,    # 添加图标签         # 标签的 Y 位置（可调）
               font.label = list(size = 14, face = "bold")  # 字体样式
)

p  <- p  + theme(plot.margin = margin(0,0,0,0))
pp <- pp + theme(plot.margin = margin(0,0,0,0))
ppp <- ggarrange(pp, p,
                ncol =1, nrow = 2, 
                labels = c("c", "d"),     # 添加图标签
    
                font.label = list(size = 14, face = "bold")  # 字体样式
)

ggsave("/home/wangduo/my_plot.png",plot = ppp,bg = "transparent", width = 15, height = 10, dpi = 300)



