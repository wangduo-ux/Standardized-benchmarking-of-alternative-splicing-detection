library(tidyr)
library(factoextra)

df <- "isoform-level or event-level expression matrix"

df <- log2(df+0.01)
df <- df[!apply(df, 1, function(x) length(unique(x))) == 1, ]
df <- na.omit(df)

data <- t(df)  #PCA分析需要将表达矩阵转置  
data.pca <- prcomp(data,retx=T,scale. = T)  #这是后续作图的文件 #输出特征向量

group=c(rep("M8",126),rep("F7",126),rep("D5",126),rep("T1",126),rep("T2",126),rep("D6",126))

#group=c(rep("M8",57),rep("F7",57),rep("D5",57),rep("T1",57),rep("T2",57),rep("D6",57)) ## for high-quality labs

color_vector <- c("#4AC5D9",'#7AC7A4', "#FFC55D", "#F36644",'#CC9749','#95B786')

p <- fviz_pca_ind(data.pca, 
                  fill.ind=group,
                  mean.point=F,  
                  label = "none", 
                  col.ind = "#404040",
                  geom.ind = "point",
                  pointshape = 21,
                  pointsize = 5,
                  palette = color_vector)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        panel.grid = element_blank(),axis.line = element_blank(),
        text = element_text(family = "Arial", size = 18)
  )

p
