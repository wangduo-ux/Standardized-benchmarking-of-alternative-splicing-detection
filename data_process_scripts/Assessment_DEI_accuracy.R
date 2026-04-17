library(edgeR)
library(dplyr)
library(parallel)
library(BiocParallel)
library(NOISeq)
library(limma)
library(DESeq2)
library(sleuth)
library(BiocParallel)
library(precrec)



#########  This script evaluates all 11 differential isoform expression analysis methods and 21 isoform quantification pipeline combinations.
#########  It requires the Quartet reference dataset and the RT-qPCR reference dataset as input, along with the raw counts files for all analysis pipeline combinations.
#########  These data are available upon request. (18801232285@163.com)



calculate_metrics <- function(df) {
  
  df$combined <- paste(df$Final, df$sig, sep = "_")
  TP <- sum(df$combined == "up-regulate_up-regulate" | df$combined == "down-regulate_down-regulate")
  FP1 <- sum(df$combined == "up-regulate_down-regulate" | df$combined == "down-regulate_up-regulate")
  FP2 <- sum(df$combined == "non-DEI_down-regulate" | df$combined == "non-DEI_up-regulate")
  TN <- sum(df$combined == "non-DEI_non-DEI")
  TPR <- TP / 1077
  FN <- 1077 - TP - FP1
  precision <- TP / (TP + FP1 + FP2)
  FPR <- (FP1 + FP2) / (FP1 + FP2 + TN)
  mcc1 <- (TP * TN - (FP1 + FP2) * FN) / sqrt((TP + (FP1 + FP2)) * (TP + FN) * (TN + (FP1 + FP2)) * (TN + FN))
  
  metrics <- list(TP = TP, FP1 = FP1, FP2 =FP2,TN=TN,FN=FN,TPR=TPR,FPR=FPR,mcc = mcc1,precision = precision)
  return(metrics)
}

calculate_metrics_binary <- function(df) {
  
  df$true_bin <- df$Final != "non-DEI"
  df$pred_bin <- df$sig != "non-DEI"
  
  TP <- sum(df$true_bin & df$pred_bin)
  FP <- sum(!df$true_bin & df$pred_bin)
  FN <- sum(df$true_bin & !df$pred_bin)
  TN <- sum(!df$true_bin & !df$pred_bin)
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  FPR <- FP/(FP+TN)
  
  mcc <- (TP*TN - FP*FN) /
    sqrt(prod(c(TP+FP, TP+FN, TN+FP, TN+FN)))
  
  list(
    precision = precision,
    recall = recall,
    mcc = mcc
  )
}
calculate_metrics_binary_pcr <- function(df) {
  
  df$true_bin <- df$DEG != "NonDEG"
  df$pred_bin <- df$sig != "non-DEI"
  
  TP <- sum(df$true_bin & df$pred_bin)
  FP <- sum(!df$true_bin & df$pred_bin)
  FN <- sum(df$true_bin & !df$pred_bin)
  TN <- sum(!df$true_bin & !df$pred_bin)
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  FPR <- FP/(FP+TN)
  
  mcc <- (TP*TN - FP*FN) /
    sqrt(prod(c(TP+FP, TP+FN, TN+FP, TN+FN)))
  
  list(
    precision = precision,
    recall = recall,
    mcc = mcc
  )
}
calculate_metrics_pcr <- function(df) {
  
  df$combined <- paste(df$DEG, df$sig, sep = "_")
  TP <- sum(df$combined == "up-regulate_up-regulate" | df$combined == "down-regulate_down-regulate")
  FP1 <- sum(df$combined == "up-regulate_down-regulate" | df$combined == "down-regulate_up-regulate")
  FP2 <- sum(df$combined == "NonDEG_down-regulate" | df$combined == "NonDEG_up-regulate")
  TN <- sum(df$combined == "NonDEG_non-DEI")
  TPR <- TP / 78
  FN <- 78 - TP - FP1
  precision <- TP / (TP + FP1 + FP2)
  FPR <- (FP1 + FP2) / (FP1 + FP2 + TN)
  mcc1 <- (TP * TN - (FP1 + FP2) * FN) / sqrt((TP + (FP1 + FP2)) * (TP + FN) * (TN + (FP1 + FP2)) * (TN + FN))
  
  metrics <- list(TP = TP, FP1 = FP1, FP2 =FP2,TN=TN,FN=FN,TPR=TPR,FPR=FPR,mcc = mcc1,precision = precision)
  return(metrics)
}

calculate_AUCPR <- function(df,tool) {
  
  #df$y_true <- ifelse(df$Final == "non-DEI", 0, 1)
  if (tool == "edgeR"){
    df$score <- -log10(df$FDR)
  }
  else if (tool == "limma"){
    df$score <- -log10(df$adj.P.Val)
  }
  else if (tool == "noiseq"){
    df$score <- df$prob
  }
  else if (tool == "deseq"){
    df$score <- -log10(df$padj)
  }
  else if (tool == "sleuth"){
    df$score <- -log10(df$qval)
  }
  mm <- mmdata(scores = df$score, labels = df$y_true)
  res <- evalmod(mm)
  auc_table <- auc(res)
  aucs <- auc_table$aucs
  return(aucs)
}

calculate_AUCPR_pcr  <- function(df,tool) {
  
  #df$y_true <- ifelse(df$DEG == "NonDEG", 0, 1)
  if (tool == "edgeR"){
    df$score <- -log10(df$FDR)
  }
  else if (tool == "limma"){
    df$score <- -log10(df$adj.P.Val)
  }
  else if (tool == "noiseq"){
    df$score <- df$prob
  }
  else if (tool == "deseq"){
    df$score <- -log10(df$padj)
  } 
  else if (tool == "sleuth"){
    df$score <- -log10(df$qval)
  }
  mm <- mmdata(scores = df$score, labels = df$y_true)
  res <- evalmod(mm)
  auc_table <- auc(res)
  aucs <- auc_table$aucs
  return(aucs)
}

#################################
##### truth #####################

truth <- read.table("Quartet_isoform_reference_datasets", sep=",", header = T)
truth <- truth[, c("isoform", "compare","log2FC","Final")]
truth2 <- read.table("qPCR_isoform_reference_datasets", sep="\t", header = T)
truth2 <- truth2[, c("transcript_id", "sample_pair","log2FC","DEG")]

truth$isoform_compare <- paste0(truth$compare,"_",truth$isoform)
truth2$isoform_compare <- paste0(truth2$sample_pair,"_",truth2$transcript_id)

truth$y_true <- ifelse(abs(as.numeric(as.character(truth$log2FC))) >= 1, 1, 0)
truth2$y_true <- ifelse(abs(as.numeric(as.character(truth2$log2FC))) >= 1, 1, 0)
##################################
##### output #####################
file_path <- "output.txt" 
file <- file(file_path, "w")


#######################################
########### edgeR #####################
#######################################

processData_edgeR <- function(rawdata,mode) {
  group <- factor(c(rep("test",3),rep("control",3)))
  lrt <- NULL
  rawdata <- na.omit(rawdata)
  y <- DGEList(counts = rawdata, genes = rownames(rawdata), group = group)
  #keep <- rowSums(y$counts) >= x
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- normLibSizes(y)
  
  if (mode=="qCML"){
    ############################
    ##      classic edgeR     ##
    y <- estimateDisp(y)
    et <- exactTest(y)
  }
  else if (mode=="glm"){
    ############################
    ##         glm edgeR      ##
    design <- model.matrix(~group)
    y <- estimateDisp(y, design)
    fit  <- glmFit(y, design)
    et <- glmLRT(fit)
  }
  else if (mode=="QL"){
    ############################
    ## quasi-likelihood edgeR ##
    design <- model.matrix(~group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    et <- glmQLFTest(fit)
  }
  
  
  DEG <- summary(de <- decideTests(et, adjust.method="BH", p.value = 0.05, lfc = 1))
  down <- DEG[1, ]
  NotDEG <- DEG[2, ]
  up <- DEG[3, ]
  lrt <- topTags(et, n = nrow(y$counts))
  lrt <- as.data.frame(lrt)
  lrt$genes <- sapply(strsplit(lrt$genes, split = "|", fixed = TRUE), "[", 1)
  lrt[which(lrt$logFC >= 1 & lrt$FDR < 0.05), 'sig'] <- 'up-regulate'
  lrt[which(lrt$logFC <= -1 & lrt$FDR < 0.05), 'sig'] <- 'down-regulate'
  lrt[which(abs(lrt$logFC) < 1 | lrt$FDR >= 0.05), 'sig'] <- 'non-DEI'
  result <- list(lrt = lrt, down = down, up = up, NotDEG = NotDEG)
  return(result)
}

#######################################
########### limma #####################
#######################################

processData_limma <- function(countData,mode) {
  group <- factor(c(rep("test",3),rep("control",3)))
  design <- model.matrix(~group)
  colnames(design) <- levels(group)
  countData <- na.omit(countData)
  dge <- DGEList(countData)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  if (mode=="limma-trend"){
    # limma-trend  num
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
    x <- topTable(fit, coef=ncol(design),n=Inf)
  }
  else if (mode=="voom"){
    # voom When the library sizes are quite variable between samples
    v <- voom(dge, design, plot=F)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    x <- topTable(fit, coef=ncol(design),n=Inf)
  }
  #print(nrow(res))
  up <- nrow(subset(x, adj.P.Val < 0.05 & logFC >= 1))
  down <- nrow(subset(x, adj.P.Val < 0.05 & logFC <= -1))
  row_names <- rownames(x)
  rownames(x) <- sapply(strsplit(row_names, split = "|", fixed = TRUE), "[", 1)
  x[which(x$logFC >= 1 & x$adj.P.Val < 0.05), 'sig'] <- 'up-regulate'
  x[which(x$logFC <= -1 & x$adj.P.Val < 0.05), 'sig'] <- 'down-regulate'
  x[which(abs(x$logFC) < 1 | x$adj.P.Val >= 0.05), 'sig'] <- 'non-DEI'
  NotDEG <- sum(x$sig == "non-DEI", na.rm = TRUE)
  result <- list(lrt = x, down = down, up = up, NotDEG = NotDEG)
  return(result)
}

########################################
########### MOISeq #####################
########################################


processData_noiseq <- function(countData,mode) {
  myfactors = data.frame(factor = c( "test", "test", "test","control", "control", "control"))
  mydata <- readData(data = countData,factors= myfactors)
  #print(nrow(res))
  if (mode=="noiseqbio"){
    mynoiseqbio = noiseqbio(mydata, k = 0.5, norm = "rpkm", factor = "factor",lc = 1, r = 50, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter = 1,cv.cutoff = 500, cpm = 1)
    mynoiseq.deg = degenes(mynoiseqbio, q = 0.8, M = NULL)
    x <- mynoiseqbio@results[[1]] 
    #x <- mynoiseq.deg
    up <- nrow(subset(x, prob > 0.8 & log2FC >= 1))
    down <- nrow(subset(x, prob > 0.8 & log2FC <= -1))
    x[which(x$log2FC >= 1 & x$prob > 0.8), 'sig'] <- 'up-regulate'
    x[which(x$log2FC <= -1 & x$prob > 0.8), 'sig'] <- 'down-regulate'
    x[which(abs(x$log2FC) < 1 | x$prob <= 0.8), 'sig'] <- 'non-DEI'
  }
  else if (mode=="noiseqreal"){
    mynoiseq = noiseq(mydata, norm = "rpkm", factor = "factor",  replicates = "technical")
    x <- mynoiseq@results[[1]] 
    up <- nrow(subset(x, prob > 0.8 & M >= 1))
    down <- nrow(subset(x, prob > 0.8 & M <= -1))
    x[which(x$M >= 1 & x$prob > 0.8), 'sig'] <- 'up-regulate'
    x[which(x$M <= -1 & x$prob > 0.8), 'sig'] <- 'down-regulate'
    x[which(abs(x$M) < 1 | x$prob <= 0.8), 'sig'] <- 'non-DEI'
  }
  
  NotDEG <- sum(x$sig == "non-DEI", na.rm = TRUE)
  result <- list(lrt = x, down = down, up = up, NotDEG = NotDEG)
  return(result)
}

########################################
############ DESeq #####################
########################################

processData_deseq <- function(countData) {
  condition <- factor(c(rep("test",3), rep("control",3)))
  countData <- na.omit(countData)
  countData[]  <- lapply(countData , round)
  coldata <- data.frame(row.names = colnames(countData), condition)
  dds <- DESeqDataSetFromMatrix(countData = countData,colData = coldata, design = ~ condition)
  #keep <- rowSums(counts(dds) >= 10) >= 3
  #dds <- dds[keep,]
  #value <- quantile(apply(cpm(counts(dds)),1,function(row) median(row, na.rm = TRUE)), probs = 0.7)
  #print(value)
  n <- nrow(dds)
  dds <- DESeq(dds)
  res <- results(dds) ### independent filtering 
  res <- res[!is.na(res$padj), ]
  #print(nrow(res))
  up <- nrow(subset(res, padj < 0.05 & log2FoldChange >= 1))
  down <- nrow(subset(res, padj < 0.05 & log2FoldChange <= -1))
  NotDEG <- nrow(res)-up-down
  row_names <- rownames(res)
  rownames(res) <- sapply(strsplit(row_names, split = "|", fixed = TRUE), "[", 1)
  res[which(res$log2FoldChange >= 1 & res$padj < 0.05), 'sig'] <- 'up-regulate'
  res[which(res$log2FoldChange <= -1 & res$padj < 0.05), 'sig'] <- 'down-regulate'
  res[which(abs(res$log2FoldChange) < 1 | res$padj >= 0.05), 'sig'] <- 'non-DEI'
  result <- list(lrt = res, down = down, up = up, NotDEG = NotDEG,n=n)
  return(result)
}


##############################################
################# sleuth #####################
##############################################

sleuth<-function(sample_id,path){
  
  s2c <- data.frame(
    sample = sample_id,
    condition   = rep(c("control", "case"), each = 3),
    path = paste0(path , "/",sample_id))
  s2c$condition <- factor(
    s2c$condition,
    levels = c("control", "case")   # control → 参考，case → 对比
  )
  so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  #so <- sleuth_lrt(so, 'reduced', 'full') 输出结果没有方向
  so <- sleuth_wt(so, which_beta = 'conditioncase', which_model = "full")
  sleuth_table <- sleuth_results(  so,  'conditioncase',test_type = 'wt')
  sleuth_table <- sleuth_table[!is.na(sleuth_table$qval), ]
  #sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE) 适用于lrt拟合
  sleuth_table[which(sleuth_table$b >= 1 & sleuth_table$qval < 0.05), 'sig'] <- 'up-regulate'
  sleuth_table[which(sleuth_table$b <= -1 & sleuth_table$qval < 0.05), 'sig'] <- 'down-regulate'
  sleuth_table[which(abs(sleuth_table$b) < 1 | sleuth_table$qval >= 0.05), 'sig'] <- 'non-DEI'
  up <- nrow(subset(sleuth_table, sig == "up-regulate"))
  down <- nrow(subset(sleuth_table, sig == "down-regulate"))
  NotDEG <- nrow(subset(sleuth_table, sig == "non-DEI"))
  result <- list(lrt = sleuth_table, down = down, up = up, NotDEG = NotDEG)
  return(result)
}


# "star_stringtie_trans","hisat2_stringtie_trans","subjunc_stringtie_trans","star_featurecounts","hisat2_featurecounts","subjunc_featurecounts","star_rsem","hisat2_rsem","bowtie2_rsem","star_express","bowtie2_express","star_salmon","bowtie_salmon","salmon","kallisto","sailfish"
param <- MulticoreParam(workers=30)
for (tool in c("star_stringtie_trans","hisat2_stringtie_trans","subjunc_stringtie_trans","star_featurecounts","hisat2_featurecounts","subjunc_featurecounts","star_rsem","hisat2_rsem","bowtie2_rsem","star_express","bowtie2_express","star_salmon","bowtie_salmon","salmon","kallisto","sailfish")){ # c("bowtie_salmon","salmon","kallisto","sailfish")){
  if (tool %in% c("star_stringtie_trans","hisat2_stringtie_trans","subjunc_stringtie_trans")) {
    counts <- read.table(paste0("input_dir",tool,".txt"),row.names=1,sep = '\t', header = TRUE,check.names=FALSE)
  }
  else {
    counts <- read.table(paste0("input_dir",tool,".txt"),row.names=1,sep = '\t', header = TRUE,check.names=FALSE)
  }
  
  for (i in  c(1,4,5,8,11,13,16,17,21,22,23,24,30,31,34,36,38,45)) {
    lab <- paste("lab", i, "_",sep = "")
    
    cols_to_select <- grep(lab, names(counts), value = TRUE)
    if (length(cols_to_select) > 0) {
      raw_data <- counts[, cols_to_select]
      M8_D6 <- raw_data[,c(paste0(lab,"202207"),paste0(lab,"202208"),paste0(lab,"202209"),paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"))]
      F7_D6 <- raw_data[,c(paste0(lab,"202210"),paste0(lab,"202211"),paste0(lab,"202212"),paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"))]
      D5_D6 <- raw_data[,c(paste0(lab,"202213"),paste0(lab,"202214"),paste0(lab,"202215"),paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"))]
      
      counts_list <- list(
        M8_D6 = M8_D6,
        F7_D6 = F7_D6,
        D5_D6 = D5_D6
      )
      
      
      ####################################
      ##############  edgeR ##############    
      for (mode in c("qCML","glm","QL")){ 
        cat(tool,lab,"edgeR",mode)
        # 并行执行 processData
        qlf_results <- bplapply(names(counts_list), function(sample_name) {
          counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
          result <- processData_edgeR(counts_matrix, mode)  # 运行 edgeR 差异分析
          return(result)  # 返回结果
        }, BPPARAM = param)
        
        # 给结果命名，方便后续访问
        names(qlf_results) <- names(counts_list)
        
        lrt1 <- qlf_results$M8_D6$lrt
        lrt2 <- qlf_results$F7_D6$lrt
        lrt3 <- qlf_results$D5_D6$lrt
        lrt1$genes <- sub("\\..*$", "", lrt1$genes)
        lrt2$genes <- sub("\\..*$", "", lrt2$genes)
        lrt3$genes <- sub("\\..*$", "", lrt3$genes)
        lrt1$genes <- paste0("M8/D6_", lrt1$genes)
        lrt2$genes <- paste0("F7/D6_", lrt2$genes)
        lrt3$genes <- paste0("D5/D6_", lrt3$genes)
        lrt <- bind_rows(lrt1, lrt2, lrt3)
        #write.csv(lrt,paste0("/cold_data/wangduo/isoform_benchmark/DEG_file/",tool,"_edgeR_",mode,"_",lab,".txt"),index=F)
        df1 <- merge(truth,lrt,by.x="isoform_compare",by.y = "genes")
        df2 <- merge(truth2,lrt,by.x="isoform_compare",by.y = "genes")
        res1 <- calculate_metrics(df1)
        res2 <- calculate_metrics_pcr(df2)
        calculate_metrics_binary(df1)
        calculate_metrics_binary_pcr(df2)
        auc1 <- calculate_AUCPR(df1,"edgeR")
        auc2 <- calculate_AUCPR_pcr(df2,"edgeR")
        cat(res1$mcc,res2$mcc,auc1,auc2)
        writeLines(paste(tool,mode,lab,nrow(df1),
                         qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                         qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                         qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                         "quartet:",res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      }
      
      
      ####################################
      ##############  limma ##############    
      for (mode in c("limma-trend","voom")){ 
        cat(tool,lab,"limma",mode)
        # 并行执行 processData
        qlf_results <- bplapply(names(counts_list), function(sample_name) {
          counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
          result <- processData_limma(counts_matrix, mode)  # 运行 edgeR 差异分析
          return(result)  # 返回结果
        }, BPPARAM = param)
        
        # 给结果命名，方便后续访问
        names(qlf_results) <- names(counts_list)
        
        lrt1 <- as.data.frame(qlf_results$M8_D6$lrt)
        lrt2 <- as.data.frame(qlf_results$F7_D6$lrt)
        lrt3 <- as.data.frame(qlf_results$D5_D6$lrt)
        rownames(lrt1) <- sub("\\..*$", "", rownames(lrt1))
        rownames(lrt2) <- sub("\\..*$", "", rownames(lrt2))
        rownames(lrt3) <- sub("\\..*$", "", rownames(lrt3))
        rownames(lrt1) <- paste0("M8/D6_", rownames(lrt1),sep = "")
        rownames(lrt2) <- paste0("F7/D6_", rownames(lrt2),sep = "")
        rownames(lrt3) <- paste0("D5/D6_", rownames(lrt3),sep = "")
        lrt <- bind_rows(as.data.frame(lrt1), as.data.frame(lrt2), as.data.frame(lrt3))
        #write.csv(lrt,paste0("/cold_data/wangduo/isoform_benchmark/DEG_file/",tool,"_limma_",mode,"_",lab,".txt"),index=F)
        df1 <- merge(truth,lrt,by.x="isoform_compare" ,by.y= "row.names")
        df2 <- merge(truth2,lrt,by.x = "isoform_compare",by.y="row.names")
        res1 <- calculate_metrics(df1)
        res2 <- calculate_metrics_pcr(df2)
        auc1 <- calculate_AUCPR(df1,"limma")
        auc2 <- calculate_AUCPR_pcr(df2,"limma")
        cat(res1$mcc,res2$mcc,auc1,auc2)
        writeLines(paste(tool,mode,lab,nrow(df1),
                         qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                         qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                         qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                         "quartet:",res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      }
      
      ####################################
      ###### DESeq #######################
      
      cat(tool,lab,"DESeq")
      # 并行执行 processData
      qlf_results <- bplapply(names(counts_list), function(sample_name) {
        counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
        result <- processData_deseq(counts_matrix)  # 运行 edgeR 差异分析
        return(result)  # 返回结果
      }, BPPARAM = param)
      
      # 给结果命名，方便后续访问
      names(qlf_results) <- names(counts_list)
      
      
      lrt1 <- as.data.frame(qlf_results$M8_D6$lrt)
      lrt2 <- as.data.frame(qlf_results$F7_D6$lrt)
      lrt3 <- as.data.frame(qlf_results$D5_D6$lrt)
      rownames(lrt1) <- sub("\\..*$", "", rownames(lrt1))
      rownames(lrt2) <- sub("\\..*$", "", rownames(lrt2))
      rownames(lrt3) <- sub("\\..*$", "", rownames(lrt3))
      rownames(lrt1) <- paste0("M8/D6_", rownames(lrt1))
      rownames(lrt2) <- paste0("F7/D6_", rownames(lrt2))
      rownames(lrt3) <- paste0("D5/D6_", rownames(lrt3))
      lrt <- bind_rows(as.data.frame(lrt1), as.data.frame(lrt2), as.data.frame(lrt3))
      #write.csv(lrt,paste0("/cold_data/wangduo/isoform_benchmark/DEG_file/",tool,"_DEseq_",lab,".txt"),index=F)
      df1 <- merge(truth,lrt,by.x = "isoform_compare",by.y="row.names")
      df2 <- merge(truth2,lrt,by.x = "isoform_compare",by.y="row.names")
      
      res1 <- calculate_metrics(df1)
      res2 <- calculate_metrics_pcr(df2)
      auc1 <- calculate_AUCPR(df1,"deseq")
      auc2 <- calculate_AUCPR_pcr(df2,"deseq")
      cat(res1$mcc,res2$mcc,auc1,auc2)
      writeLines(paste(tool,"_",lab,nrow(df1),
                       qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                       qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                       qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                       "quartet:", res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      
      
      #####################################
      ##############  NOIseq ############## 
      M8_D6 <- raw_data[,c(paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"),paste0(lab,"202207"),paste0(lab,"202208"),paste0(lab,"202209"))]
      F7_D6 <- raw_data[,c(paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"),paste0(lab,"202210"),paste0(lab,"202211"),paste0(lab,"202212"))]
      D5_D6 <- raw_data[,c(paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"),paste0(lab,"202213"),paste0(lab,"202214"),paste0(lab,"202215"))]
      
      counts_list <- list(
        M8_D6 = M8_D6,
        F7_D6 = F7_D6,
        D5_D6 = D5_D6
      )
      
      for (mode in c("noiseqbio","noiseqreal")){ 
        cat(tool,lab,"NOISeq",mode)
        # 并行执行 processData
        qlf_results <- bplapply(names(counts_list), function(sample_name) {
          counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
          result <- processData_noiseq(counts_matrix, mode)  # 运行 edgeR 差异分析
          return(result)  # 返回结果
        }, BPPARAM = param)
        
        # 给结果命名，方便后续访问
        names(qlf_results) <- names(counts_list)
        lrt1 <- as.data.frame(qlf_results$M8_D6$lrt)
        lrt2 <- as.data.frame(qlf_results$F7_D6$lrt)
        lrt3 <- as.data.frame(qlf_results$D5_D6$lrt)
        rownames(lrt1) <- sub("\\..*$", "", rownames(lrt1))
        rownames(lrt2) <- sub("\\..*$", "", rownames(lrt2))
        rownames(lrt3) <- sub("\\..*$", "", rownames(lrt3))
        rownames(lrt1) <- paste0("M8/D6_", rownames(lrt1),sep = "")
        rownames(lrt2) <- paste0("F7/D6_", rownames(lrt2),sep = "")
        rownames(lrt3) <- paste0("D5/D6_", rownames(lrt3),sep = "")
        lrt <- bind_rows(as.data.frame(lrt1), as.data.frame(lrt2), as.data.frame(lrt3))
        lrt <- lrt[rowSums(is.na(lrt)) != ncol(lrt), ]
        #write.csv(lrt,paste0("/cold_data/wangduo/isoform_benchmark/DEG_file/",tool,"_NOIseq_",mode,"_",lab,".txt"),index=F)
        df1 <- merge(truth,lrt,by.x="isoform_compare" ,by.y= "row.names")
        df2 <- merge(truth2,lrt,by.x = "isoform_compare",by.y="row.names")
        res1 <- calculate_metrics(df1)
        res2 <- calculate_metrics_pcr(df2)
        auc1 <- calculate_AUCPR(df1,"noiseq")
        auc2 <- calculate_AUCPR_pcr(df2,"noiseq")
        cat(res1$mcc,res2$mcc,auc1,auc2)
        writeLines(paste(tool,mode,lab,nrow(df1),
                         qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                         qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                         qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                         "quartet:",res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      }
    }
  }
}



for (tool in c("assembly_star_stringtie","assembly_hisat2_stringtie","assembly_subjunc_stringtie")) { #c("assembly_star_stringtie","assembly_hisat2_stringtie","assembly_subjunc_stringtie")){
  
  for (i in  1:48) {
    lab <- paste("lab", i,sep = "")
    path <- paste0("input_dir",tool,"_",lab,".txt")
    if (file.exists(path)) {
      raw_data <- read.table(path,row.names=1,sep = ',', header = TRUE,check.names=FALSE)
      M8_D6 <- raw_data[,c(paste0(lab,"_202207"),paste0(lab,"_202208"),paste0(lab,"_202209"),paste0(lab,"_202222"),paste0(lab,"_202223"),paste0(lab,"_202224"))]
      F7_D6 <- raw_data[,c(paste0(lab,"_202210"),paste0(lab,"_202211"),paste0(lab,"_202212"),paste0(lab,"_202222"),paste0(lab,"_202223"),paste0(lab,"_202224"))]
      D5_D6 <- raw_data[,c(paste0(lab,"_202213"),paste0(lab,"_202214"),paste0(lab,"_202215"),paste0(lab,"_202222"),paste0(lab,"_202223"),paste0(lab,"_202224"))]
      
      counts_list <- list(
        M8_D6 = M8_D6,
        F7_D6 = F7_D6,
        D5_D6 = D5_D6
      )
      
      ####################################
      ##############  edgeR ##############    
      for (mode in c("qCML","glm","QL")){ 
        cat(tool,lab,"edgeR",mode)
        # 并行执行 processData
        qlf_results <- bplapply(names(counts_list), function(sample_name) {
          counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
          result <- processData_edgeR(counts_matrix, mode)  # 运行 edgeR 差异分析
          return(result)  # 返回结果
        }, BPPARAM = param)
        
        # 给结果命名，方便后续访问
        names(qlf_results) <- names(counts_list)
        
        lrt1 <- qlf_results$M8_D6$lrt
        lrt2 <- qlf_results$F7_D6$lrt
        lrt3 <- qlf_results$D5_D6$lrt
        lrt1$genes <- paste0("M8/D6_", lrt1$genes)
        lrt2$genes <- paste0("F7/D6_", lrt2$genes)
        lrt3$genes <- paste0("D5/D6_", lrt3$genes)
        lrt <- bind_rows(lrt1, lrt2, lrt3)
        df1 <- merge(truth,lrt,by.x="isoform_compare",by.y = "genes")
        df2 <- merge(truth2,lrt,by.x="isoform_compare",by.y = "genes")
        res1 <- calculate_metrics(df1)
        res2 <- calculate_metrics_pcr(df2)
        auc1 <- calculate_AUCPR(df1,"edgeR")
        auc2 <- calculate_AUCPR_pcr(df2,"edgeR")
        cat(res1$mcc,res2$mcc)
        writeLines(paste(tool,mode,lab,nrow(df1),
                         qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                         qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                         qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                         "quartet:",res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      }
      
      
      ####################################
      ##############  limma ##############    
      for (mode in c("limma-trend","voom")){ 
        cat(tool,lab,"limma",mode)
        # 并行执行 processData
        qlf_results <- bplapply(names(counts_list), function(sample_name) {
          counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
          result <- processData_limma(counts_matrix, mode)  # 运行 edgeR 差异分析
          return(result)  # 返回结果
        }, BPPARAM = param)
        
        # 给结果命名，方便后续访问
        names(qlf_results) <- names(counts_list)
        
        lrt1 <- as.data.frame(qlf_results$M8_D6$lrt)
        lrt2 <- as.data.frame(qlf_results$F7_D6$lrt)
        lrt3 <- as.data.frame(qlf_results$D5_D6$lrt)
        
        rownames(lrt1) <- paste0("M8/D6_", rownames(lrt1),sep = "")
        rownames(lrt2) <- paste0("F7/D6_", rownames(lrt2),sep = "")
        rownames(lrt3) <- paste0("D5/D6_", rownames(lrt3),sep = "")
        lrt <- bind_rows(as.data.frame(lrt1), as.data.frame(lrt2), as.data.frame(lrt3))
        
        df1 <- merge(truth,lrt,by.x="isoform_compare" ,by.y= "row.names")
        df2 <- merge(truth2,lrt,by.x = "isoform_compare",by.y="row.names")
        res1 <- calculate_metrics(df1)
        res2 <- calculate_metrics_pcr(df2)
        auc1 <- calculate_AUCPR(df1,"limma")
        auc2 <- calculate_AUCPR_pcr(df2,"limma")
        cat(res1$mcc,res2$mcc)
        writeLines(paste(tool,mode,lab,nrow(df1),
                         qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                         qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                         qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                         "quartet:",res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      }
      
      ####################################
      ###### DESeq #######################
      
      cat(tool,lab,"DESeq")
      # 并行执行 processData
      qlf_results <- bplapply(names(counts_list), function(sample_name) {
        counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
        result <- processData_deseq(counts_matrix)  # 运行 edgeR 差异分析
        return(result)  # 返回结果
      }, BPPARAM = param)
      
      # 给结果命名，方便后续访问
      names(qlf_results) <- names(counts_list)
      
      
      lrt1 <- as.data.frame(qlf_results$M8_D6$lrt)
      lrt2 <- as.data.frame(qlf_results$F7_D6$lrt)
      lrt3 <- as.data.frame(qlf_results$D5_D6$lrt)
      
      rownames(lrt1) <- paste0("M8/D6_", rownames(lrt1))
      rownames(lrt2) <- paste0("F7/D6_", rownames(lrt2))
      rownames(lrt3) <- paste0("D5/D6_", rownames(lrt3))
      lrt <- bind_rows(as.data.frame(lrt1), as.data.frame(lrt2), as.data.frame(lrt3))
      df1 <- merge(truth,lrt,by.x = "isoform_compare",by.y="row.names")
      df2 <- merge(truth2,lrt,by.x = "isoform_compare",by.y="row.names")
      
      res1 <- calculate_metrics(df1)
      res2 <- calculate_metrics_pcr(df2)
      auc1 <- calculate_AUCPR(df1,"deseq")
      auc2 <- calculate_AUCPR_pcr(df2,"deseq")
      cat(res1$mcc,res2$mcc)
      writeLines(paste(tool,"_",lab,nrow(df1),
                       qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                       qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                       qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                       "quartet:", res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      
      
      #####################################
      ##############  NOIseq ############## 
      M8_D6 <- raw_data[,c(paste0(lab,"_202222"),paste0(lab,"_202223"),paste0(lab,"_202224"),paste0(lab,"_202207"),paste0(lab,"_202208"),paste0(lab,"_202209"))]
      F7_D6 <- raw_data[,c(paste0(lab,"_202222"),paste0(lab,"_202223"),paste0(lab,"_202224"),paste0(lab,"_202210"),paste0(lab,"_202211"),paste0(lab,"_202212"))]
      D5_D6 <- raw_data[,c(paste0(lab,"_202222"),paste0(lab,"_202223"),paste0(lab,"_202224"),paste0(lab,"_202213"),paste0(lab,"_202214"),paste0(lab,"_202215"))]
      
      counts_list <- list(
        M8_D6 = M8_D6,
        F7_D6 = F7_D6,
        D5_D6 = D5_D6
      )
      
      for (mode in c("noiseqbio","noiseqreal")){ 
        cat(tool,lab,"NOISeq",mode)
        # 并行执行 processData
        qlf_results <- bplapply(names(counts_list), function(sample_name) {
          counts_matrix <- counts_list[[sample_name]]  # 获取对应的 counts 矩阵
          result <- processData_noiseq(counts_matrix, mode)  # 运行 edgeR 差异分析
          return(result)  # 返回结果
        }, BPPARAM = param)
        
        # 给结果命名，方便后续访问
        names(qlf_results) <- names(counts_list)
        lrt1 <- as.data.frame(qlf_results$M8_D6$lrt)
        lrt2 <- as.data.frame(qlf_results$F7_D6$lrt)
        lrt3 <- as.data.frame(qlf_results$D5_D6$lrt)
        
        rownames(lrt1) <- paste0("M8/D6_", rownames(lrt1),sep = "")
        rownames(lrt2) <- paste0("F7/D6_", rownames(lrt2),sep = "")
        rownames(lrt3) <- paste0("D5/D6_", rownames(lrt3),sep = "")
        lrt <- bind_rows(as.data.frame(lrt1), as.data.frame(lrt2), as.data.frame(lrt3))
        lrt <- lrt[rowSums(is.na(lrt)) != ncol(lrt), ]
        df1 <- merge(truth,lrt,by.x="isoform_compare" ,by.y= "row.names")
        df2 <- merge(truth2,lrt,by.x = "isoform_compare",by.y="row.names")
        res1 <- calculate_metrics(df1)
        res2 <- calculate_metrics_pcr(df2)
        auc1 <- calculate_AUCPR(df1,"noiseq")
        auc2 <- calculate_AUCPR_pcr(df2,"noiseq")
        cat(res1$mcc,res2$mcc)
        writeLines(paste(tool,mode,lab,nrow(df1),
                         qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                         qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                         qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                         "quartet:",res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
      }
    }
  }
}


for (tool in c("salmon/ensembl_salmon","salmon/star-salmon","salmon/bowtie2-salmon","kallisto/ensembl_kallisto","sailfish/sailfish")) { #c("salmon/ensembl_salmon","salmon/star-salmon","salmon/bowtie2-salmon","kallisto/ensembl_kallisto","sailfish/sailfish")){
  path = paste0("input_dir",tool)
  for (i in c(1,4,5,8,11,13,16,17,21,22,23,24,30,31,34,36,38,45)){
    lab <- paste0("lab",i,"_")
    
    sample_list <- list(
      M8_D6 = c(paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"),paste0(lab,"202207"),paste0(lab,"202208"),paste0(lab,"202209")),
      F7_D6 = c(paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"),paste0(lab,"202210"),paste0(lab,"202211"),paste0(lab,"202212")),
      D5_D6 = c(paste0(lab,"202222"),paste0(lab,"202223"),paste0(lab,"202224"),paste0(lab,"202213"),paste0(lab,"202214"),paste0(lab,"202215"))
    )
    
    # 并行执行 processData
    qlf_results <- bplapply(names(sample_list), function(sample_name) {
      sample_id <- sample_list[[sample_name]]  # 获取对应的 counts 矩阵
      result <- sleuth(sample_id, path)  
      return(result)  # 返回结果
    }, BPPARAM = param)
    
    # 给结果命名，方便后续访问
    names(qlf_results) <- names(sample_list)
    lrt1 <- qlf_results$M8_D6$lrt
    lrt2 <- qlf_results$F7_D6$lrt
    lrt3 <- qlf_results$D5_D6$lrt
    lrt1$target_id  <- sub("\\..*$", "", lrt1$target_id )
    lrt2$target_id  <- sub("\\..*$", "", lrt2$target_id )
    lrt3$target_id  <- sub("\\..*$", "", lrt3$target_id )
    lrt1$target_id  <- paste0("M8/D6_", lrt1$target_id )
    lrt2$target_id  <- paste0("F7/D6_", lrt2$target_id )
    lrt3$target_id  <- paste0("D5/D6_", lrt3$target_id )
    lrt <- bind_rows(lrt1, lrt2, lrt3)
    df1 <- merge(truth,lrt,by.x="isoform_compare",by.y = "target_id")
    df2 <- merge(truth2,lrt,by.x="isoform_compare",by.y = "target_id")
    res1 <- calculate_metrics(df1)
    res2 <- calculate_metrics_pcr(df2)
    
    auc1 <- calculate_AUCPR(df1,"sleuth")
    auc2 <- calculate_AUCPR_pcr(df2,"sleuth")
    writeLines(paste(tool,lab,nrow(df1),
                     qlf_results$M8_D6$down,qlf_results$M8_D6$up,qlf_results$M8_D6$NotDEG,
                     qlf_results$F7_D6$down,qlf_results$F7_D6$up,qlf_results$F7_D6$NotDEG,
                     qlf_results$D5_D6$down,qlf_results$D5_D6$up,qlf_results$D5_D6$NotDEG,
                     "quartet:",
                     res1$TPR,res1$precision,res1$FPR,res1$mcc,auc1[1],auc1[2],"pcr:",res2$TPR,res2$precision,res2$FPR,res2$mcc,auc2[1],auc2[2]), con = file)
    
  }
}
close(file)