# calculate_snr.R
# --------------------------------------------------
# Calculate SNR from expression matrix and metadata
# --------------------------------------------------

#' Calculate SNR from input files
#'
#' @param expr_file path to expression matrix file (rows = features, cols = samples)
#' @param meta_file path to metadata file containing `sample` and `library`
#' @param sep field separator
#' @param verbose logical
#'
#' @return data.frame with number of features and SNR
#'
#'

snrdb_function<-function(exprMat,group){
  
  library(data.table)
  
  IDs<- colnames(exprMat)
  IDs.group.mat<-data.table(
    IDs=IDs,
    group=group) 
  
  pca_prcomp <- prcomp(t(exprMat),retx=T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  dt.perc.pcs <- data.table(PCX=1:nrow(pcs),
                            Percent=summary(pca_prcomp)$importance[2,],
                            AccumPercent=summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs,each=length(IDs)),
                        ID.B = rep(IDs,time=length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A,IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B,IDs.group.mat$IDs)]$group
  
  dt.dist[,Type:=ifelse(ID.A==ID.B,'Same',
                        ifelse(group.A==group.B,'Intra','Inter'))]
  
  dt.dist[,Dist:=(dt.perc.pcs[1]$Percent*(pcs[ID.A,1]-pcs[ID.B,1])^2+dt.perc.pcs[2]$Percent*(pcs[ID.A,2]-pcs[ID.B,2])^2)]
  
  dt.dist.stats <- dt.dist[,.(Avg.Dist=mean(Dist)),by=.(Type)]
  setkey(dt.dist.stats,Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist/dt.dist.stats['Intra']$Avg.Dist  
  
  signoise_db <- 10*log10(signoise)
  return(signoise_db)
  
}


calculate_snr_from_file <- function(expr_file,
                                    meta_file,
                                    snrdb_function,
                                    sep = "\t",
                                    verbose = TRUE) {
  
  ## -------------------------------
  ## 1. Read input files
  ## -------------------------------
  df <- read.table(
    expr_file,
    header = TRUE,
    row.names = 1,
    sep = sep,
    check.names = FALSE
  )
  
  meta_batch <- read.table(
    meta_file,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE
  )
  
  required_cols <- c("sample", "library")
  if (!all(required_cols %in% colnames(meta_batch))) {
    stop("meta_file must contain columns: sample, library")
  }
  
  ## -------------------------------
  ## 2. Filtering
  ## -------------------------------
  df <- df[!apply(df, 1, function(x) any(is.nan(x))), , drop = FALSE]
  df <- na.omit(df)
  
  n_features <- nrow(df)
  
  ## -------------------------------
  ## 3. Z-score normalization
  ## -------------------------------
  dat_z <- t(apply(df, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  
  ## -------------------------------
  ## 4. Sample annotation
  ## -------------------------------
  sample <- meta_batch$sample[
    match(colnames(dat_z), meta_batch$library)
  ]
  
  ## -------------------------------
  ## 5. SNR calculation
  ## -------------------------------
  snr_value <- round(
    snrdb_function(dat_z, as.factor(sample)),
    1
  )
  
  if (verbose) {
    message(
      "File: ", basename(expr_file),
      " | Features: ", n_features,
      " | SNR: ", snr_value
    )
  }
  
  return(data.frame(
    n_features = n_features,
    snr = snr_value
  ))
}
