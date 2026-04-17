#!/usr/bin/env Rscript

## ============================================================
## Junction truth set construction from long-read RNA-seq
## ============================================================

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
})

## ===============================
## CONFIGURATION
## ===============================

SJ_PATH <- "/cold_data/wangduo/lr_data/minimap2/"
FC_PATH <- "/cold_data/wangduo/lr_data/minimap2/featurecounts/"
GTF_FILE <- "/hot_warm_data/wangduo/reference/Homo_sapiens.GRCh38.109.gtf.gz"
KNOWN_SJ_FILE <- "/hot_warm_data/wangduo/transcriptome/STAR/genomeIndex/sjdbList.fromGTF.out.tab"
GENE_INFO_FILE <- "/hot_warm_data/wangduo/transcriptome/STAR/genomeIndex/geneInfo.tab"

LABS <- c(
  "M_PAB_LN_B1","P_ONT_LG_B1","P_ONT_LN_B2",
  "D_ONT_LG_B1","M_PAB_LG_B2","P_ONT_LG_B2","D_ONT_LW_B1"
)

SAMPLES <- c("M8","F7","D5","D6")

MIN_READS_PER_REP <- 1
MIN_REP_SUPPORT   <- 2
MIN_LAB_SUPPORT   <- 5
MIN_GENE_LABS     <- 5

## ===============================
## FUNCTIONS
## ===============================

read_junction_replicates <- function(path, lab, sample) {
  
  dfs <- lapply(1:3, function(i) {
    df <- read.table(
      paste0(path, lab, "_", sample, "_", i, ".SJ.out.motif.SJ.tab"),
      sep = "\t", header = FALSE
    )
    colnames(df) <- c("chr","start","end","r1","r2","r3","r4","motif")
    df$ID <- paste(df$chr, df$start, df$end, sep = "_")
    df[, c("ID", "r2")]
  })
  
  merged <- Reduce(function(x,y) merge(x,y,by="ID",all=TRUE), dfs)
  colnames(merged) <- c("ID","rep1","rep2","rep3")
  
  merged$status <- ifelse(
    rowSums(merged[,2:4] >= MIN_READS_PER_REP, na.rm = TRUE) >= MIN_REP_SUPPORT,
    "D","U"
  )
  
  merged
}

gene_coverage_filter <- function(labs, sample) {
  
  gene_hits <- c()
  
  for (lab in labs) {
    dfs <- lapply(1:3, function(i) {
      df <- read.table(
        paste0(FC_PATH, lab, "_", sample, "_", i, ".txt"),
        header = TRUE, sep = "\t", check.names = FALSE
      )
      df[, c("Geneid", ncol(df))]
    })
    
    merged <- Reduce(function(x,y) merge(x,y,by="Geneid",all=TRUE), dfs)
    keep <- rowSums(merged[,2:4] >= 2, na.rm = TRUE) >= 2
    gene_hits <- c(gene_hits, merged$Geneid[keep])
  }
  
  as.data.frame(table(gene_hits))
}

build_truth_single_sample <- function(sample) {
  
  message("Processing sample: ", sample)
  
  dfs <- lapply(LABS, function(lab) {
    df <- read_junction_replicates(SJ_PATH, lab, sample)
    colnames(df)[2:5] <- paste0(lab, c("_1","_2","_3","_status"))
    df
  })
  
  merged <- Reduce(function(x,y) merge(x,y,by="ID",all=TRUE), dfs)
  
  status_cols <- grep("_status$", colnames(merged), value = TRUE)
  
  merged$status <- case_when(
    rowSums(merged[,status_cols]=="D", na.rm=TRUE) >= MIN_LAB_SUPPORT ~ "Y_H",
    rowSums(is.na(merged[,status_cols])) == length(status_cols) ~ "N_H",
    TRUE ~ "not_determined"
  )
  
  ## annotate known junctions
  known <- read.table(KNOWN_SJ_FILE, sep="\t", header=FALSE)
  colnames(known) <- c("seqnames","start","end","strand","gene")
  known$ID <- paste(known$seqnames, known$start, known$end, sep="_")
  known$anno <- "known"
  
  merged <- merge(merged, known, by="ID", all=TRUE)
  merged$sample <- sample
  
  ## gene annotation
  merged <- merged %>% separate_rows(gene, sep=",")
  
  gene_info <- read.table(GENE_INFO_FILE, sep="\t", skip=1)
  colnames(gene_info) <- c("gene_id","gene_name","type")
  gene_info$index <- row.names(gene_info)
  
  merged <- left_join(merged, gene_info, by=c("gene"="index"))
  
  ## gene coverage QC
  cov <- gene_coverage_filter(LABS, sample)
  low_cov_genes <- setdiff(
    merged$gene_id[merged$status=="N_H"],
    cov$gene_hits[cov$Freq >= MIN_GENE_LABS]
  )
  
  merged <- merged %>% filter(!gene_id %in% low_cov_genes)
  merged
}

## ===============================
## MAIN
## ===============================

all_truth <- lapply(SAMPLES, build_truth_single_sample)
all_truth <- do.call(rbind, all_truth)

truth_known <- all_truth %>%
  filter(!is.na(strand), status != "not_determined")

truth_novel <- all_truth %>%
  filter(is.na(strand), status == "Y_H")

write.table(
  truth_known,
  "junction_annotated_truth_c1.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

write.table(
  truth_novel,
  "junction_novel_truth_c1.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

message("Junction truth construction finished.")
