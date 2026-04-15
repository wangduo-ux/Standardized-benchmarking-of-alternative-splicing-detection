library(dplyr)
library(reshape2)
library(ggridges)
library(ggplot2)

lab <- "lab1_"
dir <- "~/tpm"

df_list <- list()
n <- 1

for (tool in c("star_stringtie_assemebly","hisat2_stringtie_assemebly","subjunc_stringtie_assemebly","star_cuffdiff","hisat2_cuffdiff","star_stringtie","hisat2_stringtie","subjunc_stringtie","star_featurecounts","hisat2_featurecounts","subjunc_featurecounts","star_rsem","hisat2_rsem","bowtie2_rsem","star_express","bowtie2_express","star_salmon","bowtie2_salmon","salmon","kallisto","sailfish")){
  tpm <- read.table(paste0(dir,tool,".tpm"),sep = '\t', header = TRUE,check.names=FALSE)
  
  tpm <- tpm[, grep("transcript_id|202207|202208|202209", colnames(tpm), value = TRUE)]
  tpm$transcript_id <- sub("\\..*$", "", tpm$transcript_id)
  cols_to_select <- c("transcript_id", grep(lab, names(tpm), value = TRUE))
  if (length(cols_to_select) > 0) {
      
      raw_data <- tpm[, cols_to_select]
      raw_data$mean <- rowMeans(raw_data[, c(paste0(lab,"202207"), paste0(lab,"202208"), paste0(lab,"202209"))], na.rm = TRUE)
      df_merge <- df %>%
        left_join(raw_data, by = "transcript_id")   
      
      df_cv <- df_merge %>%
        group_by(gene_id) %>%
        summarise(
          mean_expr = mean(mean, na.rm = TRUE),       
          sd_expr   = sd(mean, na.rm = TRUE),       
          cv   = sd_expr / mean_expr            
        )
      df_cv <- data.frame(df_cv)
      df_cv <- df_cv[,c("gene_id","cv")]
      df_list[[n]] <- df_cv
    }
  }

merged_df <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), df_list)
colnames(merged_df) <- c("gene_id","STA_Stri (A)",	"HIS_Stri (A)",	"Sub_Stri (A)","STA_Cuff",	"HIS_Cuff",	"STA_Stri",	"HIS_Stri",	"Sub_Stri",	"STA_feat",	"HIS_feat",	"Sub_feat",	"STA_RS"	,"HIS_RS",	"Bow_RS",	"STA_eXpr",	"Bow_eXpr",	"STA_Sal",	"Bow_Sal",	"Sal",	"Kall",	'Sail')


df_long <- reshape2::melt(merged_df, 
                          id.vars = "gene_id", 
                          variable.name = "Group", 
                          value.name = "Value")

my_colors <- c("#669999","#669999","#669999","#EEB8C3","#EEB8C3","#9ECAE1","#9ECAE1","#9ECAE1","#CBB7FF","#CBB7FF","#CBB7FF","#FFEBAC","#FFEBAC","#FFEBAC","#FEB381","#FEB381","#99DFB9","#99DFB9","#99DFB9","#FB6501","#8A6800")

ggplot(df_long, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "CV", y = "Group") +
  scale_fill_manual(values = my_colors) +
  scale_y_discrete(limits = rev)

