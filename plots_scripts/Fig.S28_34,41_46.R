library(tidyverse)
library(ggpubr)

# -------------------------------
# Read the raw stratified analysis data
# -------------------------------
df <- read.table("input data")

### input data are RMSE, CV, or SD results from stratified analysis of 21 isoform quantification pipelines and SUPPA2-based event analysis pipelines. These data are available in the Source Data file in our paper 

# -------------------------------
# Define software factor levels
# -------------------------------

### for isoform_level analysis
software_levels <- c("STAR_StringTie (Assembly)","HISAT2_StringTie (Assembly)","Subjunc_StringTie (Assembly)","STAR_Cuffdiff","HISAT2_Cuffdiff","STAR_StringTie","HISAT2_StringTie","Subjunc_StringTie",
                     "STAR_featureCounts","HISAT2_featureCounts","Subjunc_featureCounts",
                     "STAR_RSEM","HISAT2_RSEM","Bowtie2_RSEM",
                     "STAR_eXpress","Bowtie2_eXpress",
                     "STAR_Salmon","Bowtie2_Salmon","Salmon",
                     "Kallisto","Sailfish")

### for SUPPA2-based event_level analysis
software_levels <- c(
  "STAR_StringTie (Assembly)_SUPPA2", "HISAT2_StringTie (Assembly)_SUPPA2", "Subjunc_StringTie (Assembly)_SUPPA2",
  "STAR_Cuffdiff_SUPPA2", "HISAT2_Cuffdiff_SUPPA2", "STAR_StringTie_SUPPA2", "HISAT2_StringTie_SUPPA2", "Subjunc_StringTie_SUPPA2",
  "STAR_featureCounts_SUPPA2", "HISAT2_featureCounts_SUPPA2", "Subjunc_featureCounts_SUPPA2",
  "STAR_RSEM_SUPPA2", "HISAT2_RSEM_SUPPA2", "Bowtie2_RSEM_SUPPA2",
  "STAR_eXpress_SUPPA2", "Bowtie2_eXpress_SUPPA2",
  "STAR_Salmon_SUPPA2", "Bowtie2_Salmon_SUPPA2", "Salmon_SUPPA2",
  "Kallisto_SUPPA2", "Sailfish_SUPPA2"
)



# -------------------------------
# Define group levels for plotting
# -------------------------------

############  isoform level analysis #################

# (1) k value
group_order <- c("<4",	"4-5.9"	,"5.9-8.5",	"8.5-13.8",	">13.8") ## Quartet reference datasets based accuracy evaluation
group_order <- c("<4",	"4-8",	"8-12",	"12-24",	">24")  ## reproducibility evaluation using CV


# (2) isoform length
group_order <- c("<1520bp",	"1520-2500bp",	"2500-3697bp",	"3697-5622bp",	">5622bp")  ## Quartet reference datasets based accuracy evaluation
group_order <- c("<600bp",	"600-1500bp",	"1500-3000bp",	"3000-5000bp"	,">5000bp")  ## reproducibility evaluation using CV

# (3) isoformnumber
group_order <- c("1",	"2-5",	"6-10",	"11-18",	">18")## Quartet reference datasets based accuracy evaluation
group_order <- c("<5",	"6-10",	"11-15",	"16-24",	">24")  ## reproducibility evaluation using CV

# (4) exon number
group_order <- c("<4",	"4-5",	"6-8",	"9-13",	">13") ## Quartet reference datasets based accuracy evaluation
group_order <- c("<3",	"3-4",	"5-7",	"8-11",	">11")  ## reproducibility evaluation using CV

# (5) exon length
group_order <- c("<172bp",	"172-243bp",	"243-361bp",	"361-646bp",	">646bp") ## Quartet reference datasets based accuracy evaluation
group_order <- c("<140bp",	"140-200bp",	"200-300bp",	"300-500bp",	">450bp")  ## reproducibility evaluation using CV

# (6) GC content
group_order <- c("<40.6%",	"40.6-45.0%",	"45.0-50.5%",	"50.5-56.7%",	">56.7%") ## Quartet reference datasets based accuracy evaluation
group_order <- c("<42%",	"42-47%",	"47-52%",	"52-58%",	">58%")  ## reproducibility evaluation using CV

# (7) isoform expression 
group_order <- c("<q20",	"q20-q40",	"q40-q60",	"q60-q80",	">q80") 


colnames(df) <- group_order
############  event level analysis #################

# (1) overlapped event number
group_order <- c("<4",	"4-8",	"9-15",	"16-35",	">35")

# (2) isoform number
group_order <- c("<3",	"3",	"4-5",	"6-9",	">9")

# (3) GC content
group_order <- c("<43%",	"43%-49%",	"49%-55%",	"55%-61%",	">61%")

# (4) Mappability
group_order <- c("<0.2",	"0.2-0.5"	,"0.5-0.7",	"0.7-0.9",	">0.9")

# (5) PSI and isoform expression 
group_order <- c("<q20",	"q20-q40",	"q40-q60",	"q60-q80",	">q80") 


colnames(df) <- group_order



# Assign software names repeated for each group
df$software <- rep(software_levels, each = 18)
# -------------------------------
# Transform from wide to long format
# -------------------------------
df_long <- df %>%
  pivot_longer(
    cols = all_of(group_order),  # pivot only the five groups
    names_to = "Group",
    values_to = "RMSE"
  ) %>%
  mutate(
    Group = factor(Group, levels = group_order),          # ensure correct order
    software = factor(software, levels = software_levels) # correct factor order
  )

# -------------------------------
# Create boxplot with jittered points
# -------------------------------
p1 <- ggplot(df_long, aes(x = Group, y = RMSE, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    shape = 21,
    color = "black",
    stroke = 0.1,
    width = 0.2,
    size = 1.5,
    alpha = 0.7
  ) +
  facet_wrap(~ software, ncol = 8, scales = "free_y", labeller = label_wrap_gen(width = 8)) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    axis.title.x = element_blank()
  )

# -------------------------------
# Combine plots (example p2 placeholder)
# -------------------------------
# Replace `p2` with your second plot if available, and replace "RMSE" with "cv" or "SD" according to your inuput data

p <- ggarrange(
  p1, p2,  # NULL is a placeholder for p2
  ncol = 1, nrow = 2,
  labels = c("a", "b"),       # subplot labels
  label.x = 0, label.y = 0.99,
  font.label = list(size = 14, face = "bold")
)

