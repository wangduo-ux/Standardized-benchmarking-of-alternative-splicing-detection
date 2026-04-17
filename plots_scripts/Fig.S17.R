library(ggplot2)
library(ggpubr)

### DEI_file is input, which is available in the Source Data file of our paper

df <- read.table("DEI_file")

# -------------------------------
# Define a consistent theme for all plots
# -------------------------------
arial_theme <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major = element_line(color = "grey80"),  # Major grid lines
    panel.grid.minor = element_blank(),                 # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)  # Add border
  )

# -------------------------------
# Define colors for each sample comparison
# -------------------------------
cols <- c(
  "M8/D6"="#F36644",
  "F7/D6"="#FFC55D",
  "D5/D6"="#4AC5D9",
  "T1/D6"="#CC9749",
  "T2/D6"="#95B786"
)

# -------------------------------
# Generate a list of ggplot objects
# -------------------------------
plots <- list()
i <- 1
for (prefix in names(cols)) {
  depth_col <- paste0(prefix, "_depth")  # Column for sequencing depth
  deg_col   <- paste0(prefix, "_DEG")    # Column for number of DEIs
  
  # Scatter plot with linear fit and Pearson correlation
  p <- ggplot(df, aes_string(x = depth_col, y = deg_col)) +
    geom_point(color = cols[prefix]) +  # Points
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = cols[prefix]) +  # Linear fit
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4, family = "Arial") +  # Pearson correlation
    labs(
      x = paste("Depth in exonic regions of", prefix), 
      y = ifelse(i == 1, "No. of DEIs", "")  # Only first plot shows y-axis label
    ) +
    arial_theme
  
  plots[[i]] <- p
  i <- i + 1
}

# -------------------------------
# Arrange plots horizontally with labels
# -------------------------------
final_plot <- ggarrange(
  plotlist = plots,
  ncol = 5, nrow = 1,
  labels = c("b","c","d","e","f"),  # Subfigure labels
  label.x = 0, label.y = 0.99,
  font.label = list(size = 14, face = "bold")
)