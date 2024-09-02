library(ggplot2)
library(reshape2)
library(patchwork)  
library(cowplot)
library(data.table)
library(ggpubr)
library(fmsb)

rm(list = ls())
BRCA <- read.csv("outputs/tripletsBRCA.csv", row.names = NULL)
LIHC <- read.csv("outputs/tripletsLIHC.csv", row.names = NULL)
PRAD <- read.csv("outputs/tripletsPRAD.csv", row.names = NULL)

# Boxplots
color_palette <- c("#8FB5A3",  # Light green
                   "#A6C8DD",  # Light blue
                   "#DC7972",  # Light red
                   "#E3AE7E")  # Gold/yellow
datasets <- list(BRCA = BRCA, LIHC = LIHC, PRAD = PRAD)
plot_list <- list()
for (name in names(datasets)) {
  data <- datasets[[name]]
  data_box <- data[, c("Synergy", "Unique_X", "Unique_Y", "Redundancy")]
  df_long <- melt(data_box, measure.vars = c("Synergy", "Unique_X", "Unique_Y", "Redundancy"))
  
  plot <- ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = color_palette) +  # Set color palette to match the screenshot
    theme_classic() +
    ggtitle(name) +  
    theme(legend.position = "none")
  # Save plot to list
  plot_list[[name]] <- plot
}

# Create a plot with a legend
legend_plot <- ggplot(df_long, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2, aes(color = variable), size = 2, alpha = 0.6) +  # Add jitter points
  scale_fill_manual(values = color_palette) +  # Set fill color palette to match the screenshot
  scale_color_manual(values = color_palette) +  # Set point color to match fill colors
  theme_classic() +
  theme(
    legend.position = "right",  # Place legend on the right side, vertically aligned
    legend.title = element_text(size = 12),  # Set legend title font size
    legend.text = element_text(size = 10),   # Set legend text font size
    legend.key.size = unit(0.5, "cm")        # Set legend key size
  )
# Extract legend
legend <- cowplot::get_legend(legend_plot)
# Create a plot with only the legend
legend_only <- ggdraw(legend)

# Combine the three plots and the legend into a single row
final_plot <- plot_list$BRCA + plot_list$LIHC + plot_list$PRAD + legend_only + 
  plot_layout(ncol = 4, widths = c(1, 1, 1, 0.5))
# Save the combined plot as a PDF
ggsave("./plotoutputs/combined_plot_with_legend.pdf", plot = final_plot, width = 12, height = 6)


# Radar plots
miRNAx1 = "hsa-miR-19a-3p"
mRNAx1 = "RHOB"
RBPx1 = "ELAVL1"
miRNAx2 = "hsa-miR-93-5p"
mRNAx2 = "CALD1"
RBPx2 = "NOP58"

for(name in names(datasets))
{
  data <- datasets[[name]]
  result1 <- subset(data, miRNA == miRNAx1 & RBP == RBPx1 & mRNA == mRNAx1)
  result2 <- subset(data, miRNA == miRNAx2 & RBP == RBPx2 & mRNA == mRNAx2)
  
  if (nrow(result1) != 0) {
    var1 <- c(result1$Synergy, result1$Unique_X, result1$Unique_Y, result1$Redundancy)
    data <- matrix(var1,nrow=1)
    data <- as.data.frame(data)
    colnames(data) <- c("Synergy", "UniqueX","UniqueY","Redundancy")
    var2 <- rep(1.5,4)
    var3 <- rep(0,4)
    data <- rbind(var2,var3,data)
    
    pdf(file = paste0("./plotoutputs/Radar_Plot1_", name, ".pdf"))
    radarchart(data,
               axistype=1,seg=5,
               pty=16,pcol=rgb(0.2,0.5,0.5,0.9),plty=7,pfcol=rgb(0.2,0.5,0.5,0.5),plwd=3,
               cglcol="grey",cglty=1,axislabcol="black",cglwd=1.5,caxislabels=seq(0,1.5,0.3),
               vlcex=0.8,title=paste("Radar -", name))
    dev.off()
  }
  
  if (nrow(result2) != 0) {
    var1 <- c(result2$Synergy, result2$Unique_X, result2$Unique_Y, result2$Redundancy)
    data <- matrix(var1,nrow=1)
    data <- as.data.frame(data)
    colnames(data) <- c("Synergy", "UniqueX","UniqueY","Redundancy")
    var2 <- rep(1.5,4)
    var3 <- rep(0,4)
    data <- rbind(var2,var3,data)
    
    pdf(file = paste0("./plotoutputs/Radar_Plot2_", name, ".pdf"))
    radarchart(data,
               axistype=1,seg=5,
               pty=16,pcol=rgb(0.2,0.5,0.5,0.9),plty=7,pfcol=rgb(0.2,0.5,0.5,0.5),plwd=3,
               cglcol="grey",cglty=1,axislabcol="black",cglwd=1.5,caxislabels=seq(0,1.5,0.3),
               vlcex=0.8,title=paste("Radar -", name))
    dev.off()
  }
}

# Scatter plots
# Load necessary libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(cowplot)
rm(list = ls())
# Define the cancer types to be processed
cancers <- c("BRCA", "LIHC", "PRAD")

# Define miRNA, RBP, and mRNA names
miRNA_name <- "hsa-miR-19a-3p"
RBP_name <- "ELAVL1"
mRNA_name <- "RHOB"

# Initialize a list to store all plots
plot_list <- list()

# Loop through each cancer type
for (cancer in cancers) {
  # Load the data
  load(paste0("outputs/", cancer, "_data.RData"))
  
  # Log2 transform the expression data
  mRNAexpression <- log2(mRNAexpression + 1)
  miRNAexpression <- log2(miRNAexpression + 1)
  RBPexpression <- log2(RBPexpression + 1)
  
  # Get the expression data for the specified genes and miRNA
  miRNA_index <- which(rownames(miRNAexpression) == miRNA_name)
  RBP_index <- which(rownames(RBPexpression) == RBP_name)
  mRNA_index <- which(rownames(mRNAexpression) == mRNA_name)
  
  # Generate the first scatter plot: miRNA vs RBP
  miRNA_pair_df <- data.frame(
    miRNA_expr = as.numeric(miRNAexpression[miRNA_index, ]),
    RBP_expr = as.numeric(RBPexpression[RBP_index, ])
  )
  sp1 <- ggscatter(miRNA_pair_df, 
                   x = "miRNA_expr", 
                   y = "RBP_expr", 
                   color = "#999999", 
                   palette = "red", 
                   size = 1,
                   add = "reg.line",  
                   add.params = list(color = "red", fill = "lightgray"), 
                   conf.int = FALSE, 
                   rug = FALSE 
  ) + 
    stat_cor(method = "pearson", digits = 4) +  
    xlab(miRNA_name) + 
    ylab(RBP_name) + 
    border("black")
  
  # Generate the second scatter plot: RBP vs mRNA
  miRNA_mRNA_df <- data.frame(
    RBP_expr = as.numeric(RBPexpression[RBP_index, ]),
    mRNA_expr = as.numeric(mRNAexpression[mRNA_index, ])
  )
  sp2 <- ggscatter(miRNA_mRNA_df, 
                   x = "RBP_expr", 
                   y = "mRNA_expr", 
                   color = "#999999", 
                   palette = "red", 
                   size = 1,
                   add = "reg.line",  
                   add.params = list(color = "red", fill = "lightgray"), 
                   conf.int = FALSE, 
                   rug = FALSE 
  ) + 
    stat_cor(method = "pearson", digits = 4) + 
    xlab(RBP_name) + 
    ylab(mRNA_name) + 
    border("black")
  
  # Generate the third scatter plot: miRNA vs mRNA
  miRNA_mRNA_df <- data.frame(
    miRNA_expr = as.numeric(miRNAexpression[miRNA_index, ]),
    mRNA_expr = as.numeric(mRNAexpression[mRNA_index, ])
  )
  sp3 <- ggscatter(miRNA_mRNA_df, 
                   x = "miRNA_expr", 
                   y = "mRNA_expr", 
                   color = "#999999", 
                   palette = "red", 
                   size = 1,
                   add = "reg.line",  
                   add.params = list(color = "red", fill = "lightgray"), 
                   conf.int = FALSE, 
                   rug = FALSE 
  ) + 
    stat_cor(method = "pearson", digits = 4) + 
    xlab(miRNA_name) + 
    ylab(mRNA_name) + 
    border("black")
  
  # Add the generated plots to the list
  plot_list[[paste0(cancer, "_sp1")]] <- sp1
  plot_list[[paste0(cancer, "_sp2")]] <- sp2
  plot_list[[paste0(cancer, "_sp3")]] <- sp3
}

combined_plot <- (plot_list$BRCA_sp1 + plot_list$BRCA_sp2 + plot_list$BRCA_sp3) / 
  (plot_list$LIHC_sp1 + plot_list$LIHC_sp2 + plot_list$LIHC_sp3) / 
  (plot_list$PRAD_sp1 + plot_list$PRAD_sp2 + plot_list$PRAD_sp3)

# Save the combined plot as a PDF file
ggsave("./plotoutputs/Scatter_Plots1.pdf", plot = combined_plot, width = 10, height = 10)
rm(list = ls())
# Define the cancer types to be processed
cancers <- c("BRCA", "LIHC", "PRAD")

# Define miRNA, RBP, and mRNA names
miRNA_name <- "hsa-miR-93-5p"
RBP_name <- "NOP58"
mRNA_name <- "CALD1"
# Initialize a list to store all plots
plot_list <- list()

# Loop through each cancer type
for (cancer in cancers) {
  # Load the data
  load(paste0("outputs/", cancer, "_data.RData"))
  
  # Log2 transform the expression data
  mRNAexpression <- log2(mRNAexpression + 1)
  miRNAexpression <- log2(miRNAexpression + 1)
  RBPexpression <- log2(RBPexpression + 1)
  
  # Get the expression data for the specified genes and miRNA
  miRNA_index <- which(rownames(miRNAexpression) == miRNA_name)
  RBP_index <- which(rownames(RBPexpression) == RBP_name)
  mRNA_index <- which(rownames(mRNAexpression) == mRNA_name)
  
  # Generate the first scatter plot: miRNA vs RBP
  miRNA_pair_df <- data.frame(
    miRNA_expr = as.numeric(miRNAexpression[miRNA_index, ]),
    RBP_expr = as.numeric(RBPexpression[RBP_index, ])
  )
  sp1 <- ggscatter(miRNA_pair_df, 
                   x = "miRNA_expr", 
                   y = "RBP_expr", 
                   color = "#999999", 
                   palette = "red", 
                   size = 1,
                   add = "reg.line",  
                   add.params = list(color = "red", fill = "lightgray"), 
                   conf.int = FALSE, 
                   rug = FALSE 
  ) + 
    stat_cor(method = "pearson", digits = 4) +  
    xlab(miRNA_name) + 
    ylab(RBP_name) + 
    border("black")
  
  # Generate the second scatter plot: RBP vs mRNA
  miRNA_mRNA_df <- data.frame(
    RBP_expr = as.numeric(RBPexpression[RBP_index, ]),
    mRNA_expr = as.numeric(mRNAexpression[mRNA_index, ])
  )
  sp2 <- ggscatter(miRNA_mRNA_df, 
                   x = "RBP_expr", 
                   y = "mRNA_expr", 
                   color = "#999999", 
                   palette = "red", 
                   size = 1,
                   add = "reg.line",  
                   add.params = list(color = "red", fill = "lightgray"), 
                   conf.int = FALSE, 
                   rug = FALSE 
  ) + 
    stat_cor(method = "pearson", digits = 4) + 
    xlab(RBP_name) + 
    ylab(mRNA_name) + 
    border("black")
  
  # Generate the third scatter plot: miRNA vs mRNA
  miRNA_mRNA_df <- data.frame(
    miRNA_expr = as.numeric(miRNAexpression[miRNA_index, ]),
    mRNA_expr = as.numeric(mRNAexpression[mRNA_index, ])
  )
  sp3 <- ggscatter(miRNA_mRNA_df, 
                   x = "miRNA_expr", 
                   y = "mRNA_expr", 
                   color = "#999999", 
                   palette = "red", 
                   size = 1,
                   add = "reg.line",  
                   add.params = list(color = "red", fill = "lightgray"), 
                   conf.int = FALSE, 
                   rug = FALSE 
  ) + 
    stat_cor(method = "pearson", digits = 4) + 
    xlab(miRNA_name) + 
    ylab(mRNA_name) + 
    border("black")
  
  # Add the generated plots to the list
  plot_list[[paste0(cancer, "_sp1")]] <- sp1
  plot_list[[paste0(cancer, "_sp2")]] <- sp2
  plot_list[[paste0(cancer, "_sp3")]] <- sp3
}

combined_plot <- (plot_list$BRCA_sp1 + plot_list$BRCA_sp2 + plot_list$BRCA_sp3) / 
  (plot_list$LIHC_sp1 + plot_list$LIHC_sp2 + plot_list$LIHC_sp3) / 
  (plot_list$PRAD_sp1 + plot_list$PRAD_sp2 + plot_list$PRAD_sp3)

# Save the combined plot as a PDF file
ggsave("./plotoutputs/Scatter_Plots2.pdf", plot = combined_plot, width = 10, height = 10)
