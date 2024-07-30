library(ggplot2)
library(data.table)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(EnvStats)
library(hrbrthemes)
library(stringr)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

#Directory of the full table comprising all results
final_table_dir = "~/rprojects/review_data/results/new_generated_res_test/final_dataframe_humna.txt"
final_table = read.table(file = final_table_dir, sep = "\t",
                         header = TRUE)

df_cell = unique(final_table[,c("cell_type","trait")])

# Calculate the common cell types between diseases
disease_cell_type <- df_cell %>%
  group_by(trait) %>%
  summarize(cell_types = list(cell_type)) %>%
  ungroup()

# Create a matrix to store the names of common cell types
disease_names <- unique(df_cell$trait)
common_matrix <- matrix("", nrow = length(disease_names), ncol = length(disease_names))
rownames(common_matrix) <- disease_names
colnames(common_matrix) <- disease_names

# Matrix to store the counts of common cell types
count_matrix <- matrix(0, nrow = length(disease_names), ncol = length(disease_names))
rownames(count_matrix) <- disease_names
colnames(count_matrix) <- disease_names

for (i in 1:length(disease_names)) {
  for (j in 1:length(disease_names)) {
    if (i < j) {
      cell_types_i <- unlist(disease_cell_type$cell_types[disease_cell_type$trait == disease_names[i]])
      cell_types_j <- unlist(disease_cell_type$cell_types[disease_cell_type$trait == disease_names[j]])
      common_cell_types <- intersect(cell_types_i, cell_types_j)
      common_matrix[i, j] <- paste(common_cell_types, collapse = "\n")
      common_matrix[j, i] <- common_matrix[i, j]
      count_matrix[i, j] <- length(common_cell_types)
      count_matrix[j, i] <- length(common_cell_types)
    }
  }
}



# Replace empty strings with NA for better visualization
common_matrix[common_matrix == ""] <- NA

# Create a color function for the heatmap
col_fun <- colorRamp2(c(0, max(count_matrix)), c("white", "red"))

max_count <- max(count_matrix)

# Plot the heatmap with common cell types names
heatmap_plot <- Heatmap(count_matrix,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.rect(x, y, width - unit(1, "mm"), height - unit(1, "mm"),
                                    gp = gpar(col = "black", fill = col_fun(count_matrix[i, j])))
                          if (!is.na(common_matrix[i, j])) {
                            grid.text(common_matrix[i, j], x, y, gp = gpar(fontsize = 28))
                          }
                        },
                        cluster_rows = FALSE, cluster_columns = FALSE,
                        name = "Number of Common Cell Types",
                        col = col_fun,
                        row_names_side = "left", column_names_side = "top",
                        show_heatmap_legend = TRUE,
                        row_names_gp = gpar(fontsize = 28, fontface = "bold"),  # Adjust row names font
                        column_names_gp = gpar(fontsize = 28, fontface = "bold"),  # Adjust column names font
                        heatmap_legend_param = list(title = "Common \nCell Types\n",
                                                    title_gp = gpar(fontsize = 28, fontface = "bold"),
                                                    at = 0:max_count,  # Adjust legend ticks
                                                    labels = 0:max_count,
                                                    labels_gp = gpar(fontsize = 28)))  # Adjust legend title font


# Draw the heatmap
#draw(heatmap_plot)

#Generating Figure 4 for common cell types

#Address of the cell type part of figure 4
heat_dir = "~/rprojects/review_data/results/new_generated_res_test/traitVStrait_cell.png"
file.remove(heat_dir)
png(heat_dir,width = 8400, height = 7300, res = 300)
#png(heat_dir,width = 12800, height = 8600, res = 300)
draw(heatmap_plot)

dev.off()

#####################################################################################

df_tf = unique(final_table[,c("TF","trait")])

# Calculate the common TFs between diseases
disease_tf <- df_tf %>%
  group_by(trait) %>%
  summarize(TFs = list(TF)) %>%
  ungroup()

# Create a matrix to store the names of common TFs
disease_names <- unique(df_tf$trait)
common_matrix <- matrix("", nrow = length(disease_names), ncol = length(disease_names))
rownames(common_matrix) <- disease_names
colnames(common_matrix) <- disease_names

# Matrix to store the counts of common TFs
count_matrix <- matrix(0, nrow = length(disease_names), ncol = length(disease_names))
rownames(count_matrix) <- disease_names
colnames(count_matrix) <- disease_names

for (i in 1:length(disease_names)) {
  for (j in 1:length(disease_names)) {
    if (i < j) {
      TFs_i <- unlist(disease_tf$TFs[disease_tf$trait == disease_names[i]])
      TFs_j <- unlist(disease_tf$TFs[disease_tf$trait == disease_names[j]])
      common_TFs <- intersect(TFs_i, TFs_j)
      common_matrix[i, j] <- paste(common_TFs, collapse = "\n")
      common_matrix[j, i] <- common_matrix[i, j]
      count_matrix[i, j] <- length(common_TFs)
      count_matrix[j, i] <- length(common_TFs)
    }
  }
}



# Replace empty strings with NA for better visualization
common_matrix[common_matrix == ""] <- NA

# Create a color function for the heatmap
col_fun <- colorRamp2(c(0, max(count_matrix)), c("white", "red"))

max_count <- max(count_matrix)

# Plot the heatmap with common TFs names
heatmap_plot <- Heatmap(count_matrix,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.rect(x, y, width - unit(1, "mm"), height - unit(1, "mm"),
                                    gp = gpar(col = "black", fill = col_fun(count_matrix[i, j])))
                          if (!is.na(common_matrix[i, j])) {
                            grid.text(common_matrix[i, j], x, y, gp = gpar(fontsize = 28))
                          }
                        },
                        cluster_rows = FALSE, cluster_columns = FALSE,
                        name = "Number of Common \nTranscription Factors",
                        col = col_fun,
                        row_names_side = "left", column_names_side = "top",
                        show_heatmap_legend = TRUE,
                        row_names_gp = gpar(fontsize = 28, fontface = "bold"),  # Adjust row names font
                        column_names_gp = gpar(fontsize = 28, fontface = "bold"),  # Adjust column names font
                        heatmap_legend_param = list(title = "Common \nTranscription Factors\n",
                                                    title_gp = gpar(fontsize = 28, fontface = "bold"),
                                                    at = 0:max_count,  # Adjust legend ticks
                                                    labels = 0:max_count,
                                                    labels_gp = gpar(fontsize = 28)))  # Adjust legend title font


# Draw the heatmap
#draw(heatmap_plot)

#Generating Figure 4 for common cell types

#Address of the TF part of figure 4
heat_dir = "~/rprojects/review_data/results/new_generated_res_test/traitVStrait_TF.png"
file.remove(heat_dir)
png(heat_dir,width = 8400, height = 7300, res = 300)
#png(heat_dir,width = 12800, height = 8600, res = 300)
draw(heatmap_plot)

dev.off()
