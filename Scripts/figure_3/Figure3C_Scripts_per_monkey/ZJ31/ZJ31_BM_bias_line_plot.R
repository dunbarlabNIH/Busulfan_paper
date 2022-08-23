library(reshape2)
library(readxl)
library(barcodetrackR)
library(ggplot2)
library(tidyverse)
library(SummarizedExperiment)
library(mclust)
library(factoextra)
library(cluster)
library(ClusterR)
library(FactoMineR)
library(factoextra)
library(readxl)
library(dplyr)
library(ggpubr)


source("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/github_submission/helper_scripts/barcode_ggheatmap_annie_v2.1.R")
source("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/github_submission/helper_scripts/create_bias_line_graph.R")

# input information about monkey and samples
monkey_ID <- "ZJ31"
condition_type <- "TBI"
analysis_type <- "BM_bias"
process_date <- 20210302

# other parameters to change
foldchange_threshold <- 10
num_clones_checked <- 10

# load sample information and data 
# file paths
read_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/", monkey_ID, "_", condition_type, sep = "")
save_dir <- paste(read_dir, analysis_type, sep = "/")

# read in metadata file 
metadata_df <- read.delim(file = file.path(read_dir,paste(monkey_ID, "_",process_date,"_metadata_edited.txt", sep="")),
                          sep = "\t")

# read in counts file 
count_data <- read.delim(paste(read_dir, "/", monkey_ID, "_combined_", process_date, "_counts.txt", sep = ''), 
                         header = T, row.names = 1, sep = "\t")
count_data <- count_data[,metadata_df$SAMPLENAME]

# read in sample list 
sample_list <- read.delim(file = file.path(save_dir,"BM_bias_sample_list.txt"),
                          sep = "\t", header=F)
colnames(sample_list) <- "sample_name"
sample_list <- sample_list$sample_name

# find the corresponding sample indices
sample_index <- rep(NA, length(sample_list))
for (i in 1:length(sample_list)){
  sample_index[i] <- which(metadata_df$SAMPLENAME == sample_list[i])
}

metadata_df[["Cell.Type"]] <- metadata_df$simple_cell_type
metadata_df[sample_index[c(1,3,5,7)],]$Cell.Type <- "LBM"
metadata_df[sample_index[c(2,4,6,8)],]$Cell.Type <- "RBM"

# Subsetting data by only including the samples of interest 
metadata_subset <- metadata_df[sample_index,]
count_data_subset <- count_data[, sample_index]

# double checking to make sure subsetting worked
plot_labels_old <- paste( paste(metadata_df[sample_index, ]$months, "m", sep = ""), metadata_df[sample_index, ]$Cell.Type)
plot_labels <- paste( paste(metadata_subset$months, "m", sep = ""), metadata_subset$Cell.Type)
print(sum(plot_labels_old == plot_labels)) # this number should equal to the number of samples being analyzed 

your_SE <- barcodetrackR::create_SE(count_data_subset, metadata_subset, threshold = 0)

############### GETTING DATA TO CREATE BIAS LINE GRAPHS ###############
breaks_list <- c(10, 5, 2, 1.5, 0)

cell_type1 <- "LBM"
cell_type2 <- "RBM"
timepoints <- sort(intersect(your_SE$months[which(your_SE$Cell.Type == cell_type1)], 
                             your_SE$months[which(your_SE$Cell.Type == cell_type2)]))
bias_data <- bias_histogram(your_SE = your_SE,
                            split_bias_on = "Cell.Type",
                            bias_1 = cell_type1,
                            bias_2 = cell_type2,
                            split_bias_over = "months",
                            return_table = T, breaks = breaks_list)

png(file = paste(save_dir, "/", monkey_ID, "_biased_line_plot_", cell_type2, 
                 "_vs_", cell_type1, ".png", sep = ""), width=950, height=700)
create_bias_line_graph(bias_data, timepoints, paste(cell_type2, "(left)"), 
                       paste(cell_type1, "(right)"), textsize = 30) 
dev.off()

cell_type1 <- "RBM"
cell_type2 <- "LBM"
timepoints <- sort(intersect(your_SE$months[which(your_SE$Cell.Type == cell_type1)], 
                             your_SE$months[which(your_SE$Cell.Type == cell_type2)]))
bias_data <- bias_histogram(your_SE = your_SE,
                            split_bias_on = "Cell.Type",
                            bias_1 = cell_type1,
                            bias_2 = cell_type2,
                            split_bias_over = "months",
                            return_table = T, breaks = breaks_list)

png(file = paste(save_dir, "/", monkey_ID, "_biased_line_plot_", cell_type2, 
                 "_vs_", cell_type1, ".png", sep = ""), width=950, height=700)
create_bias_line_graph(bias_data, timepoints, paste(cell_type2, "(left)"), 
                       paste(cell_type1, "(right)"), textsize = 30)
dev.off()


