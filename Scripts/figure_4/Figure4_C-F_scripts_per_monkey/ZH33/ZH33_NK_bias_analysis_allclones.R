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

source("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/github_submission/helper_scripts/barcode_ggheatmap_annie_v2.1.R")
source("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/github_submission/helper_scripts/barcode_ggheatmap_annie_v2.2.R")

# input information about monkey and samples
monkey_ID <- "ZH33"
condition_type <- "TBI"
analysis_type <- "NK_bias"
process_date <- 20220426

# other parameters to change
foldchange_threshold <- 10
min_proportion_threshold <- 0.001

# load sample information and data 
# file paths
read_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/", monkey_ID, "_", condition_type, sep = "")
save_dir <- paste(read_dir, analysis_type, sep = "/")

# read in metadata file 
metadata_df <- read.delim(file = file.path(read_dir,paste(monkey_ID, "_combined_",process_date,"_metadata_cleaned.txt", sep="")),
                          sep = "\t")

# read in counts file 
count_data <- read.delim(paste(read_dir, "/", monkey_ID, "_combined_", process_date, "_counts_cleaned.txt", sep = ''), 
                         header = T, row.names = 1, sep = "\t")
count_data <- count_data[,metadata_df$SAMPLENAME]

# read in sample list 
sample_list <- read.delim(file = file.path(save_dir,"NK_bias_sample_list.txt"),
                          sep = "\t", header=F)
colnames(sample_list) <- "sample_name"
sample_list <- sample_list$sample_name

# find the corresponding sample indices
sample_index <- rep(NA, length(sample_list))
for (i in 1:length(sample_list)){
  sample_index[i] <- which(metadata_df$SAMPLENAME == sample_list[i])
}

# Subsetting data by only including the samples of interest 
metadata_subset <- metadata_df[sample_index,]
count_data_subset <- count_data[, sample_index]

# double checking to make sure subsetting worked
plot_labels_old <- paste( paste(metadata_df[sample_index, ]$Month_timepoint, "m", sep = ""), metadata_df[sample_index, ]$simple_cell_type)
plot_labels <- paste( paste(metadata_subset$Month_timepoint, "m", sep = ""), metadata_subset$simple_cell_type)
print(sum(plot_labels_old == plot_labels)) # this number should equal to the number of samples being analyzed 

timepoints_unique <- unique(metadata_subset$Month_timepoint)

timepoint_index_list <- list()

for (i in 1:length(timepoints_unique)){
  timepoint_temp <- timepoints_unique[i]
  timepoint_indices <- which(metadata_subset$Month_timepoint == timepoint_temp)
  timepoint_indices_df <- data.frame(timepoint_indices)
  colnames(timepoint_indices_df) <- timepoint_temp
  #print(timepoint_indices_df)
  timepoint_index_list <- c(timepoint_index_list, timepoint_indices_df)
  #print(timepoint_index_list)
}

your_SE <- barcodetrackR::create_SE(count_data_subset, metadata_subset, threshold = 0)

NK_bias_count_all <- list()
NK_bias_barcode_all <- list()
NK_bias_contribution_all <- list()

for (i in 1:length(timepoints_unique)){
  
  timepoint_temp <- toString(timepoints_unique[i])
  SE_subset <- your_SE[, timepoint_index_list[[timepoint_temp]]]
  CD16_index <- which(SE_subset$simple_cell_type == "CD16+")
  
  SE_subset_proportions_assay <- SummarizedExperiment::assays(SE_subset)[["proportions"]]
  colnames(SE_subset_proportions_assay) <- SE_subset$simple_cell_type
  
  proportions_assay_filtered <- SE_subset_proportions_assay %>% 
    tibble::rownames_to_column('barcodes') %>%
    dplyr::filter(.[["CD16+"]] >= min_proportion_threshold) %>%
    tibble::column_to_rownames('barcodes')
  
  CD16_proportions_filtered <- proportions_assay_filtered$`CD16+`
  
  biased_clones_all <- c()
  
  for (j in 1:length(SE_subset$simple_cell_type)){
    
    if (SE_subset$simple_cell_type[j] != "CD16+"){
      #print(j)
      temp_assay <- data.frame(proportions_assay_filtered[, j])
      rownames(temp_assay) <- rownames(proportions_assay_filtered)
      #temp_logs_assay <- temp_logs_assay[rownames(CD16_logs_assay), ]
      fold_change_indiv_comp <- CD16_proportions_filtered / temp_assay
      biased_clone_temp_index <- which(fold_change_indiv_comp >= foldchange_threshold)
      biased_clone_ID <- rownames(proportions_assay_filtered)[biased_clone_temp_index]
      #print(biased_clone_ID)
      biased_clones_all <- c(biased_clones_all, biased_clone_ID)
      
    }
  }
  
  biased_clones_table <- table(biased_clones_all)
  #print(biased_clones_table)
  biased_clones_final_barcode <- names(which(biased_clones_table >= length(unique(SE_subset$simple_cell_type)) - 1))
  biased_clones_final_num <- length(biased_clones_final_barcode)
  #print(biased_clones_final_num)
  
  biased_clones_final_barcode_list <- list(biased_clones_final_barcode)
  names(biased_clones_final_barcode_list) <- timepoint_temp
  NK_bias_barcode_all <- c(NK_bias_barcode_all, biased_clones_final_barcode_list)
  
  biased_clones_final_list <- list(biased_clones_final_num)
  names(biased_clones_final_list) <- timepoint_temp
  NK_bias_count_all <- c(NK_bias_count_all, biased_clones_final_list)
  
}

NK_bias_cum_contribution_CD16 <- data.frame(matrix(nrow = 0, ncol = 3))

for (i in 1:length(timepoints_unique)){
  
  timepoint_temp <- toString(timepoints_unique[i])
  SE_subset <- your_SE[, timepoint_index_list[[timepoint_temp]]]
  CD16_index <- which(SE_subset$simple_cell_type == "CD16+")
  
  biased_clones_barcode <- NK_bias_barcode_all[[timepoint_temp]]
  cum_contrib <- sum(SummarizedExperiment::assays(SE_subset)[["proportions"]][biased_clones_barcode, CD16_index])
  NK_bias_cum_contribution_CD16 <- rbind(NK_bias_cum_contribution_CD16, c(cum_contrib, timepoint_temp, monkey_ID))
  colnames(NK_bias_cum_contribution_CD16) <- c("Cum_Contribution", "Timepoint", "Monkey")
  
}

saveRDS(NK_bias_cum_contribution_CD16, file = paste(save_dir, "/", monkey_ID, "_NK_biased_cum_contribution_threshold", min_proportion_threshold, ".RDS", sep = ""))

saveRDS(NK_bias_barcode_all, file = paste(save_dir, "/", monkey_ID, "_NK_biased_clone_barcode_threshold", min_proportion_threshold,".RDS", sep = ""))
saveRDS(NK_bias_count_all, file = paste(save_dir, "/", monkey_ID, "_NK_biased_clone_count_threshold", min_proportion_threshold, ".RDS", sep = ""))




