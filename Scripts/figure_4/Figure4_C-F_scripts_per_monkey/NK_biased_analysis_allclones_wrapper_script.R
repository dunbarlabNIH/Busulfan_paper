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
library(rlang)

analysis_type <- "NK_bias"

TBI_monkeys <- c("ZJ31", "ZG66", "ZH33", "ZK22")
Busulfan_monkeys <- c("H84D", "M10U004", "M11021142")

min_proportion_threshold <- 0.001

######################### TBI MONKEY ANALYSIS #########################
biased_clone_cum_contrib_TBI <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in 1:length(TBI_monkeys)){
  monkey_ID <- TBI_monkeys[i]
  condition_type <- "TBI"
  read_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/", 
                    monkey_ID, "_", condition_type, "/", analysis_type, sep = "")
  biased_clone_contribution <- readRDS(paste(read_dir, "/", monkey_ID, "_NK_biased_cum_contribution_threshold", 
                                             min_proportion_threshold, ".RDS", sep = ""))
  
  biased_clone_cum_contrib_TBI <- rbind(biased_clone_cum_contrib_TBI, biased_clone_contribution)
  
  # biased_clone_barcode <- readRDS(paste(read_dir, "/", monkey_ID, 
  #                                       "_NK_biased_clone_barcode_manual.RDS", sep = ""))
  # assign(paste(monkey_ID, "_biased_clone_num", sep = ""), biased_clone_num)
  # assign(paste(monkey_ID, "_biased_clone_barcode", sep = ""), biased_clone_barcode)
}

biased_clone_cum_contrib_TBI$Timepoint <- as.numeric(biased_clone_cum_contrib_TBI$Timepoint)
biased_clone_cum_contrib_TBI$Cum_Contribution <- as.numeric(biased_clone_cum_contrib_TBI$Cum_Contribution)
biased_clone_cum_contrib_TBI[["Condition_type"]] <- rep("TBI", length(biased_clone_cum_contrib_TBI$Timepoint ))

ggplot(data=biased_clone_cum_contrib_TBI, aes(x=Timepoint, y=Cum_Contribution, group=Monkey)) + 
  geom_line(aes(color=Monkey))+ geom_point(aes(color=Monkey))+theme_classic()

######################### Busulfan MONKEY ANALYSIS #########################

biased_clone_cum_contrib_busulfan <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in 1:length(Busulfan_monkeys)){
  monkey_ID <- Busulfan_monkeys[i]
  condition_type <- "Busulfan"
  read_dir <- paste("/Volumes/dirgroup/TSCBB/LAB/Joint_Diana_Joy/For_Annie/monkey_data/", 
                    monkey_ID, "_", condition_type, "/", analysis_type, sep = "")
  biased_clone_contribution <- readRDS(paste(read_dir, "/", monkey_ID, "_NK_biased_cum_contribution_threshold", 
                                             min_proportion_threshold, ".RDS", sep = ""))
  biased_clone_cum_contrib_busulfan <- rbind(biased_clone_cum_contrib_busulfan, biased_clone_contribution)
  
  # biased_clone_barcode <- readRDS(paste(read_dir, "/", monkey_ID, 
  #                                       "_NK_biased_clone_barcode_manual.RDS", sep = ""))
  # assign(paste(monkey_ID, "_biased_clone_num", sep = ""), biased_clone_num)
  # assign(paste(monkey_ID, "_biased_clone_barcode", sep = ""), biased_clone_barcode)
}

biased_clone_cum_contrib_busulfan$Timepoint <- as.numeric(biased_clone_cum_contrib_busulfan$Timepoint)
biased_clone_cum_contrib_busulfan$Cum_Contribution <- as.numeric(biased_clone_cum_contrib_busulfan$Cum_Contribution)
biased_clone_cum_contrib_busulfan[["Condition_type"]] <- rep("Busulfan", length(biased_clone_cum_contrib_busulfan$Timepoint ))

ggplot(data=biased_clone_cum_contrib_busulfan, aes(x=Timepoint, y=Cum_Contribution, group=Monkey)) + 
  geom_line(aes(color=Monkey))+ geom_point(aes(color=Monkey))+theme_classic()

######################### COMBINED ANALYSIS #########################

biased_clone_cum_contrib_all <- rbind(biased_clone_cum_contrib_TBI, biased_clone_cum_contrib_busulfan)

biased_clone_cum_contrib_all$Condition_type <- factor(biased_clone_cum_contrib_all$Condition_type, levels = c("TBI", "Busulfan"))

monkey_list <- unique(biased_clone_cum_contrib_all$Monkey)
color_list <- c("ZJ31" = "black", "ZG66" = "black", "ZH33" = "black", "ZK22" = "black", "H84D" = "green", "M10U004" = "red", "M11021142" = "blue")
shape_list <- c("ZJ31" = 15, "ZG66" = 19, "ZH33" = 17, "ZK22" = 18, "H84D" = 19, "M10U004" = 19, "M11021142" = 19)

ggplot(data=biased_clone_cum_contrib_all, aes(x=Timepoint, y=Cum_Contribution, group=Monkey)) + 
  geom_line(aes(color=Monkey))+ geom_point(aes(color=Monkey, shape=Monkey))+
  theme_classic() + scale_color_manual(values=color_list) + scale_shape_manual(values=shape_list) + 
  ggplot2::labs(title = paste("Cumulative Contribution of biased clones\nat each timepoint; threshold: ",  
                              min_proportion_threshold * 100, "%",sep = "")) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), text = ggplot2::element_text(size = 20)) + 
  ggplot2::scale_y_continuous(name = "Cumulative Contribution (%)", labels = function(i) (paste0(i * 100, "%")), limits = c(0, 1)) + 
  scale_x_continuous(limits = c(0, 40))


