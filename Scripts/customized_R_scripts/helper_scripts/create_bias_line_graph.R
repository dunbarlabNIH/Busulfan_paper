create_bias_line_graph <- function(temp_your_data, timepoints, bias1, bias2, linewidth = 3, pointsize = 5, textsize = 10) {
   
  added_prop_all <- data.frame(matrix(ncol = 3, nrow = 0))
  
  for (i in 1:length(temp_your_data)){
    
    timepoint_temp <- timepoints[i]
    temp_data <- temp_your_data[[i]]
    bins_list <- levels(temp_data$log2_bias_cuts)
    
    added_prop_vec <- c()
    for (j in 1:length(bins_list)){
      temp_index = which(temp_data$log2_bias_cuts == bins_list[j])
      temp_prop = sum(temp_data$added_proportions[temp_index])
      added_prop_vec <- c(added_prop_vec, temp_prop)
    }
    
    added_prop_df <- data.frame(bins_list,added_prop_vec, rep(timepoint_temp, length(bins_list)))
    
    added_prop_all <- rbind(added_prop_all, added_prop_df)
  }
  colnames(added_prop_all) <- c("bins", "added_prop", "timepoints")
  added_prop_all$timepoints <- factor(added_prop_all$timepoints, levels = timepoints)
  added_prop_all$bins <- factor(added_prop_all$bins, levels = bins_list)
  
  ggplot(data=added_prop_all, aes(x=bins, y=added_prop * 100, group=timepoints)) +
    geom_line(aes(color=timepoints), size = linewidth)+
    geom_point(aes(color=timepoints), size = pointsize) +
    theme_classic() + theme(
      # Remove panel border
      panel.border = element_blank(),  
      # Remove panel grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = textsize),
      axis.text.y = element_text(size = textsize),
      axis.text.x = ggplot2::element_text(vjust = 0.5, hjust = 1, angle = 90, size = textsize),
      legend.text=element_text(size=textsize),
      legend.title = element_text(size=textsize),
      axis.line = element_line(size = 2, linetype = "solid"), 
      axis.ticks.length = unit(0.75, "cm"), 
      axis.ticks = element_line(colour = "black", size = 2),
      plot.title = element_text(size = textsize, hjust = 0.5)
    ) + labs(color = "Timepoints", title = paste(bias1, "vs", bias2, "Population Bias"), y = "Added Proportions (%)") 

}

