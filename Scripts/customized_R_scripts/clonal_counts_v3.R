library(dplyr)
library(barcodetrackR)
library (rlang)
library(magrittr)
library(tidyr)
library(tibble)
library(SummarizedExperiment)

clonal_count_v3 <- function(your_SE,
                         percent_threshold = 0,
                         plot_over,
                         plot_over_display_choices = NULL,
                         keep_numeric = TRUE,
                         group_by,
                         group_by_choices = NULL,
                         cumulative = FALSE,
                         point_size = 3,
                         line_size = 2,
                         text_size = 12,
                         your_title = NULL,
                         return_table = FALSE, 
                         plot_total_cum = FALSE) {
  
  # Some basic error checking before running the function
  coldata_names <- colnames(SummarizedExperiment::colData(your_SE))
  
  if (!(plot_over %in% coldata_names)) {
    stop("plot_over must match a column name in colData(your_SE)")
  }
  if (!(group_by %in% coldata_names)) {
    stop("group_by must match a column name in colData(your_SE)")
  }
  
  if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
    plot_over_display_choices <- plot_over_display_choices %||% sort(unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
    plot_over_display_choices <- as.numeric(as.character(plot_over_display_choices))
  } else if (is.factor(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
    plot_over_display_choices <- plot_over_display_choices %||% factor(SummarizedExperiment::colData(your_SE)[[plot_over]], levels = levels(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  } else {
    plot_over_display_choices <- plot_over_display_choices %||% factor(SummarizedExperiment::colData(your_SE)[[plot_over]], levels = unique(SummarizedExperiment::colData(your_SE)[[plot_over]]))
  }
  
  group_by_choices <- group_by_choices %||% levels(as.factor(SummarizedExperiment::colData(your_SE)[[group_by]]))
  
  # More error handling
  if (!all(plot_over_display_choices %in% SummarizedExperiment::colData(your_SE)[, plot_over])) {
    stop("All elements of plot_over_display_choices must match cvalues in plot_over column")
  }
  if (!all(group_by_choices %in% SummarizedExperiment::colData(your_SE)[, group_by])) {
    stop("All elements of group_by_choices must match values in group_by column")
  }
  
  # extract metadata
  temp_subset <- your_SE[, (your_SE[[plot_over]] %in% plot_over_display_choices)]
  # Keep only the data included in group_by_choices
  temp_subset <- temp_subset[, (temp_subset[[group_by]] %in% group_by_choices)]
  temp_subset_coldata <- SummarizedExperiment::colData(temp_subset) %>% tibble::as_tibble()
  your_data <- SummarizedExperiment::assays(temp_subset)[["proportions"]]
  your_data <- your_data[rowSums(your_data) > 0, , drop = FALSE]
  
  # calculate measure for each sample
  if (!cumulative) {
    calculated_index <- colSums(your_data > percent_threshold) %>%
      tibble::enframe(name = "SAMPLENAME", value = "index") %>%
      dplyr::mutate(index_type = "unique barcode count")
    ylabel <- "unique barcode count"
  } else {
    
    # Get count data into tidy format (and remove zeros)
    tidy_counts <- your_data %>%
      dplyr::mutate(barcode_seq = rownames(your_data)) %>%
      dplyr::select(.data$barcode_seq, everything()) %>%
      tidyr::pivot_longer(cols = seq_len(ncol(your_data)) + 1, names_to = "sample", values_to = "count") %>%
      dplyr::filter(count > 0)
    
    # Add group by variable into tidy counts
    tidy_counts$group_by <- temp_subset_coldata[match(tidy_counts$sample, temp_subset_coldata$SAMPLENAME), ][[group_by]]
    
    # Add plot over variable into tidy counts
    tidy_counts$plot_over <- temp_subset_coldata[match(tidy_counts$sample, temp_subset_coldata$SAMPLENAME), ][[plot_over]]
    # tidy_counts$plot_over <- factor(tidy_counts$plot_over, levels = levels(plot_over_display_choices))
    
    # Order tidy counts by desired order
    if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
      tidy_counts_ordered <- tidy_counts %>%
        dplyr::arrange(group_by, plot_over)
    } else {
      tidy_counts_ordered <- tidy_counts %>%
        dplyr::arrange(group_by, factor(plot_over, levels = levels(plot_over_display_choices)))
    }
    
    # Only keep the first occurence of each barcode within each group_by category
    tidy_counts_filtered <- tidy_counts_ordered %>%
      dplyr::group_by(group_by) %>%
      dplyr::distinct(.data$barcode_seq, .keep_all = T)
    
    
    # Summarize number of new barcodes and cumulative barcodes at each timepoint

    if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
      summarized_data <- tidy_counts_filtered %>%
        dplyr::group_by(.data$group_by, plot_over, .data$sample) %>%
        dplyr::summarise(new_count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(group_by) %>%
        dplyr::mutate(cumulative_count = cumsum(.data$new_count))
    } else {
      summarized_data <- tidy_counts_filtered %>%
        dplyr::group_by(group_by, factor(plot_over, levels = levels(plot_over_display_choices)), sample) %>%
        dplyr::summarise(new_count = dplyr::n(), .groups = "drop") %>%
        dplyr::group_by(group_by) %>%
        dplyr::mutate(cumulative_count = cumsum(.data$new_count))
      colnames(summarized_data)[2] <- "plot_over"
    }

    
    if (length(unique(tidy_counts_ordered$sample)) != length(unique(summarized_data$sample))){
      print("during filtering step, samples without unique barcodes were removed...manually adding those samples back in")
      samples_removed <- setdiff(tidy_counts_ordered$sample, summarized_data$sample)
      
      for (i in 1:length(samples_removed)){
        index_sample_removed <- which(your_SE$SAMPLENAME == samples_removed[i])
        summarized_data[dim(summarized_data)[1] + 1, 3] <- samples_removed[i]
        summarized_data[dim(summarized_data)[1], 4] <- 0 
        summarized_data[dim(summarized_data)[1], 1] <- your_SE[[group_by]][index_sample_removed]
        summarized_data[dim(summarized_data)[1], 2] <- your_SE[[plot_over]][index_sample_removed]
        summarized_data[dim(summarized_data)[1], 5] <- 0
        max_count <- max(summarized_data$cumulative_count[summarized_data$group_by == your_SE[[group_by]][index_sample_removed]])
        summarized_data[dim(summarized_data)[1], 5] <- summarized_data$new_count[dim(summarized_data)[1]] + as.numeric(max_count)
      }
      
    }
    
    if (plot_total_cum == T){
      ## Calculating total cumulative counts
      # Only keep the first occurence of each barcode within each group_by category
      tidy_counts_filtered_total <- tidy_counts_ordered %>%
        dplyr::ungroup() %>%
        dplyr::distinct(.data$barcode_seq, .keep_all = T)
      
      if (is.numeric(SummarizedExperiment::colData(your_SE)[[plot_over]])) {
        summarized_data_total <- tidy_counts_filtered_total %>%
          dplyr::group_by(plot_over) %>%
          dplyr::summarise(new_count = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(cumulative_count = cumsum(.data$new_count)) %>%
          dplyr::mutate(group_by = "total", .before = plot_over) %>%
          dplyr::mutate(sample = paste("total_clonal_count", .data$plot_over, sep = "_"), .before = new_count)
        
      } else {
        summarized_data_total <- tidy_counts_filtered_total %>%
          dplyr::group_by(plot_over) %>%
          dplyr::summarise(new_count = dplyr::n(), .groups = "drop") %>%
          dplyr::mutate(cumulative_count = cumsum(.data$new_count))%>%
          dplyr::mutate(group_by = "total", .before = plot_over) %>%
          dplyr::mutate(sample = paste("total_clonal_count", .data$plot_over, sep = "_"), .before = new_count)
        colnames(summarized_data_total)[2] <- "plot_over"
      }
      
      summarized_data <- rbind(summarized_data, summarized_data_total)
    }else{
      summarized_data <- summarized_data
    }
    
    # Put into proper structure
    calculated_index <- as.data.frame(summarized_data) %>%
      dplyr::select(.data$sample, .data$cumulative_count) %>%
      dplyr::rename(SAMPLENAME = .data$sample, index = .data$cumulative_count) %>%
      dplyr::mutate(index_type = "unique barcodes") %>%
      dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME))

    # Fix the fact that certain samples will be dropped if they have 0 new clones
    # # Hacky way to do it. Need to go back and make it better later.
    # if (length(unique(tidy_counts_filtered$sample)) < length(unique(tidy_counts_ordered$sample))) {
    #   samples_not_found <- unique(tidy_counts_ordered$sample)[unique(tidy_counts_ordered$sample) %in% unique(tidy_counts_filtered$sample) == FALSE]
    #   calculated_index <- rbind(calculated_index, data.frame(
    #     SAMPLENAME = samples_not_found,
    #     index = rep(0, length(samples_not_found)),
    #     index_type = rep("unique barcodes", length(samples_not_found))
    #   ))
    #   # Reorder
    #   calculated_index <- calculated_index %>%
    #     dplyr::arrange(factor(.data$SAMPLENAME), levels = unique(tidy_counts_ordered$sample))
    #   
    #   # Give the missing samples the same cumulative count as above samples
    #   for (i in seq_along(nrow(calculated_index))) {
    #     if (calculated_index[i, "SAMPLENAME"] %in% samples_not_found) {
    #       calculated_index[i, "index"] <- calculated_index[i - 1, "index"]
    #     }
    #   }
    # }
    ylabel <- "cumulative barcode count"
  }
  
  if (plot_total_cum == T){
    
    for (i in 1:length(unique(temp_subset_coldata[[plot_over]]))){
      temp_subset_coldata <- temp_subset_coldata %>% 
        dplyr::add_row(SAMPLENAME = paste("total_clonal_count", as.character(unique(temp_subset_coldata[[plot_over]])[i]), sep="_"))    
      temp_subset_coldata[[plot_over]][dim(temp_subset_coldata)[1]] <- unique(temp_subset_coldata[[plot_over]])[i]
      temp_subset_coldata[[group_by]][dim(temp_subset_coldata)[1]] <-"total"
    }
    plotting_data <- temp_subset_coldata %>%
      dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME)) %>%
      dplyr::left_join(calculated_index, by = "SAMPLENAME")
  }else{
    # merge measures with colData
    plotting_data <- temp_subset_coldata %>%
      dplyr::mutate(SAMPLENAME = as.character(.data$SAMPLENAME)) %>%
      dplyr::left_join(calculated_index, by = "SAMPLENAME")
  }
  
  

  
  # Make sure plot over is a factor if not numeric or specified to not keep numeric.
  if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric) {
  } else if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric == FALSE) {
    plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = unique(plot_over_display_choices))
  } else {
    plotting_data[[plot_over]] <- factor(plotting_data[[plot_over]], levels = levels(plot_over_display_choices))
  }
  
  if (return_table) {
    return(calculated_index)
  }
  
  plotting_data$x_value <- plotting_data[[plot_over]]
  plotting_data$group_by <- plotting_data[[group_by]]
  
  # Create ggplot
  g <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = .data$x_value, y = .data$index, group = .data$group_by, colour = .data$group_by)) +
    ggplot2::geom_line(size = line_size) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::labs(x = plot_over, col = group_by, y = ylabel) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = text_size)) +
    ggplot2::ggtitle(your_title)
  
  if (is.numeric(temp_subset_coldata[[plot_over]]) & keep_numeric) {
    g + ggplot2::scale_x_continuous(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  } else {
    g + ggplot2::scale_x_discrete(paste0(plot_over), breaks = plot_over_display_choices, labels = plot_over_display_choices)
  }
  
  

}