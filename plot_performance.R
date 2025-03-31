# script to plot performance metric of naive bayes classifier in predicting blood type causal variants

# MUST BE IN the directory Blood-type-GWAS and not any other subdirectory when launching R to activate Renv
rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(purrr)
library(tidyr)
library(data.table)
library(dplyr)
library(reshape2)
#library(modules)

# avoid using attach() as it will pollute the environment attach("R")

# set working directory on HPC or codespace
setwd("/home/svu/e0726996/Blood-type-GWAS/")

# load utility functions
# source("../../../R/utilityfun.R")

# --- functions ----
# function to add gene column to chr_results_dt
# Description: This function adds a gene column to the input data frame based on the Target_Column values.
# Inputs:
#   - data: A data frame containing the results of the analysis.
# Outputs:
#   - A data frame with additional columns 'Gene' and 'Gene_Coordinate'.
add_gene_column <- function(data){
  PATH_TABLE <- "./R/erythrogene_coordinate_fixed_with_chromosome.tsv"
  # read in table and remove trailing white space
  df_allele <- read.delim(PATH_TABLE, header = TRUE)

  # Trim trailing whitespace from all character columns
  df_allele[] <- lapply(df_allele, function(x) {
    if (is.character(x)) {
      trimws(x)  # Remove leading and trailing whitespace
    } else {
      x  # Return non-character columns unchanged
    }
  })

  data$Gene <- df_allele$Gene[match(data$Target_Column, df_allele$GRCh37_Start)] %>% trimws() # should check for same chromosome but unlikely 
  data$HGVS <- df_allele$Nucleotide.Change[match(data$Target_Column, df_allele$GRCh37_Start)] %>% trimws() # HGVS notation
  data <- as.data.frame(data)
# create new column by concatenating Gene and Target_Column
  data$Gene_Coordinate <- paste(data$Gene, data$Target_Column, sep = "_")

  return(data)
}

#---- PLOTS ENTIRE DATAFRAME ----

# chromosome by population by metric by gene
# function to plot chromosome by population by metric by gene
# Description: This function generates and saves a plot for a specific chromosome, population, and metric by gene.
# Inputs:
#   - df_all_res: A data frame containing the results of the analysis.
#   - chr_num: The chromosome number to plot.
#   - metric: The performance metric to plot (default is "Test_F1").
# Outputs:
#   - A ggplot object representing the plot.
plot_chromosome_population_metric_gene <- function(df_all_res, chr_num, metric = "Test_F1"){
  df <- df_all_res[df_all_res$Chromosome == chr_num, ]

  p <- ggplot(df, aes(x = HGVS, y = !!sym(metric), color = MAF, shape = Superpopulation)) +
    geom_point(size = 3) + 
    scale_shape_manual(values = c("AFR" = 15, "AMR" = 16, "EAS" = 17, "EUR" = 18, "SAS" = 9, "ALL" = 25)) + 
    scale_color_gradient(low = "powderblue", high = "gold") +
    labs(
      title = paste0("Chr ", chr_num, " ", metrics_label[metric]),
      x = "Variant",
      y = metrics_label[metric],
      color = "MAF",
      shape = "Population"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12)
    ) +
    ylim(c(0, 1))
  
  # return the plot
  ggsave(paste0("./R/results/Plots/NB_CV/chr", chr_num, "_", metric,".png") , plot = p, width = 8, height = 6, dpi = 300)
  return(p)
}


# Function to plot metrics by chromosome
# Description: This function generates and saves facet plots for all specified metrics by chromosome.
# Inputs:
#   - df_all_res: A data frame containing the results of the analysis.
#   - metrics: A vector of performance metrics to plot (default includes various metrics).
# Outputs:
#   - None (plots are saved to files).
plot_all_metric_by_chromosome <- function(df_all_res ,metrics = c("Test_Prevalence", 'Test_Accuracy', 'Test_Balanced_Accuracy', 'Test_Sensitivity', 'Test_Specificity', 'Test_Pos_Pred_Value', 'Test_Neg_Pred_Value', 'Test_Precision', 'Test_F1', "AUC" )) {
  # Chromosomes to skip
  chr_to_skip <- c(5, 8, 10, 13, 14, 16, 20, 21)

  for (chr_num in 1:23) {
    if (chr_num %in% chr_to_skip) {
      next
    }
    
    # Subset df_all_res by chromosome
    chr_results_dt_subset <- df_all_res[df_all_res$Chromosome == chr_num, ]

    # Convert 'Target_Column' to factor if it's not already
    chr_results_dt_subset$Target_Column <- as.factor(chr_results_dt_subset$Target_Column)

    # Specify the metrics columns to reshape (modify if needed)
    # metrics <- colnames(df_all_res)[!colnames(df_all_res) %in% c("Chromosome", "Gene", "Gene_Coordinate", "Target_Column", "MAF", "Superpopulation")]

    # Reshape the data into long format
    chr_results_dt_long <- chr_results_dt_subset %>%
      pivot_longer(cols = all_of(metrics), 
                    names_to = "Metric", 
                    values_to = "Value")

    # Generate the facet plot
    fp <- ggplot(chr_results_dt_long, aes(x = HGVS, y = Value, color = MAF, shape = Superpopulation)) +
      geom_point(size = 3) +
      scale_shape_manual(values = c("AFR" = 15, "AMR" = 16, "EAS" = 17, "EUR" = 18, "SAS" = 9, "ALL" = 25)) + 
      scale_color_gradient(low = "powderblue", high = "gold") +
      facet_wrap(~ Metric, nrow = 5, ncol = 2, scales = "free_y") +
      labs(title = paste0("Chr ", chr_num),
        x = "Variant",
        y = "Performance Metric",
        color = "MAF",
        shape = "Population") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12)) + 
      ylim(c(0, 1)) 


    # Save the plot for each chromosome
    ggsave(paste0("./R/results/Plots/NB_CV/chr", chr_num, "_all_performance_metrics.png"), plot = fp, width = 16, height = 12, dpi = 300)
  }
}

# Function to plot cumulative count by metrics 
# Description: This function generates and saves cumulative count plots for all specified metrics by population.
# Inputs:
#   - df_all_res: A data frame containing the results of the analysis.
#   - metrics: A vector of performance metrics to plot (default includes various metrics).
# Outputs:
#   - A list of ggplot objects representing the cumulative count plots.
plot_cumulative_count <- function(df_all_res, metrics = c("Test_Prevalence", 'Test_Accuracy', 'Test_Balanced_Accuracy', 'Test_Sensitivity', 'Test_Specificity', 'Test_Pos_Pred_Value', 'Test_Neg_Pred_Value', 'Test_Precision', 'Test_F1', "AUC" )){
  # Initialize plot list
  plots <- list()
  
  for (metric in metrics){
    # Calculate cumulative counts grouped by superpopulation
    cum_counts <- df_all_res %>%
      group_by(Superpopulation, !!sym(metric)) %>%
      summarise(Count = n()) %>%
      group_by(Superpopulation) %>%
      mutate(Cumulative_Counts = cumsum(Count))

    # Plot cumulative counts for all populations
    cp <- ggplot(cum_counts, aes(x = !!sym(metric), y = Cumulative_Counts, color = Superpopulation)) + 
      geom_line(size = 1) + 
      geom_point(size = 2) + 
      labs(title = paste0("Cumulative Count for ", metric, " by Population"), 
          x = metrics_label[[metric]], 
          y = "Number of Variants") + 
      theme(legend.position = "bottom") +
      xlim(c(0, 1)) 

    
    # Save plot
    ggsave(file.path("./R/results/Plots/NB_CV/", paste0("cumulative_count_", metric, ".png")), plot = cp, width = 16, height = 12, dpi = 300)
    
    # Store plot in list
    plots[[metric]] <- cp
  }
  
  # Return plot list
  return(plots)
}

# function to plot number of variants more than threshold for metric
# Description: This function generates and saves plots showing the number of variants exceeding a threshold for each specified metric.
# Inputs:
#   - df_all_res: A data frame containing the results of the analysis.
#   - metrics: A vector of performance metrics to plot (default includes various metrics).
# Outputs:
#   - A list of ggplot objects representing the threshold plots.
plot_threshold_metric <- function(df_all_res, metrics = c("Test_Prevalence", 'Test_Accuracy', 'Test_Balanced_Accuracy', 'Test_Sensitivity', 'Test_Specificity', 'Test_Pos_Pred_Value', 'Test_Neg_Pred_Value', 'Test_Precision', 'Test_F1', "AUC" )){
  
  # Initialize plot list
  plots <- list()

  # Create a vector of threshold values
  threshold_values <- seq(0,1, by = 0.01)
  for (metric in metrics){

    ttt <- df_all_res[,c(metric,"Superpopulation")] # subset data

    # make data frame with columns for population, threshold values, count of variants more than threshold value
    plot_df <- map_dfr(threshold_values, function(threshold) {
      ttt %>%
        group_by(Superpopulation) %>%
        summarise(Counts = sum(!!sym(metric) >= threshold)) %>%
        mutate(Threshold = threshold)
    }, .id = NULL)

    cp <- ggplot(plot_df, aes(x = Threshold, y = Counts, color = Superpopulation)) + 
      geom_line(size = 1) + 
      geom_point(size = 2) + 
      labs(title = paste0("Number of Variants by ", metrics_label[[metric]]), 
        x = metrics_label[[metric]], 
        y = "Number of Variants (MAF > 0.01)",
        color = "Superpopulation") + 
      theme(legend.position = "bottom",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12)) +
      xlim(c(0, 1)) +
      ylim(c(0, max(plot_df$Counts)))
    
    # Save plot
    ggsave(file.path("./R/results/Plots/NB_CV/", paste0("number_of_variants_by_threshold_", metric, ".png")), plot = cp, dpi = 300)
    
    # Store plot in list
    plots[[metric]] <- cp
  }

  # Return plot list
  return(plots)
}

# function to plot MAF by prevalence by metric threshold
# Description: This function generates and saves a scatter plot of MAF by prevalence for variants exceeding a specified metric threshold.
# Inputs:
#   - df_all_res: A data frame containing the results of the analysis.
#   - metric_threshold: The threshold value for the metric (default is 0.9).
#   - metric: The performance metric to filter by (default is "Test_F1").
# Outputs:
#   - A filtered data frame containing variants that meet the threshold criteria.
plot_MAF_by_prevalence_by_metric_threshold <- function(df_all_res, metric_threshold =0.9, metric = "Test_F1"){
  # filter by metric threshold
  df_all_res <- df_all_res %>% filter(!!sym(metric) >= metric_threshold, Superpopulation == "ALL")
  # Scatter plot with labels and color by Superpopulation
  sp <- ggplot(df_all_res, aes(x = Train_Prevalence, y = !!sym(metric), color = MAF)) +
        geom_point(aes(color = MAF), size = 3) +  # Color points by MAF
        geom_text(aes(label = Gene_Coordinate), 
                color = "black", vjust = -0.05, hjust = 0.05, size = 6) +  # Black text labels with increased size
        labs(title = paste0("Variants with ", metric , " over ", metric_threshold, " by Reference Allele Frequency"), 
          y = metrics_label[[metric]],
          color = "MAF",
          x = "Reference Allele Frequency") + 
        scale_color_gradient(low = "powderblue", high = "gold") +
        theme(#legend.position = "bottom",
              # axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 20, face = "bold"),
              axis.title = element_text(size = 18),  # Increased axis title size
              axis.text = element_text(size = 16),  # Increased axis tick labels size
              legend.title = element_text(size = 16),  # Increased legend title size
              legend.text = element_text(size = 14))  # Increased legend text size
  
  ggsave(file.path("./R/results/Plots/NB_CV/", paste0("Prevalence_by_", metric, "_threshold_", metric_threshold, ".png")), plot = sp, width = 16, height = 12, dpi = 300)
  return(df_all_res)
}

# ---- MAIN SCRIPT ----

PATH_RESULT <- paste0("./R/results/chr_1_23_NB_CV.Rdata")
load(PATH_RESULT)

metrics <- c("Test_Prevalence", 'Test_Accuracy', 'Test_Balanced_Accuracy', 'Test_Sensitivity', 
              'Test_Specificity', 'Test_Pos_Pred_Value', 
              'Test_Neg_Pred_Value', 'Test_Precision', 'Test_F1', "AUC" )

metrics_label <- c(Test_Prevalence = "Prevalence", 
                    Test_Accuracy = "Accuracy", 
                    Test_Balanced_Accuracy = "Balanced Accuracy", 
                    Test_Sensitivity = "Sensitivity", 
                    Test_Specificity = "Specificity", 
                    Test_Pos_Pred_Value = "Positive Predictive Value", 
                    Test_Neg_Pred_Value = "Negative Predictive Value", 
                    Test_Precision = "Precision", 
                    Test_F1 = "F1", 
                    AUC = "AUC"
                  )

# add gene column to chr_results_dt
df_all_res <- as.data.frame(add_gene_column(chr_results_dt))
dim(chr_results_dt)
dim(df_all_res)
# str(df_all_res) # not all 43 genes in erythrogene table present, double check ld.ld result table
df_all_res$Target_Column <- as.numeric(df_all_res$Target_Column) 

# replace NA by 0 in df_all_res
df_all_res[is.na(df_all_res)] <- 0 

chr_to_skip <- c(5, 8, 10, 13, 14, 16, 20, 21)

for (chr_num in 1:23) {
  if (chr_num %in% chr_to_skip) {
    next
  }
  for (metric in metrics) {
    plot_chromosome_population_metric_gene(df_all_res, chr_num, metric)
  }
}


plot_all_metric_by_chromosome(df_all_res)

#cpp<-plot_cumulative_count(df_all_res)
plots <- plot_threshold_metric(df_all_res)

# list of variants over threshold 
for (mm in metrics){
  if (mm != "Test_Prevalence") {
    plot_MAF_by_prevalence_by_metric_threshold(df_all_res, metric = mm)
  }
}           
# plot_MAF_by_prevalence_by_metric_threshold(df_all_res)
