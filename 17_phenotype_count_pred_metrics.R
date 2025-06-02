# R script to generate summary table of inferred phenotype for each sample in each blood group
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)
library(caret)

calculate_metrics <- function(predictions, reference) {
  conf <- confusionMatrix(data = predictions, reference = reference, mode = "everything")
  num_class <- length(levels(predictions))

  # Ensure 'conf$byClass' is treated as a matrix for consistent indexing
  if (num_class == 2) {
    conf_by_class <- as.data.frame(t(conf$byClass))  # Convert named vector to data frame
  } else {
    conf_by_class <- conf$byClass  # Already a matrix for multiclass
  }

  # Get class prevalence (proportion of each class in the data)
  weights <- conf_by_class[, "Prevalence"] / sum(conf_by_class[, "Prevalence"], na.rm = TRUE)

  # Compute weighted macro-averaged metrics
  metrics <- list(
    Accuracy = conf$overall["Accuracy"],
    Kappa = conf$overall["Kappa"],
    Sensitivity = sum(conf_by_class[, "Sensitivity"] * weights, na.rm = TRUE),
    Specificity = sum(conf_by_class[, "Specificity"] * weights, na.rm = TRUE),
    Pos_Pred_Value = sum(conf_by_class[, "Pos Pred Value"] * weights, na.rm = TRUE),
    Neg_Pred_Value = sum(conf_by_class[, "Neg Pred Value"] * weights, na.rm = TRUE),
    Precision = sum(conf_by_class[, "Precision"] * weights, na.rm = TRUE),
    Recall = sum(conf_by_class[, "Recall"] * weights, na.rm = TRUE),
    F1 = sum(conf_by_class[, "F1"] * weights, na.rm = TRUE),
    Prevalence = mean(conf_by_class[, "Prevalence"], na.rm = TRUE),  # Mean, not weighted
    Detection_Rate = sum(conf_by_class[, "Detection Rate"] * weights, na.rm = TRUE),
    Detection_Prevalence = sum(conf_by_class[, "Detection Prevalence"] * weights, na.rm = TRUE),
    Balanced_Accuracy = sum(conf_by_class[, "Balanced Accuracy"] * weights, na.rm = TRUE)
  )

  return(metrics)
}

# Set predicted results file path
PREDICT_PATH <- "./R/results/predict"
pred_files <- list.files(PREDICT_PATH, pattern = "final_summary_table.*\\.Rdata$", full.names = TRUE)

# Load final summary table generated from original df_genotype and convert each columns to factors
# to prepare for calculating classificaiton metrics
load("./R/results/final_summary_table.Rdata")
original_df_obs_table <- df_sum_table
original_df_obs_table[] <- lapply(original_df_obs_table, function(x) as.factor(trimws(as.character(x))))

chr_results <- list()

# Calculate classification metrics for each final summary table generated from using predicted genotype values
for (file in pred_files) {
  df_obs_table <- original_df_obs_table
  print(file)
  load(file)
  chr_num <- sub(".*/chr([0-9]+)_.*", "\\1", file)
  target_column <- sub(".*_([0-9]+)_.*", "\\1", file)
  pop <- sub(".*_([A-Z]+)\\.Rdata$", "\\1", file)
  col_name <- colnames(df_sum_table)
  df_sum_table[[col_name]] <- as.factor(trimws(as.character(df_sum_table[[col_name]])))
  keep_rows <- rownames(df_sum_table)
  df_obs_table <- df_obs_table[rownames(df_obs_table) %in% keep_rows, ]

  common_levels <- union(levels(df_sum_table[[col_name]]), levels(df_obs_table[[col_name]]))
  df_sum_table[[col_name]] <- factor(df_sum_table[[col_name]], levels = common_levels)
  df_obs_table[[col_name]] <- factor(df_obs_table[[col_name]], levels = common_levels)
  print(chr_num)
  print(col_name)
  print(pop)
  metrics <- calculate_metrics(df_sum_table[[col_name]], df_obs_table[[col_name]])

  metrics_df <- data.frame(
      Target_Column = target_column,

      # Performance metrics
      Test_Accuracy = metrics$Accuracy,
      Test_Kappa = metrics$Kappa,
      Test_Sensitivity = metrics$Sensitivity,
      Test_Specificity = metrics$Specificity,
      Test_Pos_Pred_Value = metrics$Pos_Pred_Value,
      Test_Neg_Pred_Value = metrics$Neg_Pred_Value,
      Test_Precision = metrics$Precision,
      Test_Recall = metrics$Recall,
      Test_F1 = metrics$F1,
      Test_Prevalence = metrics$Prevalence,
      Test_Detection_Rate = metrics$Detection_Rate,
      Test_Detection_Prevalence = metrics$Detection_Prevalence,
      Test_Balanced_Accuracy = metrics$Balanced_Accuracy,
      Chromosome = chr_num,
      Superpopulation = pop,
      row.names = NULL
      )

  chr_results[[length(chr_results) + 1]] <- metrics_df
}

chr_results_dt <- data.table::rbindlist(chr_results, use.names = TRUE, fill = TRUE)
save(chr_results_dt, chr_results, file = "./R/results/chr_1_23_Pred_Pheno_Metrics_XG_100_WMA.Rdata")