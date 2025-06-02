rm(list = ls())

# Set random seed
set.seed(0)

# set working directory on HPC or codespace
setwd("/home/svu/e0726996/Blood-type-GWAS/")

# load utility functions
source("./R/utilityfun.R")

options(warn = 1)

# --- Main script ---
chr_results <- list()
RF_results <- list()
chr_to_skip <- c(5,8,10,13,14,16,20,21)

# path to file with sample ID and superpopulation label
PATH_population_label <- "./biobank/1KG/igsr-1000_genomes_phase_3_release.tsv"

df_population_label <- read.table(PATH_population_label, header = TRUE, sep = "\t")
df_population_label$ID <- rownames(df_population_label)

# ALL to stand in for using whole df_genotype
superpopulation_list <- c("ALL", "EUR", "EAS", "AMR", "SAS","AFR")

chr_num <- 23

for (chr_num in 1:23) {
  if (chr_num %in% chr_to_skip) {
    next
  }
  
  print(chr_num)
  path_rdata <- paste0(getwd(), "/biobank/1KG/phase3_grch37/predict/results/chr", chr_num, "_df_genotype.Rdata")
  if (file.exists(path_rdata)) {
    load(path_rdata)
  } else {
    warning(paste("File not found:", path_rdata))
    next
  }
  
  cv_MAF <- 1 - sapply(causal_variants, function(x) sum(df_genotype[[as.character(x)]] == 0) / length(df_genotype[[as.character(x)]]))
  cat("\nNumber of variants MAF > 0.01: ", sum(cv_MAF > 0.01), "\n")

  for (superpopulation in superpopulation_list) {
    df_genotype_subset <- split_by_superpopulation(df_genotype, df_population_label, superpopulation)
    print(dim(df_genotype_subset))
    print(superpopulation)

    for (i in seq_along(causal_variants)) {
      target_column <- as.character(causal_variants[i])
      masked_columns <- as.character(causal_variants)
      
      if (cv_MAF[i] < 0.01) {
        next
      }

      tryCatch({
        cat("target_column: ", target_column, "MAF ", 
            1 - sum(df_genotype_subset[[target_column]] == 0) / length(df_genotype[[target_column]]), "\n")
        
        if (length(unique(df_genotype_subset[[target_column]])) < 2) {
          warning(paste("Not enough data points in target column:", target_column))
          next
        }

        df_results <- train_and_evaluate_model_XG(df_genotype_subset, target_column, masked_columns)
        if (is.null(df_results)) {
          warning(paste("Model training failed for target column:", target_column))
          next
        }

        df_results$Chromosome <- chr_num
        df_results$Superpopulation <- superpopulation
        chr_results[[length(chr_results) + 1]] <- df_results
      }, error = function(e) {
        warning(paste("Error processing target column:", target_column, "in superpopulation:", superpopulation, " - ", e$message))
      })
    }
  }
}

chr_results_dt <- data.table::rbindlist(chr_results, use.names = TRUE, fill = TRUE)
save(chr_results_dt, chr_results, file = "./R/results/predict/chr_1_23_XG_100_ONE_WMA.Rdata")

dim(chr_results_dt)