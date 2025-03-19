# R script to infer phenotype based on variant positions
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Load necessary libraries
library(dplyr)
library(data.table)

# load utility functions
source("./R/utilityfun.R")

# Set the path to the data table
dir.create("./R/results/predict", showWarnings = FALSE)
PATH_DATA_TABLE <- "./R/Blood_Allele_Table_Separated.tsv"
df_database <- read.delim(PATH_DATA_TABLE, header = TRUE)
df_database$GRCh37.VCF.Position <- as.character(df_database$GRCh37.VCF.Position)
df_database$Gene <- as.character(df_database$Gene)

# Load predicted results from Machine Learning
load(paste0("./R/results/chr_1_23_NB_CV.Rdata"))
chr_results_dt[, Pred_Row_Index := lapply(strsplit(Pred_Row_Index, ",\\s*"), as.numeric)]
chr_results_dt[, Pred_Class := lapply(strsplit(Pred_Class, ",\\s*"), as.factor)]

chr_to_skip <-c(5,8,10,13,14,16,20,21)

for (chr_num in 1:23) {
  if (chr_num %in% chr_to_skip) {
    next
  }

  print(paste0("Processing Chromosome ", chr_num))

  # Create a temporary data.table to store rows with current chromosome number
  temp_chr_results_dt <- chr_results_dt[Chromosome == chr_num]

  # Store unique variant positions of current chromosome number found in df_database in a vector
  pos_to_keep <- df_database %>%
    filter(Chromosome == chr_num) %>%
    pull(GRCh37.VCF.Position) %>%
    unique()

  # Load df_genotype
  load(paste0("./biobank/1KG/phase3_grch37/predict/results/chr", chr_num, "_df_genotype.Rdata"))
  df_genotype <- df_genotype %>% dplyr::select(-contains("."))

  # Keep only positions in cols_to_keep to decrease run time
  df_genotype <- df_genotype %>%
    dplyr::select(intersect(names(df_genotype), pos_to_keep))

  # Change genotype values from factor to character type
  df_genotype <- df_genotype %>%
    mutate(across(everything(), as.character))
  
  # Get variant positions from df_genotype
  variant_pos <- as.character(colnames(df_genotype))
  
  for (i in seq_len(nrow(temp_chr_results_dt))) {
    temp_df_genotype <- df_genotype
    
    # Update genotype values of target column in df_genotype with results from Machine Learning
    col <- temp_chr_results_dt[i, Target_Column]
    fold_num <- temp_chr_results_dt[i, Fold]
    row_idx <- unlist(temp_chr_results_dt[i, Pred_Row_Index])
    val <- temp_chr_results_dt[i, Pred_Class]
    temp_df_genotype[row_idx, col] <- val

    # Keep only the rows in temp_df_genotype where genotype values are updated
    temp_df_genotype <- temp_df_genotype[row_idx, ]

    # Get potential phenotypes for each variant position
    pheno_results <- get_phenotypes(variant_pos, chr_num, df_database)
    print(paste0("Assigning phenotypes to Chromosome ", chr_num, " for ", fold_num, " (", col, ")"))

    # Assign phenotypes to positions in a dataframe
    df_inferred_result <- pos_pheno_df(temp_df_genotype, pheno_results)

    OUTPUT_PATH <- paste0("./R/results/predict/chr", chr_num, "_", col, "_", fold_num, "_inferred_phenotypes.tsv")
    write.table(df_inferred_result, file = OUTPUT_PATH, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    save(temp_df_genotype, df_inferred_result, file = paste0("./R/results/predict/chr", chr_num, "_", col, "_", fold_num, "_inferred_phenotypes.Rdata"))
    print(paste0('Chromosome ', chr_num, " for fold ", fold, " (", col_name, ") done"))
  }
}