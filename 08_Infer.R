# R script to infer phenotype based on variant positions
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# load utility functions
source("./R/utilityfun.R")

# Set the path to the data table
dir.create("./R/results", showWarnings = FALSE)
PATH_DATA_TABLE <- "./R/Blood_Allele_Table_Separated.tsv"
df_database <- read.delim(PATH_DATA_TABLE, header = TRUE)
df_database$GRCh37.VCF.Position <- as.character(df_database$GRCh37.VCF.Position)
df_database$Gene <- as.character(df_database$Gene)

chr_to_skip <-c(5,8,10,13,14,16,20,21)

for (chr_num in 1:23) {
  if (chr_num %in% chr_to_skip) {
    next
  }

  if (file.exists(paste0("./R/results/chr", chr_num, "_inferred_phenotypes.Rdata"))) {
    next
  }

  print(paste0("Processing Chromosome ", chr_num))

  # Create a vector of unique variant positions of current chromosome number found in df_database
  pos_to_keep <- df_database %>%
    filter(Chromosome == chr_num) %>%
    pull(GRCh37.VCF.Position) %>%
    unique()

  # Load df_genotype data
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

  # Get potential phenotypes for each variant position
  pheno_results <- get_phenotypes(variant_pos, chr_num, df_database)
  print(paste0("Assigning phenotype to Chromosome ", chr_num))

  # Assign phenotypes to positions in a dataframe
  df_inferred_result <- pos_pheno_df(df_genotype, pheno_results)

  OUTPUT_PATH <- paste0("./R/results/chr", chr_num, "_inferred_phenotypes.tsv")
  write.table(df_inferred_result, file = OUTPUT_PATH, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  save(df_genotype, df_inferred_result, file = paste0("./R/results/chr", chr_num, "_inferred_phenotypes.Rdata"))
  print(paste0('Chromosome ', chr_num, ' done'))
}