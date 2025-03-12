# R script to count phenotype for each blood system
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Load necessary libraries
library(dplyr)

# load utility functions and summary table
source("./R/utilityfun.R")
load(paste0("./R/results/final_summary_table.Rdata"))

superpopulation_list <- c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")

# path to file with sample ID and superpopulation label
PATH_population_label <- "./biobank/1KG/igsr-1000_genomes_phase_3_release.tsv"
df_population_label <- read.table(PATH_population_label, header = TRUE, sep = "\t")
df_population_label$ID <- rownames(df_population_label)

# Count phenotype for each superpopulation
for (superpopulation in superpopulation_list) {
    df_sum_table_subset <- split_by_superpopulation(df_sum_table, df_population_label, superpopulation)
    print(superpopulation)
    count <- table(df_sum_table_subset$P1PK)
    cat(paste(names(count), "::::::", count, "\n"))
    
    # count <- sum(grepl("Lu_thirteen:", df_sum_table_subset$LU))
    # print(count/nrow(df_sum_table_subset)*100)
    
}
