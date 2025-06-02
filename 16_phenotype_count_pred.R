# R script to count phenotype for each blood system
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# load utility functions and summary table
source("./R/utilityfun.R")

superpopulation_list <- c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")

# path to file with sample ID and superpopulation label
PATH_population_label <- "./biobank/1KG/igsr-1000_genomes_phase_3_release.tsv"
df_population_label <- read.table(PATH_population_label, header = TRUE, sep = "\t")
df_population_label$ID <- rownames(df_population_label)
PREDICT_PATH <- "./R/results/predict"
pred_files <- list.files(PREDICT_PATH, pattern = "final_summary_table\\.Rdata$", full.names = TRUE)
output_file <- "./R/results/Phenotype_count_pred.txt"

sink(output_file)
# Count phenotype for each superpopulation
for (file in pred_files) {
    load(file)
    print(file)
    for (superpopulation in superpopulation_list) {
        df_sum_table_subset <- split_by_superpopulation(df_sum_table, df_population_label, superpopulation)
        print(superpopulation)
        count <- table(df_sum_table_subset)
        if (superpopulation == "ALL") {
            print(nrow(df_sum_table_subset))
            cat(paste(names(count), "|------", round(count/nrow(df_sum_table_subset) * 100, 2), "\n"))  
        }
        else {
            print(length(df_sum_table_subset))
            cat(paste(names(count), "|------", round(count/length(df_sum_table_subset) * 100, 2), "\n"))
        }
    }
}

sink()