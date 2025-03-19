# R script to generate summary table of potential phenotypes for each sample in each blood group
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# load utility functions
source("./R/utilityfun.R")

# Set paths
PATH_DATA_TABLE <- "./R/Blood_Allele_Table.tsv"
df_database <- read.delim(PATH_DATA_TABLE, header = TRUE)

# Create a column that joins the phenotype and allele name Eg. Phenotype/allele name
df_database$Pheno_Allele <- paste(df_database$Phenotype, df_database$Allele.Name, sep = "/")

# Order the blood system in ascending chromsome number
df_database <- df_database[order(df_database$Chromosome), ]

# Create a list mapping of blood system to Pheno_Allele (Key to values)
bg_to_pheno_allele_map <- split(df_database$Pheno_Allele, factor(df_database$Blood.System, levels = unique(df_database$Blood.System)))

# Drop blood system not in Erythrogene
values_to_drop = c("PEL", "ABCC1", "ER", "KANNO", "EMM", "CTL2", "MAM", "SID")
bg_to_pheno_allele_map <- bg_to_pheno_allele_map[setdiff(names(bg_to_pheno_allele_map), values_to_drop)]

chr_to_skip <-c(5,8,10,13,14,16,20,21)

for (chr_num in 1:23) {
  if (chr_num %in% chr_to_skip) {
    next
  }

  for (col_name in seq_long(target_columns)) {
    for (fold in 1:5) {
        
        # Load inferred phenotype results
        load(paste0("./R/results/chr", chr_num, "_", col_name, "_fold", fold, "_inferred_phenotypes.Rdata"))

        # Create the structure of summary table
        df_sum_table <- data.frame(matrix(ncol = length(bg_to_pheno_allele_map), nrow = nrow(temp_df_genotype)))
        colnames(df_sum_table) <- names(bg_to_pheno_allele_map)
        rownames(df_sum_table) <- rownames(temp_df_genotype)

        print(paste0("Assigning phenotypes in Chromosome ", chr_num, " for ", fold, " (", col_name, ") to summary table"))
  
        df_sum_table <- generate_sum_table(df_sum_table, df_inferred_result, temp_df_genotype, bg_to_pheno_allele_map)
        df_sum_table <- df_sum_table[, colSums(is.na(df_sum_table)) < nrow(df_sum_table)]
        
        OUTPUT_PATH <- paste0("./R/results/", colname, "_", "fold", fold, "_summary_table.tsv")
        write.table(df_sum_table, file = OUTPUT_PATH, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
        save(df_sum_table, file = paste0("./R/results/", colname, "_", fold, "_summary_table.Rdata"))
    }
    print(paste0('Chromosome ', chr_num, ' Done'))
  }
}

