# R script to generate summary table of potential phenotypes for each sample in each blood group
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# load utility functions
source("./R/utilityfun.R")

# Function to extract chromosome number in Rdata file
extract_chr <- function(filename) {
  match <- regmatches(filename, regexpr("chr[0-9XY]+", filename))  # Handles chrX, chrY too
  return(sub("chr", "", match))
}

# Set paths
PATH_DATA_TABLE <- "./R/Blood_Allele_Table.tsv"
df_database <- read.delim(PATH_DATA_TABLE, header = TRUE)
PATH_DATA_TABLE_SEP <- "./R/Blood_Allele_Table_Separated.tsv"
df_database_sep <- read.delim(PATH_DATA_TABLE_SEP, header = TRUE, fileEncoding = "UTF-8")
PREDICT_PATH <- "./R/results/predict"
pred_files <- list.files(PREDICT_PATH, pattern = "inferred_phenotypes.*\\.Rdata$", full.names = TRUE)
split_files <- split(pred_files, sapply(pred_files, extract_chr))

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

  for (file in split_files[[as.character(chr_num)]]) {
        
        col <- sub(".*chr[0-9]+_([0-9]+)_.*", "\\1", file)
        pop <- sub(".*_([A-Z]+)\\.Rdata$", "\\1", file)

        if (file.exists(paste0("./R/results/predict/chr", chr_num, "_", col, "_summary_table.Rdata"))) {
          next
        }
        
        # Load inferred phenotype results
        load(file)

        # Create the structure of summary table
        df_sum_table <- data.frame(matrix(ncol = length(bg_to_pheno_allele_map), nrow = nrow(temp_df_genotype)))
        colnames(df_sum_table) <- names(bg_to_pheno_allele_map)
        rownames(df_sum_table) <- rownames(temp_df_genotype)

        print(paste0("Assigning phenotypes in Chromosome ", chr_num, " (", col, ") to summary table (", pop, ")"))
  
        df_sum_table <- generate_sum_table(df_sum_table, df_inferred_result, temp_df_genotype, bg_to_pheno_allele_map)

        keep_bs <- df_database_sep %>%
          filter(as.character(`GRCh37.VCF.Position`) == col) %>%
          pull(`Blood.System`) %>%
          dplyr::first()

        df_sum_table <- df_sum_table[, colSums(is.na(df_sum_table)) < nrow(df_sum_table), drop = FALSE]
        df_sum_table <- df_sum_table[, colnames(df_sum_table) %in% keep_bs, drop = FALSE]
        
        save(df_sum_table, file = paste0("./R/results/predict/chr", chr_num, "_", col, "_summary_table.Rdata"))
    }
    print(paste0('Chromosome ', chr_num, ' Done'))
}


