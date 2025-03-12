rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Load necessary libraries
library(dplyr)

# Set the path to the data table
PATH_DATA_TABLE <- "./R/Blood Allele Table (Separated).tsv"
df_database <- read.delim(PATH_DATA_TABLE, header = TRUE)

# Initiate variables
chr_to_skip <-c(5,8,10,13,14,16,20,21)
remove_bs <- c("KANNO", "MAM", "CLT2", "SID", "ER", "ABCC1", "PEL", "EMM")
all_chr_not_positions <- c()

# Filter out blood system not in Erythrogene
df_database <- df_database %>%
    filter(!(Blood.System %in% remove_bs))

for (chr_num in 1:23) {
  if (chr_num %in% chr_to_skip) {
    next
  }

  # Load existing RData
  load(paste0("./biobank/1KG/phase3_grch37/predict/results/chr", chr_num, "_df_genotype.Rdata"))
  df_genotype <- df_genotype %>% dplyr::select(-contains("."))
  
  # Get ISBT variants in current chromosome number. Exclude CD59.
  df_chr <- df_database %>%
    filter(Chromosome == chr_num & Phenotype != "-" & GRCh37.VCF.Position != "-" & 
        Phenotype != "CD59") %>%
        distinct(`Nucleotide.Change`, .keep_all = TRUE)
  
  # Get positions not found in df_genotype
  positions <- df_chr$`GRCh37.VCF.Position`
  not_positions <- positions[!(positions %in% colnames(df_genotype))]
  chr_not_positions <- paste0(chr_num, ":", not_positions)
  all_chr_not_positions <- c(all_chr_not_positions, chr_not_positions)

}

write.table(all_chr_not_positions, "./R/missing_positions.tsv", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

print("Done")