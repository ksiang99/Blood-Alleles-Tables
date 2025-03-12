rm(list = ls())

library(dplyr)
library(data.table)
library(tidyverse)


# ---- -----
PATH_DATA_TABLE <- "Blood Allele Table.tsv"
df_database <- as.data.table(read.delim(PATH_DATA_TABLE))
df_database$Pheno_Allele <- paste(df_database$Phenotype, df_database$Allele.Name, sep = "/")
df_database <- df_database[order(df_database$Chromosome), ]

# group by blood group system and drop uncommon blood groups
bg_to_pheno_allele_map <- split(df_database$Pheno_Allele, factor(df_database$Blood.System, levels = unique(df_database$Blood.System)))
values_to_drop = c("PEL", "ABCC1", "ER", "KANNO", "EMM", "CTL2", "MAM")
bg_to_pheno_allele_map <- bg_to_pheno_allele_map[setdiff(names(bg_to_pheno_allele_map), values_to_drop)]





# ---- check missing positions in df_genotype ----
# load data
load("Results/chr9_inferred_phenotypes.Rdata")
df_genotype_9 <- as.data.table(df_genotype)
chr9_pos <- sort(as.numeric(colnames(df_genotype_9)))
136131594 %in% chr9_pos # FALSE
136131595 %in% chr9_pos # FALSE
136135237 %in% chr9_pos # TRUE
136135238 %in% chr9_pos # TRUE

# convert factor columns to numeric and take columnwise sum of df_genotype_9
df_genotype_9_sum <- df_genotype_9[, lapply(.SD, function(x) if(is.factor(x)) as.numeric(as.character(x)) else x)]
df_genotype_9_sum <- colSums(df_genotype_9_sum, na.rm = TRUE)
df_genotype_9_sum["136135236"] # matches IVG 
df_genotype_9_sum["136131616"] # 



# convert to tibble
df_database <- as_tibble(df_database)



# filter and keep only chromosome 23
df_database_23 <- df_database %>% filter(Chromosome == 23, Gene %in% c("XG", "XK"))
