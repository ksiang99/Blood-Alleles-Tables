# R script to infer phenotype based on variant positions
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Load necessary libraries

# Set paths
PATH_DATA_TABLE <- "./R/Blood Allele Table.tsv"
df_database <- read.delim(PATH_DATA_TABLE, header = TRUE)
df_database$Pheno_Allele <- paste(df_database$Phenotype, df_database$Allele.Name, sep = "/")
df_database <- df_database[order(df_database$Chromosome), ]
bg_to_pheno_allele_map <- split(df_database$Pheno_Allele, factor(df_database$Blood.System, levels = unique(df_database$Blood.System)))
values_to_drop = c("PEL", "ABCC1", "ER", "KANNO", "EMM", "CTL2", "MAM")
bg_to_pheno_allele_map <- bg_to_pheno_allele_map[setdiff(names(bg_to_pheno_allele_map), values_to_drop)]

chr_num <- 1
chr_to_skip <-c(5,8,10,13,14,16,20,21)
resume <- TRUE

while (chr_num <= 23) {
  if (chr_num %in% chr_to_skip) {
    chr_num <- chr_num + 1
    next
  }
  
  load(paste0("./R/results/chr", chr_num, "_inferred_phenotypes.Rdata"))
  
  if (resume) {
    TEMP_FILE <- list.files(path = "./R/results/", pattern = "^chr[0-9]+_temp\\.rds$", full.names = TRUE)
    if (length(TEMP_FILE) > 0) {
      chr_num <- as.numeric(sub(".*chr([0-9]+)_temp\\.rds", "\\1", TEMP_FILE)) + 1
      df_sum_table <- readRDS(TEMP_FILE)
      resume <- FALSE
      next
    }

    else {
      df_sum_table <- data.frame(matrix(ncol = length(bg_to_pheno_allele_map), nrow = nrow(df_genotype)))
      colnames(df_sum_table) <- names(bg_to_pheno_allele_map)
      rownames(df_sum_table) <- rownames(df_genotype)
    }
  }  
  
  print(paste0("Assigning phenotype in chromsome ", chr_num, " to summary table"))
  unique_pheno_vec <- setNames(rep(NA_character_, nrow(df_genotype)), rownames(df_genotype))
  
  for (i in seq_len(nrow(df_inferred_result))) {
    unique_pheno <- unique(as.vector(as.matrix(df_inferred_result[i, ])))
    unique_pheno <- unique(unlist(strsplit(unique_pheno, " \\| ")))
    if (length(unique_pheno) > 1) {
      unique_pheno <- unique_pheno[unique_pheno != "NA/NA"]
    }
    unique_pheno_vec[i] <- paste(unique_pheno, collapse = " | ")
  }

  for (sample in names(unique_pheno_vec)) {
      phenotype_string <- unique_pheno_vec[[sample]]
      phenotypes <- unlist(strsplit(phenotype_string, " \\| "))

      for (pheno in phenotypes) {
          matched_system <- names(bg_to_pheno_allele_map)[sapply(bg_to_pheno_allele_map, function(x) any(pheno == unlist(x)))]

          for (system in matched_system) {
            df_sum_table[sample, system] <- ifelse(is.na(df_sum_table[sample, system]), pheno, 
                                                        paste(df_sum_table[sample, system], pheno, sep=" | "))
          }
      }
  }
  
  rds_files <- list.files(path ="./R/results/", pattern = "\\.rds$", full.names = TRUE)
  file.remove(rds_files)
  TEMP_SAVE_PATH <- paste0("./R/results/chr", chr_num, "_temp.rds")
  saveRDS(df_sum_table, TEMP_SAVE_PATH)
  print(paste0('Chromosome ', chr_num, ' Done'))
  print(paste("Checkpoint saved for chromosome ", chr_num))
  resume <- FALSE

  if (chr_num == 23) {
    OUTPUT_PATH <- paste0("./R/results/summary_table.tsv")
    write.table(df_sum_table, file = OUTPUT_PATH, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    save.image(paste0("summary_table.Rdata"))
  }
  chr_num <- chr_num + 1
  }