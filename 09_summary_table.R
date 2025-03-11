# R script to generate summary table of potential phenotypes for each sample in each blood group
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Load necessary libraries

# Set paths
PATH_DATA_TABLE <- "./R/Blood Allele Table.tsv"
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


chr_num <- 1
chr_to_skip <-c(5,8,10,13,14,16,20,21)
resume <- TRUE

while (chr_num <= 23) {
  if (chr_num %in% chr_to_skip) {
    chr_num <- chr_num + 1
    next
  }

  # Load existing RData
  load(paste0("./R/results/chr", chr_num, "_inferred_phenotypes.Rdata"))
  
  
  if (resume) {
    TEMP_FILE <- list.files(path = "./R/results/", pattern = "^chr[0-9]+_temp\\.rds$", full.names = TRUE)
    
    # Resume generating the summary table for chr_num if a save file exists from the last R termination
    if (length(TEMP_FILE) > 0) {
      chr_num <- as.numeric(sub(".*chr([0-9]+)_temp\\.rds", "\\1", TEMP_FILE)) + 1
      df_sum_table <- readRDS(TEMP_FILE)
      resume <- FALSE
      next
    }

    # Otherwise, create the structure of summary table
    else {
      df_sum_table <- data.frame(matrix(ncol = length(bg_to_pheno_allele_map), nrow = nrow(df_genotype)))
      colnames(df_sum_table) <- names(bg_to_pheno_allele_map)
      rownames(df_sum_table) <- rownames(df_genotype)
    }
  }  
  
  print(paste0("Assigning phenotype in chromsome ", chr_num, " to summary table"))

  # Create a vector to store unique phenotype for each sample
  unique_pheno_vec <- setNames(rep(NA_character_, nrow(df_genotype)), rownames(df_genotype))
  
  # Iterate through the samples to store their unique phenotypes in unique_pheno_vec
  for (i in seq_len(nrow(df_inferred_result))) {
    unique_pheno <- unique(as.vector(as.matrix(df_inferred_result[i, ])))
    unique_pheno <- unique(unlist(strsplit(unique_pheno, " \\| ")))
    if (length(unique_pheno) > 1) {
      unique_pheno <- unique_pheno[unique_pheno != "NA/NA"]
    }
    unique_pheno_vec[i] <- paste(unique_pheno, collapse = " | ")
  }

  # Assign the unique phenotypes to their respectively blood system for each sample
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
  
  # Remove previous temporary save file
  rds_files <- list.files(path ="./R/results/", pattern = "\\.rds$", full.names = TRUE)
  file.remove(rds_files)

  # Create a temporary save file to resume generating summary table in case R terminates
  TEMP_SAVE_PATH <- paste0("./R/results/chr", chr_num, "_temp.rds")
  saveRDS(df_sum_table, TEMP_SAVE_PATH)

  print(paste0('Chromosome ', chr_num, ' Done'))
  print(paste("Checkpoint saved for chromosome ", chr_num))
  resume <- FALSE

  # Save summary table once all unique phenotypes has been assigned to their respecitve blood system
  if (chr_num == 23) {
    OUTPUT_PATH <- paste0("./R/results/summary_table.tsv")
    write.table(df_sum_table, file = OUTPUT_PATH, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    save.image(paste0("summary_table.Rdata"))
  }
  chr_num <- chr_num + 1
  }