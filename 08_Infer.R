# R script to infer phenotype based on variant positions
rm(list = ls())

setwd("/home/svu/e0726996/Blood-type-GWAS/")

# Set random seed
set.seed(0)

# Load necessary libraries
library(dplyr)

# Function to assign NA to vector when no ref/alt alleles and phenotypes are found
null_na_fill <- function(x) {
  if (is.null(x)) {
    return(NA)
  }
  return(x)
}

# Set the path to the data table
dir.create("./R/results", showWarnings = FALSE)
PATH_DATA_TABLE <- "./R/Blood Allele Table (Separated).tsv"
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

  print(paste0("Processing chromosome ", chr_num))

  # Load existing RData
  load(paste0("./biobank/1KG/phase3_grch37/predict/results/chr", chr_num, "_df_genotype.Rdata"))
  df_genotype <- df_genotype %>% dplyr::select(-contains("."))
  mat_genotype <- as.matrix(df_genotype)
  
  # Get variant positions from df_genotype
  variant_pos <- as.character(colnames(df_genotype))
  
  # Create vectors to store reference/alterative phenotype/allele names
  alt_phenotypes <- vector("list", length = length(variant_pos))
  alt_allele_names <- vector("list", length = length(variant_pos))
  ref_phenotypes <- vector("list", length = length(variant_pos))
  ref_allele_names <- vector("list", length = length(variant_pos))
  coord_combi_lst <- vector("list", length = length(variant_pos))
  alt_combi_lst <- vector("list", length = length(variant_pos))

  # Loop through each variant position
  for (i in seq_along(variant_pos)) {

    # Check if the variant position is found in df_database.
    # If so, retrieve all phenotypes which are associated to the variant position.
    matched_phenotypes <- df_database %>%
      filter(GRCh37.VCF.Position == variant_pos[i] &
            Chromosome == chr_num) %>%
      dplyr::select(Phenotype, Allele.Name, Gene, Chromosome, GRCh37.VCF.Position, GRCh37.Alt.Allele)

    # If the variant position is found in df_database, look for the reference phenotype and allele name.
    if (nrow(matched_phenotypes) > 0) {

      reference <- df_database %>%
        filter(GRCh37.VCF.Position == "-" & 
              Gene == matched_phenotypes$Gene[1] & 
              chr_num == matched_phenotypes$Chromosome[1]) %>%
        dplyr::select(Phenotype, Allele.Name) %>%
        dplyr::slice_head(n = 1)

      ref_phenotypes[[i]] <- reference$Phenotype
      ref_allele_names[[i]] <- reference$Allele.Name

      # Retrieve variant positions associated to each retrieved phenotype.
      temp_coord_lst <- list()
      temp_alt_lst <- list()

      for (j in 1:nrow(matched_phenotypes)) {
        
        check_matched_phenotypes <- df_database %>%
            filter(Allele.Name == matched_phenotypes$Allele.Name[j]) %>%
            dplyr::select(Phenotype, Allele.Name, Gene, Chromosome, GRCh37.VCF.Position, GRCh37.Alt.Allele)

        # Check if all variant positions of the phenotype are in df_genotype. If so, store them in a list.
        if (all(check_matched_phenotypes$GRCh37.VCF.Position %in% variant_pos)) {
          alt_phenotypes[[i]] <- c(alt_phenotypes[[i]], check_matched_phenotypes$Phenotype[1])
          alt_allele_names[[i]] <- c(alt_allele_names[[i]], check_matched_phenotypes$Allele.Name[1])
          temp_coord_lst <- c(temp_coord_lst, list(check_matched_phenotypes$GRCh37.VCF.Position))
          temp_alt_lst <- c(temp_alt_lst, list(check_matched_phenotypes$GRCh37.Alt.Allele))
        } 
      }
      coord_combi_lst[[i]] <- temp_coord_lst
      alt_combi_lst[[i]] <- temp_alt_lst
    }

    # Fill up list with NA when no matches are found.
    alt_phenotypes[[i]] <- null_na_fill(alt_phenotypes[[i]])
    alt_allele_names[[i]] <- null_na_fill(alt_allele_names[[i]])
    ref_phenotypes[[i]] <- null_na_fill(ref_phenotypes[[i]])
    ref_allele_names[[i]] <- null_na_fill(ref_allele_names[[i]])
    coord_combi_lst[[i]] <- null_na_fill(coord_combi_lst[[i]])
    alt_combi_lst[[i]] <- null_na_fill(alt_combi_lst[[i]])
  }

  print(paste0("Assigning phenotype to chromosome ", chr_num))

  TEMP_SAVE_PATH <- paste0("./R/results/chr", chr_num, "_temp.rds")

  if (file.exists(TEMP_SAVE_PATH)) {
    df_inferred_result <- readRDS(TEMP_SAVE_PATH)
    last_processed_row <- which(df_inferred_result[, 1] != "") %>% max(na.rm = TRUE)
    print(paste0("Resuming from row: ", last_processed_row + 1, " for chromosome ", chr_num))
  }
  else {
    # Create an empty dataframe to store the result with the IDs as rows and variant positions as columns
    df_inferred_result <- data.frame(matrix(ncol = length(variant_pos), nrow = nrow(df_genotype)))
    colnames(df_inferred_result) <- variant_pos
    rownames(df_inferred_result) <- rownames(df_genotype)
    last_processed_row <- 0
  }

  # Iterate through the matrix to assign "phenotype/allele names" for each ID at each variant position.
  # Unable to run through all chr before R terminates.
  for (i in (last_processed_row + 1):nrow(mat_genotype)) {
    for (j in 1:ncol(mat_genotype)) {
      result <-character()
      coord_combi = coord_combi_lst[[j]]
      alt_combi = alt_combi_lst[[j]]
        for (k in seq_along(coord_combi)) {
            coord <- coord_combi[[k]]
            alt_ind <- alt_combi[[k]]
            if (!any(is.na(coord))) {
              if (all(as.numeric(mat_genotype[i, unlist(coord)]) == alt_ind)) {
                result <- c(result, paste0(alt_phenotypes[[j]][[k]], "/", alt_allele_names[[j]][[k]]))
              }
            }
        }
      if (length(result) > 0) {
        df_inferred_result[i, j] <- paste0(result, collapse = " | ")
      }
      else {
        df_inferred_result[i, j] <- paste0(ref_phenotypes[[j]], "/", ref_allele_names[[j]])
      }
    }
    if (i %% 100 == 0) {
      print(paste0("Processed row: ", i, " out of ", nrow(df_genotype)))
    }

    # Create a temporary file to save df_inferred_result every 500 rows in case R terminates
    if (i %% 500 == 0) {
      saveRDS(df_inferred_result, TEMP_SAVE_PATH)
      print(paste0("Row number ", i, " saved"))
    }
  }

  OUTPUT_PATH <- paste0("./R/results/chr", chr_num, "_inferred_phenotypes.tsv")
  write.table(df_inferred_result, file = OUTPUT_PATH, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  save.image(paste0("./R/results/chr", chr_num, "_inferred_phenotypes.Rdata"))
  print(paste0('Chromosome ', chr_num, ' Done'))
  file.remove(TEMP_SAVE_PATH)
}
