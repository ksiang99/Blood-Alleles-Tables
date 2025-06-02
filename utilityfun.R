# utilityfun.R
# utility functions used for prediction tasks, loading VCF files, loading LD files

# library: Use for interactive sessions
# require: Use for conditional loading, e.g., if (require(package_name)) { ... }
# attach: Avoid using, as it pollutes the global namespace.

# Load necessary libraries, do library(modules) first in main script as import() is a function from modules
#library(modules)
library(dplyr)
# library(caret)
# library(naivebayes)
# library(pROC)
# library(themis)

#import(c(predict,createDataPartition,confusionMatrix),from= "caret")
#import(naive_bayes, from="naivebayes")
#import(roc,from="pROC") # caret can make ROC AUC calculation
# library(VariantAnnotation)
# library(GenomicRanges)

# --- Function to load and preprocess LD data ---
# This function loads the LD results file, filters variants and creates a new column for the distance between variants
#
# @param ld_file_path (character) Path to the LD data file.
# @param PATH_TABLE (character) Path to the erythrogene allele position table
# @param MAF_A_threshold (numeric) Minimum minor allele frequency for SNP A (default: 0).
# @param MAF_B_threshold (numeric) Minimum minor allele frequency for SNP B (default: 0).
# @param R2_threshold (numeric) Minimum R-squared value for LD filtering (default: 0).
#
# @return (data.frame) A data frame containing the filtered and preprocessed LD data. 
#
load_and_preprocess_ld <- function(ld_file_path, PATH_TABLE = "../tables/process_tables_scripts/erythrogene_coordinate_fixed_with_chromosome.tsv", MAF_A_threshold = 0, MAF_B_threshold = 0, R2_threshold = 0) {
  df_allele <- read.delim(PATH_TABLE, header = TRUE)
  data <- read.table(ld_file_path, header = TRUE) %>% 
    mutate(DIST = abs(BP_A - BP_B)) %>%
    filter(MAF_A >= MAF_A_threshold, MAF_B >= MAF_B_threshold, SNP_A != SNP_B, R2 >= R2_threshold)    
  
  # create a new column GENE in data according to BP_A where GENE is from erythrogene coordinate table
  gene_map <- setNames(df_allele$Gene, df_allele$GRCh37_Start)
  data$GENE_A <- gene_map[as.character(data$BP_A)]
  data$GENE_B <- gene_map[as.character(data$BP_B)] # unlikely to be from list of blood type variants 
  
  # special case for XG SNP not in erythrogene table
  # if BP_A = 2666384 and CHR_A = 23, assign GENE_A as XG
  # data$GENE_A[is.na(data$GENE_A) & 
  #       data$BP_A == 2666384 & 
  #       data$CHR_A == 23] <- "XG"    
  return(data)
}

# --- Function to load and preprocess VCF data ---
#
# @param vcf_file_path (character) Path to the VCF file.
# @param chr_num (character) Chromosome number.
# @param positions (integer vector) Vector of variant positions to extract.
# @param widths (integer vector) Vector of variant widths to extract (default: 1).
# @return chr_vcf (VariantAnnotation VCF object) VCF object containing genotype data.
load_and_preprocess_vcf <- function(vcf_file_path, chr_num, positions, widths=1) {
  # copy chr_num to chr_num_
  chr_num_ <- chr_num

  # Convert chromosome number to string if necessary
  chr_num <- ifelse(chr_num == 23, "X", as.character(chr_num))
  
  # Create a genomic range object to extract specific positions from VCF file with custom widths
  gr <- GRanges(seqnames = chr_num, ranges = IRanges(start = positions, width=widths))
  
  # Read VCF file with specified genomic ranges
  chr_vcf <- readVcf(vcf_file_path, genome = "hg19", param = ScanVcfParam(which = gr))
  
  # Handle SNPs with multiple alleles
  chr_vcf <- expand(chr_vcf)
  return(chr_vcf)  
  # # Convert VCF object to data frame
  # if (chr_num_ == 23) {
  #   df_genotype <- VCF_to_df_X(chr_vcf)
  # } else {
  # df_genotype <- VCF_to_df(chr_vcf)
  # }
  # return(df_genotype)
  
}

# --- Function to convert VariantAnnotation VCF object into dataframe ---
#
# @param chr_vcf (VariantAnnotation VCF object) VCF object containing genotype data.
#
# @return df_genotype (data.frame) A data frame containing the phased genotype data with rows representing chromosomes and columns representing variant.
VCF_to_df <-function(chr_vcf){
    # convert genotype data to a dataframe of dimension number of variants x number of subjects
    genotype_data <- geno(chr_vcf)$GT # genotype_data is a matrix of characters, not a dataframe
    genotype_phased_df <- as.data.frame(genotype_data)
    genotype_phased_df <- genotype_phased_df[ ,order(names(genotype_phased_df))] # arrange columns by column header, in case subject IDs are not sorted

    # Split phased genotypes in the form of {0,1}|{0,1} into separate columns so that each column represent variants on the same chromosome
    # dataframe of dimension number of variants x number of chromosomes (or 2 x number of subjects, not accouting for gender and X chromosome)
    genotype_split_df <- as.data.frame(
        do.call(cbind, lapply(genotype_phased_df, function(geno_col) {
        do.call(rbind, strsplit(geno_col, split = "\\|"))
        }))
    )
    
    # use position of variant as row names
    variant_positions <- start(rowRanges(chr_vcf))

    # make.unique making up new positions, avoid problem with multiple variants at same position but creates columns with name containing "."
    unique_variant_positions <- make.unique(as.character(variant_positions))
    
    # row names are the variants
    rownames(genotype_split_df) <- unique_variant_positions

    # column names are the subject ID XXX and XXX.1 for 2 autosomes per subject
    colnames(genotype_split_df) <- make.unique(sort(c(colnames(chr_vcf), colnames(chr_vcf))))# column names are subject IDs
    
    # cast as data frame and make all entries factor type (cannot use integer or binary type for ML in R)
    df_genotype <- as.data.frame(t(genotype_split_df))
    df_genotype[] <- lapply(df_genotype, as.factor)

    # remove duplicated column with names containing ".", duplicates come from indels
    df_genotype <- df_genotype %>% dplyr::select(-contains("."))
    
    return(df_genotype)
}

# --- Function to convert VariantAnnotation VCF object into dataframe for X chromosome ---
VCF_to_df_X <- function(chr_vcf) {
  PATH_ <- "../igsr-1000_genomes_phase_3_release.tsv"
  
  # Read PATH_ file and fix naming
  gender_mapping <- read.csv(PATH_, sep="\t", stringsAsFactors=FALSE, row.names=NULL)
  names(gender_mapping)[1] <- "IID"
  names(gender_mapping)[names(gender_mapping) == "IID.Sex"] <- "Sex"

  # Extract subject IDs from the VCF
  subject_ids <- colnames(chr_vcf)

  # Match IDs with Sex
  male_ids <- gender_mapping$IID[gender_mapping$Sex == "male" & gender_mapping$IID %in% subject_ids]
  female_ids <- gender_mapping$IID[gender_mapping$Sex == "female" & gender_mapping$IID %in% subject_ids]

  # Convert genotype data to a data frame
  genotype_data <- geno(chr_vcf)$GT
  genotype_phased_df <- as.data.frame(genotype_data)
  genotype_phased_df <- genotype_phased_df[, order(names(genotype_phased_df))]

  # Split phased genotypes into separate columns
  genotype_split_df <- as.data.frame(
    do.call(cbind, lapply(genotype_phased_df, function(geno_col) {
      do.call(rbind, strsplit(geno_col, split = "\\|"))
    }))
  )

  # Use positions of variants as row names
  variant_positions <- start(rowRanges(chr_vcf))
  unique_variant_positions <- make.unique(as.character(variant_positions))
  rownames(genotype_split_df) <- unique_variant_positions

  # column names are the subject ID XXX and XXX.1 for 2 autosomes per subject
  colnames(genotype_split_df) <- make.unique(sort(c(colnames(chr_vcf), colnames(chr_vcf))))# column names are subject IDs

  # drop columns that are in male_ids
  genotype_split_df <- genotype_split_df[, !colnames(genotype_split_df) %in% male_ids]

  # Transpose and convert to a data frame for ML usage
  df_genotype <- as.data.frame(t(genotype_split_df))
  df_genotype[] <- lapply(df_genotype, as.factor)

  # Remove duplicate columns (e.g., indels creating duplicates)
  df_genotype <- df_genotype %>% dplyr::select(-contains("."))

  return(df_genotype)
}


# --- Function to Calculate and Plot ROC Curve ---
# @param model (object) A trained model compatible with caret.
# @param test_data (data.frame) Test data.
# @param target_col (character) Name of the target column.
# 
# @return (numeric) AUC value of the ROC curve, or NA if calculation fails.
calculate_and_plot_roc <- function(model, test_data, target_col) {
  # Identify predictor columns
  predictor_columns <- setdiff(colnames(test_data), target_col)
  
  # Ensure model provides probabilities
  test_probs <- tryCatch({
    predict(object = model, newdata = test_data[predictor_columns], type = "prob")
  }, error = function(e) {
    warning("Model prediction for probabilities failed: ", conditionMessage(e))
    return(NULL)
  })
  
  # Check if probabilities were successfully generated
  if (is.null(test_probs)) {
    warning("Prediction probabilities could not be generated. ROC calculation aborted.")
    return(NA)
  }
  
  # Check for presence of positive samples (class '0')
  if (sum(test_data[[target_col]] == 0) == 0) {
    warning("No positive samples (class '0') in the test data. ROC cannot be calculated.")
    return(NA)
  }
  
  # Calculate ROC and AUC
  auc <- tryCatch({
    roc_test <- roc(
      response = test_data[[target_col]],
      predictor = test_probs[, 1],  # Assuming the first column corresponds to class "0"
      levels = c(0, 1)  # Specify levels explicitly
    )
    
    # Plot ROC curve
    roc_data <- data.frame(
      specificity = roc_test$specificities,
      sensitivity = roc_test$sensitivities
    )
    
    ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
      geom_line(color = "blue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = paste("ROC Curve (AUC =", round(auc(roc_test), 3), ")"),
        x = "1 - Specificity (False Positive Rate)",
        y = "Sensitivity (True Positive Rate)"
      ) +
      theme_minimal()
    
    return(auc(roc_test))  # Return AUC value
    
  }, error = function(e) {
    warning("ROC calculation failed due to an error: ", conditionMessage(e))
    return(NA)  # Return NA if ROC calculation fails
  })
  
  return(auc)
}


# --- Function to subset df_genotype according to superpopulation ---
#export("split_by_superpopulation")
split_by_superpopulation <- function(df_genotype, df_population_label, superpopulation="ALL") {
  if (superpopulation == "ALL") {
    return(df_genotype)
  }

  # ID of subjects belonging to superpopulation
  population_ID <- df_population_label[df_population_label$Superpopulation.code == superpopulation,]

  # use ID to filter rows of df_genotype, keeping also rows with row name XXX.1 where XXX is ID to keep
  df_genotype_subset <- df_genotype[rownames(df_genotype) %in% population_ID$ID | 
                                  grepl("\\.1$", rownames(df_genotype)) & 
                                  gsub("\\.1$", "", rownames(df_genotype)) %in% population_ID$ID, ]
  return(df_genotype_subset)
}


# --- Helper Function to Calculate Performance Metrics ---
# @param predictions (factor) Model predictions.
# @param reference (factor) True labels.
# 
# @return (list) A list containing the ML performance metrics.
calculate_metrics <- function(predictions, reference) {
  conf <- confusionMatrix(data = predictions, reference = reference, mode = "everything", positive = "0")
  metrics <- list(
    Accuracy = conf$overall["Accuracy"],
    Kappa = conf$overall["Kappa"],
    Sensitivity = conf$byClass["Sensitivity"],
    Specificity = conf$byClass["Specificity"],
    Pos_Pred_Value = conf$byClass["Pos Pred Value"],
    Neg_Pred_Value = conf$byClass["Neg Pred Value"],
    Precision = conf$byClass["Precision"],
    Recall = conf$byClass["Recall"],
    F1 = conf$byClass["F1"],
    Prevalence = conf$byClass["Prevalence"],
    Detection_Rate = conf$byClass["Detection Rate"],
    Detection_Prevalence = conf$byClass["Detection Prevalence"],
    Balanced_Accuracy = (conf$byClass["Sensitivity"] + conf$byClass["Specificity"]) / 2
  )
  return(metrics)
}

# --- Function to Train and Evaluate a Machine Learning Model ---
# @param df_genotype (data.frame) Genotype data.
# @param target_column (character) Name of the target column.
# @param masked_columns (character vector) Names of columns to mask during training and testing, excluding the target column.
# @param model_type (character) Type of model to train ("naive_bayes", "logistic_regression", "random_forest").
# @param proportion_training_data (numeric) Proportion of data to use for training (default: 0.8).
# @param laplace_smoothing (numeric) Laplace smoothing parameter for Naive Bayes (default: 1).
# 
# @return (data.frame) A dataframe containing the full set of performance metrics for the model.
# 
# @export
train_and_evaluate_model <- function(df_genotype, target_column, masked_columns, 
                                      model_type = "naive_bayes", 
                                      proportion_training_data = 0.8) {
  # Validate input arguments
  if (!target_column %in% colnames(df_genotype)) {
    stop("Target column not found in the genotype data.")
  }
  
  if (!model_type %in% c("naive_bayes", "logistic_regression", "random_forest")) {
    stop("Unsupported model type. Supported types are 'naive_bayes', 'logistic_regression', and 'random_forest'.")
  }
  
  # Check for single record classes
  class_counts <- table(df_genotype[[target_column]])
  if (any(class_counts == 1)) {
    warning("Some classes have a single record. Skipping this target column: ", target_column)
    return(NULL)
  }
  
  # Remove target_column from masked_columns and mask columns from data
  masked_columns <- setdiff(masked_columns, target_column)
  df_genotype <- df_genotype[, !colnames(df_genotype) %in% masked_columns]

  # Split data into training and testing sets
  sub <- createDataPartition(y = df_genotype[[target_column]], p = proportion_training_data, list = FALSE)
  df_train <- df_genotype[sub, ]
  df_test <- df_genotype[-sub, ]

  # Train the specified model
  predictor_columns <- setdiff(colnames(df_train), target_column)
  model <- switch(model_type,
                "naive_bayes" = naive_bayes(x = df_train[predictor_columns], y = df_train[[target_column]], laplace = 1),
                "logistic_regression" = train(df_train[predictor_columns], df_train[[target_column]], method = "glm", family = "binomial"),
                "random_forest" = train(df_train[predictor_columns], df_train[[target_column]], method = "rf"))

  # Generate predictions and calculate performance metrics
  train_predictions <- predict(model, newdata = df_train[predictor_columns], type = "class")
  test_predictions <- predict(model, newdata = df_test[predictor_columns], type = "class")

  # Ensure model provides probabilities for ROC calculation
  test_probs <- tryCatch({
    predict(object = model, newdata = df_test[predictor_columns], type = "prob")
  }, error = function(e) {
    warning("Model prediction for probabilities failed: ", conditionMessage(e))
    return(NULL)
  })
  
  train_metrics <- calculate_metrics(train_predictions, df_train[[target_column]])
  test_metrics <- calculate_metrics(test_predictions, df_test[[target_column]])

  # AUC Calculation
  auc <- tryCatch({
    if (!is.null(test_probs)) {
      calculate_and_plot_roc(model, df_test, target_column)
    } else {
      NA
    }
  }, error = function(e) {
    warning("ROC calculation failed. Setting AUC to NA.")
    return(NA)
  })
  
  # Compile metrics into a dataframe
  metrics_df <- data.frame(
    Model = model_type,
    Target_Column = target_column,
    Proportion_Training_Data = proportion_training_data,
    MAF = 1 - sum(df_genotype[[target_column]] == 0) / length(df_genotype[[target_column]]), # avoiding problems with multiallelic positions
    # Laplace_Smoothing = ifelse(model_type == "naive_bayes", laplace_smoothing, NA),
    
    # Training Metrics
    Train_Accuracy = train_metrics$Accuracy,
    Train_Kappa = train_metrics$Kappa,
    Train_Sensitivity = train_metrics$Sensitivity,
    Train_Specificity = train_metrics$Specificity,
    Train_Pos_Pred_Value = train_metrics$Pos_Pred_Value,
    Train_Neg_Pred_Value = train_metrics$Neg_Pred_Value,
    Train_Precision = train_metrics$Precision,
    Train_Recall = train_metrics$Recall,
    Train_F1 = train_metrics$F1,
    Train_Prevalence = train_metrics$Prevalence,
    Train_Detection_Rate = train_metrics$Detection_Rate,
    Train_Detection_Prevalence = train_metrics$Detection_Prevalence,
    Train_Balanced_Accuracy = train_metrics$Balanced_Accuracy,
    
    # Testing Metrics
    Test_Accuracy = test_metrics$Accuracy,
    Test_Kappa = test_metrics$Kappa,
    Test_Sensitivity = test_metrics$Sensitivity,
    Test_Specificity = test_metrics$Specificity,
    Test_Pos_Pred_Value = test_metrics$Pos_Pred_Value,
    Test_Neg_Pred_Value = test_metrics$Neg_Pred_Value,
    Test_Precision = test_metrics$Precision,
    Test_Recall = test_metrics$Recall,
    Test_F1 = test_metrics$F1,
    Test_Prevalence = test_metrics$Prevalence,
    Test_Detection_Rate = test_metrics$Detection_Rate,
    Test_Detection_Prevalence = test_metrics$Detection_Prevalence,
    Test_Balanced_Accuracy = test_metrics$Balanced_Accuracy,
    
    # ROC AUC
    AUC = as.numeric(auc)
  )
  
  # Extract logistic regression details if model type is logistic regression
  if (model_type == "logistic_regression") {
    glm_model <- model$finalModel
    glm_summary <- summary(glm_model)
    glm_details <- data.frame(
      Target = target_column,
      Predictor = rownames(glm_summary$coefficients),
      Estimate = glm_summary$coefficients[, "Estimate"],
      Std_Error = glm_summary$coefficients[, "Std. Error"],
      Z_value = glm_summary$coefficients[, "z value"],
      P_value = glm_summary$coefficients[, "Pr(>|z|)"],
      AIC = glm_model$aic
    )
    return(list(Metrics = metrics_df, Logistic_Regression_Details = glm_details))
  }
  
  return(metrics_df)
}

# --- Function to Train and Evaluate a Machine Learning Model using XGBoost ---
# @param df_genotype (data.frame) Genotype data.
# @param target_column (character) Name of the target column.
# @param masked_columns (character vector) Names of columns to mask during training and testing, excluding the target column.

# @return (data.frame) A dataframe containing the full set of performance metrics for the model.
# 
# @export
train_and_evaluate_model_XG <- function(df_genotype, target_column, masked_columns) {
  # Validate input arguments
  if (!target_column %in% colnames(df_genotype)) {
    stop("Target column not found in the genotype data.")
  }
  
  # Check for single record classes
  class_counts <- table(df_genotype[[target_column]])
  if (any(class_counts == 1)) {
    warning("Some classes have a single record. Skipping this target column: ", target_column)
    return(NULL)
  }
  
  ## Ensure valid column names
  colnames(df_genotype) <- make.names(colnames(df_genotype))
  target_column <- make.names(target_column)
  masked_columns <- make.names(masked_columns)

  # Remove target_column from masked_columns and mask columns from data
  masked_columns <- setdiff(masked_columns, target_column)
  df_genotype_named <- df_genotype[, !colnames(df_genotype) %in% masked_columns]

  # Split data into training and testing sets
  sub <- createDataPartition(y = df_genotype_named[[target_column]], p = proportion_training_data, list = FALSE)
  df_train <- df_genotype_named[sub, ]
  df_test <- df_genotype_named[-sub, ]

  # Prepare the data for XGBoost
  predictor_columns <- setdiff(colnames(df_train), target_column)
  df_train[predictor_columns] <- lapply(df_train[predictor_columns], function(col) {
    as.numeric(as.character(col))
  })

  df_test[predictor_columns] <- lapply(df_test[predictor_columns], function(col) {
    as.numeric(as.character(col))
  })
 
  # Convert to DMatrix format for XGBoost
  dtrain <- xgb.DMatrix(data = as.matrix(df_train[, predictor_columns]), label = as.numeric(df_train[[target_column]]) - 1)
  dtest <- xgb.DMatrix(data = as.matrix(df_test[, predictor_columns]), label = as.numeric(df_test[[target_column]]) - 1)

  # Set hyperparameters for XGBoost
  params <- list(
    booster = "gbtree",            # Tree-based model
    objective = "multi:softmax",   # Multi-class classification
    num_class = length(unique(df_train[[target_column]])),  # Number of classes (genotypes)
    eta = 0.1,                     # Learning rate
    max_depth = 6,                 # Maximum depth of trees
    subsample = 0.8,               # Subsample ratio of training data
    colsample_bytree = 0.8         # Subsample ratio of features for each tree
  )

  # Train the model
  model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,  # Number of boosting rounds
    verbose = 1
  )

  # Generate predictions and calculate performance metrics
  train_predictions <- predict(model, dtrain)
  test_predictions <- predict(model, dtest)

  if ("2" %in% levels(train_predictions)) {
  train_predictions <- factor(train_predictions, levels = c("0", "1", "2"), labels = c("0", "1", "2"))
  } else {
  train_predictions <- factor(train_predictions, levels = c("0", "1"), labels = c("0", "1"))
  }

  if ("2" %in% levels(test_predictions)) {
  test_predictions <- factor(test_predictions, levels = c("0", "1", "2"), labels = c("0", "1", "2"))
  } else {
  test_predictions <- factor(test_predictions, levels = c("0", "1"), labels = c("0", "1"))
  }

  # Calculate performance metrics for training and testing sets
  train_metrics <- calculate_metrics(train_predictions, df_train[[target_column]])
  test_metrics <- calculate_metrics(test_predictions, df_test[[target_column]])

  # Compile metrics into a dataframe
  metrics_df <- data.frame(
    Model = model_type,
    Target_Column = target_column,
    Proportion_Training_Data = proportion_training_data,
    MAF = 1 - sum(df_genotype[[target_column]] == 0) / length(df_genotype[[target_column]]), # Minor Allele Frequency
    
    # Training Metrics
    Train_Accuracy = train_metrics$Accuracy,
    Train_Kappa = train_metrics$Kappa,
    Train_Sensitivity = train_metrics$Sensitivity,
    Train_Specificity = train_metrics$Specificity,
    Train_Pos_Pred_Value = train_metrics$Pos_Pred_Value,
    Train_Neg_Pred_Value = train_metrics$Neg_Pred_Value,
    Train_Precision = train_metrics$Precision,
    Train_Recall = train_metrics$Recall,
    Train_F1 = train_metrics$F1,
    Train_Prevalence = train_metrics$Prevalence,
    Train_Detection_Rate = train_metrics$Detection_Rate,
    Train_Detection_Prevalence = train_metrics$Detection_Prevalence,
    Train_Balanced_Accuracy = train_metrics$Balanced_Accuracy,
    Train_Pred_Class = paste(train_predictions, collapse = ", "),   
    Train_Pred_Row = paste(rownames(df_train), collapse = ", "),
    
    # Testing Metrics
    Test_Accuracy = test_metrics$Accuracy,
    Test_Kappa = test_metrics$Kappa,
    Test_Sensitivity = test_metrics$Sensitivity,
    Test_Specificity = test_metrics$Specificity,
    Test_Pos_Pred_Value = test_metrics$Pos_Pred_Value,
    Test_Neg_Pred_Value = test_metrics$Neg_Pred_Value,
    Test_Precision = test_metrics$Precision,
    Test_Recall = test_metrics$Recall,
    Test_F1 = test_metrics$F1,
    Test_Prevalence = test_metrics$Prevalence,
    Test_Detection_Rate = test_metrics$Detection_Rate,
    Test_Detection_Prevalence = test_metrics$Detection_Prevalence,
    Test_Balanced_Accuracy = test_metrics$Balanced_Accuracy,
    Test_Pred_Class = paste(test_predictions, collapse = ", "),
    Test_Pred_Row = paste(rownames(df_test), collapse = ", ")
  )

  return(metrics_df)
}

# --- Function to retrieve Phenotypes with positions found in df_genotype ---
# @param variant_pos (character vector) A vector of positions found in df_genotype
# @param chr_num (numeric) Chromosome number
# @param df_database (data.frame) A dataframe containing information on blood alleles

# @return (list) A list containing:
#   - ref_phenotypes: Reference phenotypes for each variant position.
#   - ref_allele_names: Reference allele names corresponding to ref_phenotypes.
#   - alt_phenotypes: Alternative phenotypes associated with each variant position.
#   - alt_allele_names: Alternative allele names corresponding to alt_phenotypes.
#   - coord_combi_lst: List of lists containing variant positions associated with alternative phenotypes.
#   - geno_combi_lst: List of lists containing genotype values of each variant positions of alternative alleles.
#
# 
# @export
get_phenotypes <- function(variant_pos, chr_num, df_database) {
  
  # Create vectors to store reference/alternative phenotype/allele names
  alt_phenotypes <- vector("list", length = length(variant_pos))
  alt_allele_names <- vector("list", length = length(variant_pos))
  ref_phenotypes <- vector("list", length = length(variant_pos))
  ref_allele_names <- vector("list", length = length(variant_pos))
  coord_combi_lst <- vector("list", length = length(variant_pos))
  geno_combi_lst <- vector("list", length = length(variant_pos))

  # Loop through each variant position
  for (i in seq_along(variant_pos)) {
    # Retrieve all phenotypes associated with the variant position
    matched_phenotypes <- df_database %>%
      filter(GRCh37.VCF.Position == variant_pos[i] &
             Chromosome == chr_num) %>%
      dplyr::select(Phenotype, Allele.Name, Gene, Chromosome, GRCh37.VCF.Position, GRCh37.Alt.Allele)

    # Look for the reference phenotype and allele name
    reference <- df_database %>%
      filter(GRCh37.VCF.Position == "-" & 
             Gene == matched_phenotypes$Gene[1] & 
             chr_num == matched_phenotypes$Chromosome[1]) %>%
      dplyr::select(Phenotype, Allele.Name) %>%
      dplyr::slice_head(n = 1)

    ref_phenotypes[[i]] <- reference$Phenotype
    ref_allele_names[[i]] <- reference$Allele.Name

    # Retrieve variant positions associated with each retrieved phenotype
    temp_coord_lst <- list()
    temp_geno_lst <- list()

    for (j in seq_len(nrow(matched_phenotypes))) {
      check_matched_phenotypes <- df_database %>%
        filter(Allele.Name == matched_phenotypes$Allele.Name[j]) %>%
        dplyr::select(Phenotype, Allele.Name, Gene, Chromosome, GRCh37.VCF.Position, GRCh37.Alt.Allele)

      # Check if all variant positions of the phenotype are in df_genotype
      if (all(check_matched_phenotypes$GRCh37.VCF.Position %in% variant_pos)) {
        alt_phenotypes[[i]] <- c(alt_phenotypes[[i]], check_matched_phenotypes$Phenotype[1])
        alt_allele_names[[i]] <- c(alt_allele_names[[i]], check_matched_phenotypes$Allele.Name[1])
        temp_coord_lst <- c(temp_coord_lst, list(check_matched_phenotypes$GRCh37.VCF.Position))
        temp_geno_lst <- c(temp_geno_lst, list(check_matched_phenotypes$GRCh37.Alt.Allele))
      } 
    }
    coord_combi_lst[[i]] <- temp_coord_lst
    geno_combi_lst[[i]] <- temp_geno_lst
  }

  return(list(
    ref_phenotypes = ref_phenotypes,
    ref_allele_names = ref_allele_names,
    alt_phenotypes = alt_phenotypes,
    alt_allele_names = alt_allele_names,
    coord_combi_lst = coord_combi_lst,
    geno_combi_lst = geno_combi_lst
  ))
}

# --- Function to assign phenotypes to each position in df_genotype ---
# @param df_genotype (data.frame) Genotype data.
# @param pheno_results (list) A list containing:
#   - ref_phenotypes: Reference phenotypes for each variant position.
#   - ref_allele_names: Reference allele names corresponding to ref_phenotypes.
#   - alt_phenotypes: Alternative phenotypes associated with each variant position.
#   - alt_allele_names: Alternative allele names corresponding to alt_phenotypes.
#   - coord_combi_lst: List of lists containing variant positions associated with alternative phenotypes.
#   - geno_combi_lst: List of lists containing genotype values of each variant positions of alternative alleles.
#
# @return (data.frame) A dataframe with the same dimensions as df_genotype, where each cell contains the inferred phenotypes for a sample at a given variant position
# 
# @export
pos_pheno_df <- function(df_genotype, pheno_results) {
  
  # Create an empty dataframe to store inferred phenotypes
  df_inferred_result <- data.frame(matrix(ncol = ncol(df_genotype), nrow = nrow(df_genotype)))
  colnames(df_inferred_result) <- variant_pos
  rownames(df_inferred_result) <- rownames(df_genotype)
  
  # Iterate through each row (sample/ID) in the genotype matrix
  for (i in 1:nrow(df_genotype)) {
    for (j in 1:ncol(df_genotype)) {
      result <- character()
      coord_combi <- pheno_results$coord_combi_lst[[j]]
      geno_combi <- pheno_results$geno_combi_lst[[j]]
    
      # Check for phenotype-allele matching
      for (k in seq_along(coord_combi)) {
        coord <- coord_combi[[k]]
        geno_ind <- geno_combi[[k]]

        if (!any(is.na(coord))) {
          if (all((df_genotype[i, unlist(coord)]) == geno_ind)) {
            result <- c(result, paste0(pheno_results$alt_phenotypes[[j]][[k]], "/", 
                                       pheno_results$alt_allele_names[[j]][[k]]))
          }
        }
      }

      # Assign phenotype/allele to the inferred result dataframe
      if (length(result) > 0) {
        df_inferred_result[i, j] <- paste0(result, collapse = " | ")
      } 
      else {
        df_inferred_result[i, j] <- paste0(pheno_results$ref_phenotypes[[j]], "/", 
                                           pheno_results$ref_allele_names[[j]])
      }
    }
    
    # Progress update every 100 rows
    if (i %% 100 == 0) {
      print(paste0("Processed row: ", i, " out of ", nrow(df_genotype)))
    }
  }
  
  return(df_inferred_result)
}

# --- Function to assign phenotypes respective blood groups for each sample ---
# @param df_sum_table (data.frame) Empty dataframe with rows as sample ID and columns blood group
# @param df_inferred_result (data.frame) Inferred phenotypes for samples at a given variant position
# @param df_genotype (data.frame) Genotype data
# @param bg_to_pheno_allele_map (list) map list of blood system to Pheno_Allele

# @return (data.frame) A filled up df_sum_table
# 
# @export
generate_sum_table <- function(df_sum_table, df_inferred_result, df_genotype, bg_to_pheno_allele_map) {
  
  # Create a vector to store unique phenotype for each sample
  unique_pheno_vec <- setNames(rep(NA_character_, nrow(df_genotype)), rownames(df_genotype))
  
  # Iterate through the samples to store their unique phenotypes in unique_pheno_vec
  for (i in seq_len(nrow(df_inferred_result))) {
    unique_pheno <- unique(as.vector(as.matrix(df_inferred_result[i, ])))
    unique_pheno <- unique(unlist(strsplit(unique_pheno, " \\| ")))
    unique_pheno_vec[i] <- paste(unique_pheno, collapse = " | ")
  }

  # Assign the unique phenotypes to their respective blood system for each sample
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
  return(df_sum_table)
}