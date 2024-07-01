library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)
library(caTools)

args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

cancer_types <- c("BLCA", "BRCA", "CESC", "HNSC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PRAD", "STAD", "THCA")

# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

# Load the constant data
# Arm Level Aneuploidies
chr_cnv <- read_tsv(
  "Data/Cancer_specific_data/Model_input/CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(
    Sample = str_replace_all(Sample, "-", "\\."),
    across(-Sample, ~ case_when(
      . == 1 ~ 2,
      . == 0 ~ 1,
      . == -1 ~ 0,
      is.na(.) ~ 1,
      TRUE ~ .
    ))
  ) %>%
  rename(
    "13q" = "13 (13q)",
    "14q" = "14 (14q)",
    "15q" = "15 (15q)",
    "21q" = "21 (21q)",
    "22q" = "22 (22q)"
  ) %>%
  column_to_rownames("Sample")

cat("\n\n CNV df: \n")
print(head(chr_cnv[, 1:5]))

# Define the features to be used
aneu_cat_feature_list <- colnames(chr_cnv)
selected_feature <- aneu_cat_feature_list[[index]]
cat(paste0("\n\n Selected feature: ", selected_feature, "\n"))
rm(aneu_cat_feature_list)

# Loop over the cancer types
for (cancer in cancer_types) {
  hyperparam_file <- paste0(
    "Data/Cancer_specific_data/Model_input/Parameters/Hyperparameters/",
    cancer,
    "/base_cat_hyperparams.csv"
  )
  rna_folder <- paste0(
    "Data/Cancer_specific_data/Model_input/RNA/",
    cancer,
    "/"
  )

  # Get the parameters from stored Parameters file
  parameters <- read.csv(hyperparam_file, header = TRUE)

  cat("\n Parameters: \n")
  print(dim(parameters))

  # select the parameters and weights corrsponding to the index
  selected_parameters <- parameters[index, ]
  rm(parameters)
  cat("\n Parameters: \n")
  print(selected_parameters)
  cat("\n\n")

  selected_feature <- selected_parameters$Feature
  selected_rna_set <- selected_parameters$RNA_set
  # selected_trees <- as.numeric(selected_parameters$Trees) + 500
  selected_trees <- 2500
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth
  selected_weights <- as.numeric(selected_parameters[c("Weight_loss", "Weight_normal", "Weight_gain")])

  rm(selected_parameters)

  # Define the extra parameters
  selected_seed <- 99

  # Get the corresonding rna set
  rna_list <- list(
    transcripts_per_million = "tpm",
    scaled_transcripts_per_million = "scld_tpm",
    log_scaled_transcripts_per_million = "log_scld_tpm",
    log_transcripts_per_million = "log_tpm",
    expected_counts = "exp",
    scaled_expected_counts = "scld_exp",
    log_expected_counts = "log_exp",
    log_scaled_expected_counts = "log_scld_exp"
  )

  rna_selection_name <- rna_list[[selected_rna_set]]

  rna_set <- read.csv(
    paste0(
      rna_folder,
      rna_selection_name,
      "_soi.csv"
    ),
    row.names = 1
  )
  cat("\n\n RNA df: \n")
  print(head(rna_set[, 1:5], 3))

  # MODELLING

  aneu_cat_metrics_df <- data.frame(
    RNA_set = character(),
    Trees = numeric(),
    Feature = character(),
    Depth = numeric(),
    Child_weight = numeric(),
    Eta = numeric(),
    Gamma = numeric(),
    Weight_loss = numeric(),
    Weight_norm = numeric(),
    Weight_gain = numeric(),
    Trained_mlogloss = numeric(),
    Test_mlogloss = numeric(),
    Seed = numeric()
  )

  # Determine the class weights for the target feature

  cat("\n\t selected weights: ")
  print(selected_weights)
  cat("\n")

  # Function to map factor levels to weights
  feature_digit_function <- function(factors) {
    sapply(factors, function(x) selected_weights[as.numeric(x)])
  }

  full_df <- merge(rna_set,
    chr_cnv,
    by = "row.names"
  )

  rm(rna_set)
  cat("\n\n full_df: \n")
  print(head(full_df[, 1:5], 3))
  cat("\n\n")

  y <- as.integer(full_df[[selected_feature]])
  cat("\n\n Target: ", selected_feature, "\n")
  print(head(y))

  X <- full_df %>% select(-c("Row.names", colnames(chr_cnv)))
  cat("\n\n Predictors: \n")

  # print(head(X[, 1:5]))

  train_y_factor <- factor(y, levels = c(0, 1, 2))
  weights <- as.numeric(feature_digit_function(train_y_factor))
  rm(train_y_factor)

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)
  rm(weights)
  rm(X)
  rm(y)

  grid <- expand.grid(
    min_child = seq(1, 9, 3),
    max_depth = seq(1, 9, 3),
  )

  for (j in 1:nrow(grid)) {
    for (param in names(grid)) {
      assign(paste0("selected_", param), grid[j, param])
    }

    set.seed(selected_seed)

    cat(paste0(
      "\t\t eta: ", selected_eta,
      "\t\t gamma: ", selected_gamma,
      "\t\t depth: ", selected_max_depth,
      "\t\t trees: ", selected_trees,
      "\t\t child_weight: ", selected_min_child,
      "\n"
    ))

    m_xgb_untuned <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "multi:softmax",
      eval_metric = "mlogloss",
      early_stopping_rounds = 250,
      nfold = 10,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      eta = selected_eta,
      gamma = selected_gamma,
      num_class = 3,
      stratified = TRUE,
      print_every_n = 5
    )

    best_iteration <- 0

    # First, check if best_iteration is valid
    if (is.null(
      m_xgb_untuned$best_iteration
    ) ||
      m_xgb_untuned$best_iteration < 1) {
      cat(paste0(
        "Warning: No valid best_iteration found.",
        " Using last iteration values instead.\n"
      ))
      # Use the last iteration if best_iteration is not valid
      best_iteration <- nrow(m_xgb_untuned$evaluation_log)
    } else {
      # Ensure that the best_iteration does not exceed the number of rows logged
      if (m_xgb_untuned$best_iteration > nrow(m_xgb_untuned$evaluation_log)) {
        cat(paste0(
          "Warning: best_iteration exceeds the number of rows in evaluation_log.",
          " Adjusting to maximum available.\n"
        ))
        best_iteration <- nrow(m_xgb_untuned$evaluation_log)
      } else {
        best_iteration <- m_xgb_untuned$best_iteration
      }


      best_mlogloss_train <- if (best_iteration > 0) {
        m_xgb_untuned$evaluation_log$train_mlogloss_mean[best_iteration]
      } else {
        NA # Or appropriate default/error value
      }

      best_mlogloss_test <- if (best_iteration > 0) {
        m_xgb_untuned$evaluation_log$test_mlogloss_mean[best_iteration]
      } else {
        NA # Or appropriate default/error value
      }

      cat(paste0(
        "The best iteration occurs with tree #: ",
        best_iteration, "\n\n"
      ))

      aneu_cat_metrics_df <- rbind(aneu_cat_metrics_df, data.frame(
        RNA_set = selected_rna_set,
        Trees = best_iteration,
        Feature = selected_feature,
        Depth = selected_max_depth,
        Child_weight = selected_min_child,
        Eta = selected_eta,
        Gamma = selected_gamma,
        Weight_loss = selected_weights[1],
        Weight_norm = selected_weights[2],
        Weight_gain = selected_weights[3],
        Trained_mlogloss = best_mlogloss_train,
        Test_mlogloss = best_mlogloss_test,
        Seed = selected_seed
      ))
    }
  }

  datetime <- Sys.time() %>%
    str_replace_all(" ", "_") %>%
    str_replace_all(":", "_") %>%
    str_replace_all("\\.", "_")

  name <- paste0(
    "Data/Cancer_specific_data/Model_output/Hyperparameters/",
    cancer, "/",
    selected_feature, "_",
    index, "_",
    datetime, ".csv"
  ) %>%
    str_replace_all(" ", "_") %>%
    str_replace_all(":", "_")

  write.csv(
    aneu_cat_metrics_df,
    file = name,
    row.names = FALSE
  )
}

cat(paste0("\n Completed processing for index: ", index, "\n"))
