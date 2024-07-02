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
## HRD scores
ori_hrd <- read_tsv("Data/Cancer_specific_data/Model_input/CIN/TCGA.HRD_withSampleID.txt",
  show_col_types = FALSE)

# Pericentromeric CNVs
peri_cnv <- read.csv("Data/Cancer_specific_data/Model_input/CIN/lim_alpha_incl_TCGA_pericentro.csv") %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID))

cat("\n\n pericentromeric data: \n")
print(head(peri_cnv[, 1:7], 3))

cat("\n\n All dfs loaded \n")

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\.")) %>%
  select(-HRD)

rownames(hrd) <- hrd$new
hrd <- hrd %>%
  select(-new)

cat("\n\n hrd data: \n")
print(head(hrd), 3)

rm(t_hrd)
rm(first_hrd)

full_cin <- merge(
  hrd,
  peri_cnv,
  by.x = "row.names",
  by.y = "sampleID"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")

cat("\n\n full cin data: \n")
print(head(full_cin[, 1:7], 3))
cat("\n\n Full cin features: \n")
print(colnames(full_cin))

rm(hrd)
rm(peri_cnv)

# Define the features to be used
aneu_reg_feature_list <- colnames(full_cin)
selected_feature <- aneu_reg_feature_list[[index]]
cat(paste0("\n\n Selected feature: ", selected_feature, "\n"))
rm(aneu_reg_feature_list)

# Loop over the cancer types
for (cancer in cancer_types) {
  cat(paste0("\n\n\t\t\t\t\t\t\t\t Processing for cancer: ", cancer, "\n"))

  hyperparam_file <- paste0(
    "Data/Cancer_specific_data/Model_input/Parameters/Hyperparameters/",
    cancer,
    "/base_reg_hyperparams.csv"
  )
  rna_folder <- paste0(
    "Data/Cancer_specific_data/Model_input/RNA/",
    cancer,
    "/"
  )

  # Get the parameters from stored Parameters file
  parameters <- read.csv(hyperparam_file, header = TRUE)

  # select the parameters and weights corresponding to the selected feature
  selected_parameters <- parameters[parameters$Feature == selected_feature,]
  rm(parameters)
  cat("\n Parameters: \n")
  print(selected_parameters)
  cat("\n\n")

  if (selected_feature != selected_parameters$Feature) {
    cat(
      "Warning: Selected feature does not match the feature in the parameters file.\n"
    )
    # stop the script for this feature
    stop("Stopping the process for this feature.")

  }

  selected_feature <- selected_parameters$Feature
  selected_rna_set <- selected_parameters$RNA_set
  # selected_trees <- as.numeric(selected_parameters$Trees) + 500
  selected_trees <- 2500
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth

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
  cat("\n RNA df: \n")
  print(head(rna_set[, 1:5], 3))

  # MODELLING

  aneu_reg_metrics_df <- data.frame(
    RNA_set = character(),
    Trees = numeric(),
    Feature = character(),
    Max_depth = numeric(),
    Child_weight = numeric(),
    Eta = numeric(),
    Gamma = numeric(),
    Trained_RMSE = numeric(),
    Test_RMSE = numeric()
  )

  full_df <- merge(rna_set,
    full_cin,
    by = "row.names"
  )

  rm(rna_set)
  cat("\n\n full_df: \n")
  print(head(full_df[, 1:5], 3))
  cat("\n\n")

  y <- as.integer(full_df[[selected_feature]])
  cat("\n\n Target: ", selected_feature, "\n")
  print(head(y))

  X <- full_df %>% select(-c("Row.names", colnames(full_cin)))
  cat("\n\n Predictors: \n")

  # print(head(X[, 1:5]))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)
  rm(X)
  rm(y)

  grid <- expand.grid(
    min_child = seq(1, 9, 3),
    max_depth = seq(1, 9, 3)
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
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = 10,
      nfold = 2,
      max_depth = selected_depth,
      min_child_weight = selected_min_child,
      eta = selected_eta,
      gamma = selected_gamma,
      print_every_n = 10
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
    }

    best_rmse_trained <- if (best_iteration > 0) {
      m_xgb_untuned$evaluation_log$train_rmse_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    best_rmse_test <- if (best_iteration > 0) {
      m_xgb_untuned$evaluation_log$test_rmse_mean[best_iteration]
    } else {
      NA # Or appropriate default/error value
    }

    cat(paste0(
      "\t\t The best iteration occurs with tree #: ",
      best_iteration, "\n\n"
    ))

    aneu_reg_metrics_df <- rbind(aneu_reg_metrics_df, data.frame(
      RNA_set = selected_rna_set,
      Trees = best_iteration,
      Feature = selected_feature,
      Max_depth = selected_depth,
      Child_weight = selected_min_child,
      Eta = selected_eta,
      Gamma = selected_gamma,
      Trained_RMSE = best_rmse_trained,
      Test_RMSE = best_rmse_test
    ))
  }

  saved_dir <- paste0(
    "Data/Cancer_specific_data/Model_output/Hyperparameters/",
    cancer, 
    "/"
  )

  if (!dir.exists(saved_dir)) {
    dir.create(saved_dir)
  }

  datetime <- Sys.time() %>%
    str_replace_all(" ", "_") %>%
    str_replace_all(":", "_") %>%
    str_replace_all("\\.", "_")

  name <- paste0(
    saved_dir,
    selected_feature, "_",
    index, "_",
    datetime, ".csv"
  ) %>%
    str_replace_all(" ", "_") %>%
    str_replace_all(":", "_")

  write.csv(
    aneu_reg_metrics_df,
    file = name,
    row.names = FALSE
  )
}

cat(paste0("\n\n\t\t\t\t Completed processing for index: ", index, "\n"))
