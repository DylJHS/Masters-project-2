# This script is the first part of the base model prediction pipeline.
# Its goal is to generate predictions for the all base models in parallel

# Load the packages
library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)

# Set the params
args <- commandArgs(trailingOnly = TRUE)
index <- as.numeric(args[1])

early_stopping <- 100
print_every <- 25
cv <- 5

input_path <- "/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project/Data/"

# Load the functions
# Function to map factor levels to weights
feature_digit_function <- function(factors) {
  sapply(factors, function(x) selected_weights[as.numeric(x)])
}

# Function to extract the fold indices from the folds
combine_all_folds <- function(folds) {
  all_indices <- c()
  for (fold in folds) {
    all_indices <- c(all_indices, fold)
  }
  return(all_indices)
}

# Set the constant variables
# The RNA list
rna_list <- list(
  transcripts_per_million = "tpm",
  scaled_transcripts_per_million = "scld_tpm",
  log_scaled_transcripts_per_million = "log_scld_tpm",
  log_transcripts_per_million = "log_tpm",
  expected_counts = "exp",
  scaled_expected_counts = "scld",
  log_expected_counts = "log_exp",
  log_scaled_expected_counts = "log_scld"
)

# print(rna_list)

# The CIN response features
# Categorical features
cat_cin <- read_tsv(
  paste0(
    input_path,
    "CIN/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
  ),
  show_col_types = FALSE
) %>%
  replace(is.na(.), 0) %>%
  select(-c("Type", "Aneuploidy Score")) %>%
  mutate(Sample = str_replace_all(Sample, "-", "\\.")) %>%
  column_to_rownames("Sample") %>%
  mutate_all(~ replace(., . == 1, 2)) %>%
  mutate_all(~ replace(., . == 0, 1)) %>%
  mutate_all(~ replace(., . == -1, 0)) %>%
  rename(
    "13q" = "13 (13q)",
    "14q" = "14 (14q)",
    "15q" = "15 (15q)",
    "21q" = "21 (21q)",
    "22q" = "22 (22q)"
  )

cat("\n Categorical CIN:  \n")
print(head(cat_cin[1:5]))

# Numerical features
# HRD scores
ori_hrd <- read_tsv(
  paste0(input_path, "CIN/TCGA.HRD_withSampleID.txt"),
  show_col_types = FALSE
)

t_hrd <- as.data.frame(t(ori_hrd))
first_hrd <- t_hrd
colnames(first_hrd) <- t_hrd[1, ]
hrd <- as.data.frame(first_hrd[-1, ]) %>%
  mutate_all(as.numeric) %>%
  rename(loh_hrd = "hrd-loh") %>%
  select(-HRD) %>%
  mutate(new = str_replace_all(rownames(.), "-", "\\."))

rm(first_hrd)
rm(ori_hrd)
rm(t_hrd)

rownames(hrd) <- hrd$new
hrd <- hrd %>%
  select(-new)

# print(head(hrd))

# Pericentromeric CNVs
peri_cnv <- read.csv(
  paste0(
    input_path,
    "CIN/lim_alpha_incl_TCGA_pericentro_cnv.csv"
  )
) %>%
  mutate_all(~ replace(., is.na(.), 0)) %>%
  mutate(sampleID = gsub("-", ".", sampleID)) %>%
  column_to_rownames("sampleID")

# print(dim(peri_cnv))

reg_cin <- merge(
  hrd,
  peri_cnv,
  by = "row.names"
) %>%
  mutate(Row.names = str_replace_all(Row.names, "-", ".")) %>%
  column_to_rownames("Row.names")
cat("\n Regression CIN:  \n")
print(head(reg_cin[1:5]))

rm(peri_cnv)
rm(hrd)

full_cin <- merge(
  cat_cin,
  reg_cin,
  by = "row.names"
) %>%
  column_to_rownames("Row.names")
cat("\n Full CIN:  \n")
print(head(full_cin[1:5]))

# All response feature names
response_features <- colnames(full_cin)
reg_features <- colnames(reg_cin)
cat_features <- colnames(cat_cin)
cat(
  "\n Response Features # :", length(response_features), "\n",
  "Regression Features # :", length(reg_features), "\n",
  "Categorical Features # :", length(cat_features), "\n"
)

rm(reg_cin)
rm(cat_cin)

# Create the folds to be used
full_data <- merge(
  read.csv(
    paste0(
      input_path,
      "mRNA/gbm_input/Train/train_tpm_soi.csv"
    ),
    row.names = 1
  ),
  full_cin,
  by = "row.names"
)
cat("\n Full Data:  \n")
print(head(full_data[1:5]))

# Save the sample ids with the index
ref_data <- full_data %>%
  mutate(act_index = as.factor(row_number())) %>%
  select(c("act_index", "Row.names"))

write.csv(
  ref_data,
  paste0(
    input_path,
    "model_output/base_predictions/XGB_pred_ref.csv"
  )
)
rm(ref_data)

# Creating the folds and returning the indices for the out-of-fold predictions only
folds <- createFolds(full_data[["1p"]], k = cv, list = TRUE, returnTrain = FALSE)
cat("\n Folds # :", length(folds), "\n")

# Create the empty dataframe that will store the out-of-fold predictions for each model
oof_predictions <- data.frame(
  matrix(
    ncol = length(response_features),
    nrow = length(combine_all_folds(folds))
  )
)
colnames(oof_predictions) <- paste0("pred_", response_features)

cat("\n Predictions Dataframe:  \n")
print(head(oof_predictions[1:5]))

# Get the actual fold labels
fold_labels <- full_data %>%
  select(all_of(intersect(response_features, colnames(.)))) %>%
  mutate(
    index = as.factor(row_number())
  ) %>%
  rename_with(~ paste0("act_", .))

# merge the fold labels with the predictions dataframe
oof_predictions <- cbind(oof_predictions, fold_labels) %>%
  arrange(act_index)


cat_parameters <- read.csv(
  paste0(
    input_path,
    "CIN/XGB_cat_parameters.csv"
  )
)
reg_parameters <- read.csv(
  paste0(
    input_path,
    "CIN/XGB_reg_parameters.csv"
  )
)


feature <- response_features[index]
cat(
  "\n index: ", index, "\n",
  "Feature: ", feature, "\n\n"
)

# Reduce the prediction dataframe to the selected feature
oof_predictions <- oof_predictions %>%
  select(c(
    "act_index",
    paste0("pred_", feature),
    paste0("act_", feature)
  ))

if (feature %in% cat_features) {
  # Get the parameters from stored Parameters file
  parameters <- cat_parameters

  # select the parameters and weights corresponding to the index
  selected_parameters <- parameters[parameters$Feature == feature, ]

  selected_feature <- selected_parameters$Feature
  selected_rna_set <- selected_parameters$RNA_set
  selected_trees <- as.numeric(selected_parameters$Trees)
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth
  selected_weights <- as.numeric(
    selected_parameters[c(
      "Weight_loss",
      "Weight_normal",
      "Weight_gain"
    )]
  )
  selected_min_child <- 1
  selected_seed <- 99

  cat(
    "\n Selected Parameters: \n",
    "Feature: ", selected_feature, "\n",
    "RNA Set: ", selected_rna_set, "\n",
    "Trees: ", selected_trees, "\n",
    "Eta: ", selected_eta, "\n",
    "Gamma: ", selected_gamma, "\n",
    "Max Depth: ", selected_max_depth, "\n",
    "Weights: ", selected_weights, "\n",
    "Min Child: ", selected_min_child, "\n",
    "Seed: ", selected_seed, "\n\n"
  )

  rna_selection_name <- rna_list[[selected_rna_set]]
  rna_set <- read.csv(
    paste0(
      input_path,
      "mRNA/gbm_input/Train/train_",
      rna_selection_name,
      "_soi.csv"
    ),
    row.names = 1
  )

  full_df <- merge(rna_set,
    full_cin,
    by = "row.names"
  )

  y <- as.integer(full_df[[selected_feature]])
  X <- full_df %>% select(-c("Row.names", all_of(response_features)))

  train_y_factor <- factor(y, levels = c(0, 1, 2))
  weights <- as.numeric(feature_digit_function(train_y_factor))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y, weight = weights)

  xgb_model <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "multi:softmax",
    eval_metric = "mlogloss",
    # early_stopping_rounds = early_stopping,
    folds = folds,
    max_depth = selected_max_depth,
    min_child_weight = selected_min_child,
    eta = selected_eta,
    gamma = selected_gamma,
    num_class = 3,
    print_every_n = print_every,
    prediction = TRUE
  )

  # Store the predictions in the corresponding column
  oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_model$pred[, 1]
} else if (feature %in% reg_features) {
  # Get the parameters from stored Parameters file
  parameters <- reg_parameters

  # select the parameters and weights corresponding to the index
  selected_parameters <- parameters[parameters$Feature == feature, ]

  selected_feature <- selected_parameters$Feature
  selected_rna_set <- selected_parameters$RNA_set
  selected_trees <- as.numeric(selected_parameters$Trees)
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth
  selected_min_child <- 1
  selected_seed <- 99

  cat(
    "\n Selected Parameters: \n",
    "Feature: ", selected_feature, "\n",
    "RNA Set: ", selected_rna_set, "\n",
    "Trees: ", selected_trees, "\n",
    "Eta: ", selected_eta, "\n",
    "Gamma: ", selected_gamma, "\n",
    "Max Depth: ", selected_max_depth, "\n",
    "Min Child: ", selected_min_child, "\n",
    "Seed: ", selected_seed, "\n\n"
  )

  rna_selection_name <- rna_list[[selected_rna_set]]
  rna_set <- read.csv(
    paste0(
      input_path,
      "mRNA/gbm_input/Train/train_",
      rna_selection_name,
      "_soi.csv"
    ),
    row.names = 1
  )

  full_df <- merge(rna_set,
    full_cin,
    by = "row.names"
  )

  y <- as.numeric(full_df[[selected_feature]])
  X <- full_df %>% select(-c("Row.names", all_of(response_features)))

  xgb_data <- xgb.DMatrix(data = as.matrix(X), label = y)

  xgb_model <- xgb.cv(
    data = xgb_data,
    nrounds = selected_trees,
    objective = "reg:squarederror",
    eval_metric = "rmse",
    # early_stopping_rounds = early_stopping,
    folds = folds,
    max_depth = selected_max_depth,
    min_child_weight = selected_min_child,
    eta = selected_eta,
    gamma = selected_gamma,
    print_every_n = print_every,
    prediction = TRUE
  )

  # Store the predictions in the corresponding column
  oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_model$pred
} else {
  cat("Feature not found")
}

write.csv(
  oof_predictions,
  paste0(
    input_path,
    "model_output/base_predictions/XGB_base_predictions_",
    feature,
    ".csv"
  )
)

cat("\n The script has finished running for feature: ", feature)
