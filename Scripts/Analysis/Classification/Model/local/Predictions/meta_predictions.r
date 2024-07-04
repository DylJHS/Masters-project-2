# Create the predictor data for the meta learners
setwd("/Users/Dyll/Documents/Education/VU_UVA/Internship/Epigenetics/Janssen_Group-UMCUtrecht/Main_Project")

library(dplyr)
library(readr)
library(tidyverse)
library(xgboost)
library(caret)

input_path <- "Data/Gen_Model_input/"

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


meta_parameters <- read.csv(
  paste0(
    input_path,
    "Hyperparameters/meta_hyperparams.csv"
  )
)

# The CIN response features
# Categorical features
cat_cin <- read_tsv(
  paste0(
    input_path,
    "CIN_Features/PANCAN_ArmCallsAndAneuploidyScore_092817.txt"
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
print(head(cat_cin[20:30], 3))

# Numerical features
# HRD scores
ori_hrd <- read_tsv(
  paste0(input_path, "CIN_Features/TCGA.HRD_withSampleID.txt"),
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
    "/CIN_Features/lim_alpha_incl_TCGA_pericentro.csv"
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

# rm(peri_cnv)
# rm(hrd)

full_cin <- merge(
  cat_cin,
  reg_cin,
  by = "row.names"
) %>%
  column_to_rownames("Row.names")
cat("\n Full CIN:  \n")
print(head(full_cin[1:5], 3))

# All response feature names
response_features <- colnames(full_cin)
reg_features <- colnames(reg_cin)
cat_features <- colnames(cat_cin)
cat(
  "\n Response Features # :", length(response_features), "\n",
  "Regression Features # :", length(reg_features), "\n",
  "Categorical Features # :", length(cat_features), "\n"
)

base_oof_predictions <- read.csv(
  paste0(
    input_path,
    "Base_predictions/Full_base_predictions.csv"
  ),
  row.names = 1
)
cat("\n Base OOF Predictions:  \n")
print(head(base_oof_predictions[1:5], 3))

# Create the training folds
meta_folds <- createFolds(base_oof_predictions[["act_1p"]], k = 10, list = TRUE, returnTrain = FALSE)

# Create the empty dataframe that will store the out-of-fold predictions for each model of the meta learners
meta_oof_predictions <- data.frame(
  matrix(
    ncol = length(response_features),
    nrow = length(combine_all_folds(meta_folds))
  )
)
colnames(meta_oof_predictions) <- paste0("pred_", response_features)

cat("\n Predictions Dataframe:  \n")
print(head(meta_oof_predictions[1:5], 3))

meta_oof_predictions <- cbind(
  meta_oof_predictions,
  labels
) %>%
  arrange(act_index) %>%
  column_to_rownames("act_index")

meta_input_data <- base_oof_predictions %>%
  select(starts_with("pred_"))

for (feature in response_features) {
  cat(
    "Feature: ", feature, "\n\n"
  )

  # Get the parameters from stored Parameters file
  parameters <- meta_parameters

  # select the parameters and weights corresponding to the index
  selected_parameters <- parameters[parameters$Feature == feature, ]

  selected_feature <- selected_parameters$Feature
  selected_trees <- as.numeric(selected_parameters$Trees)
  selected_eta <- selected_parameters$Eta
  selected_gamma <- selected_parameters$Gamma
  selected_max_depth <- selected_parameters$Max_depth
  selected_min_child <- selected_parameters$Child_weight

  cat(
    "\n Selected Parameters: \n",
    "Feature: ", selected_feature, "\t",
    "Trees: ", selected_trees, "\t",
    "Eta: ", selected_eta, "\n",
    "Gamma: ", selected_gamma, "\t",
    "Max Depth: ", selected_max_depth, "\t",
    "Min Child: ", selected_min_child, "\n"
  )

  y <- as.integer(base_oof_predictions[[paste0("act_", feature)]])

  xgb_data <- xgb.DMatrix(data = as.matrix(meta_input_data), label = y)

  if (feature %in% cat_features) {
    xgb_meta_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "multi:softmax",
      eval_metric = "mlogloss",
      folds = meta_folds,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      early_stopping_rounds = early_stop,
      eta = selected_eta,
      gamma = selected_gamma,
      num_class = 3,
      print_every_n = print_every,
      prediction = TRUE
    )

    # Store the predictions in the corresponding column
    meta_oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_meta_model$pred[, 1]
  } else if (feature %in% reg_features) {
    xgb_meta_model <- xgb.cv(
      data = xgb_data,
      nrounds = selected_trees,
      objective = "reg:squarederror",
      eval_metric = "rmse",
      early_stopping_rounds = early_stop,
      folds = meta_folds,
      max_depth = selected_max_depth,
      min_child_weight = selected_min_child,
      eta = selected_eta,
      gamma = selected_gamma,
      print_every_n = print_every,
      prediction = TRUE
    )

    # Store the predictions in the corresponding column
    meta_oof_predictions[[paste0("pred_", selected_feature)]] <- xgb_meta_model$pred
  } else {
    cat("Feature not found")
  }
}

cat("Training of the meta models complete")

# Save the out-of-fold predictions
# write.csv(
#   meta_oof_predictions,
#   "Data/Gen_model_output/Predictions/Meta_predictions/Full_meta_predictions.csv"
# )
# write.csv(
#   meta_oof_predictions,
#   "Data/Gen_model_input/Meta_predictions/Full_meta_predictions.csv"
# )
